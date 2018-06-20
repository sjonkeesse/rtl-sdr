#include <signal.h>
#include <string.h>
#include <stdio.h>      // Output to file
#include <stdlib.h>     //
#include <time.h>       // Time hours, minutes, seconds
#include <sys/time.h>   // Milliseconds
//#include <thread>
#include <unistd.h>
#include <math.h>       // Ceil
#include <pthread.h>
#include <libusb.h>

#include "rtl-sdr.h"
#include "convenience/convenience.h"


#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#define DEFAULT_BUFFER_LENGTH   (1 * 16384)
#define BUFFER_DUMP             (1<<12)

/**
 * Configuration of each tuning range
 */
struct tuning_state
{
    long int frequency_start;                // Frequency start in Hz
    long int frequency_end;                  // Frequency end in Hz
    long int frequency_center;               // Middle of the frequency
    int rate;                                // The bandwidth of the dongle in Hz
    uint8_t *buf8;                           //
    int bin_e;                               // Number 1 to 21. Needed to calculate buffer_length (1<<bin_e)
    int buffer_length;                       // The buffer size calculated with shifting
    long *avg;                               // Where samples are stored? (2^bin_e)
    int sample_count;                        // Amount of samples collected before emptying buffer?
    
    rtlsdr_dev_t *device;                    // Dedicated RTLSDR device, only here when tune hopping not supported
};

/**
  * Global variables
  */
int peak_hold = 1;                      // Use peak hold by default
int gain      = -1;                     // Gain 0 - 500, autogain = -1
int sample_rate;                        // The bandwidth of each tuning state
int total_tuning_states;                // Length of tuning_states
struct tuning_state tuning_states[50];  // Holds all the tuning_states
int device_count = 0;                   // Length of devices array
struct rtlsdr_dev_t *devices[20];       // Available SDR dongles, max 20       *?


int16_t* Sinewave;
double* power_table;
int N_WAVE, LOG2_N_WAVE;
int16_t *fft_buf;
int *window_coefs;

// TODO use struct for devices which also holds the current tuning state / position, and also a flag for if in use


void sine_table(int size)
{
    int i;
    double d;
    LOG2_N_WAVE = size;
    N_WAVE = 1 << LOG2_N_WAVE;
    Sinewave = malloc(sizeof(int16_t) * N_WAVE*3/4);
    power_table = malloc(sizeof(double) * N_WAVE);
    for (i=0; i<N_WAVE*3/4; i++)
    {
        d = (double)i * 2.0 * M_PI / N_WAVE;
        Sinewave[i] = (int)round(32767*sin(d));
        //printf("%i\n", Sinewave[i]);
    }
}

inline int16_t FIX_MPY(int16_t a, int16_t b)
{
    int c = ((int)a * (int)b) >> 14;
    b = c & 0x01;
    return (c >> 1) + b;
}

int fix_fft(int16_t iq[], int m)
{
    int mr, nn, i, j, l, k, istep, n, shift;
    int16_t qr, qi, tr, ti, wr, wi;
    n = 1 << m;
    if (n > N_WAVE)
    {return -1;}
    mr = 0;
    nn = n - 1;
    /* decimation in time - re-order data */
    for (m=1; m<=nn; ++m) {
        l = n;
        do
        {l >>= 1;}
        while (mr+l > nn);
        mr = (mr & (l-1)) + l;
        if (mr <= m)
        {continue;}
        // real = 2*m, imag = 2*m+1
        tr = iq[2*m];
        iq[2*m] = iq[2*mr];
        iq[2*mr] = tr;
        ti = iq[2*m+1];
        iq[2*m+1] = iq[2*mr+1];
        iq[2*mr+1] = ti;
    }
    l = 1;
    k = LOG2_N_WAVE-1;
    while (l < n) {
        shift = 1;
        istep = l << 1;
        for (m=0; m<l; ++m) {
            j = m << k;
            wr =  Sinewave[j+N_WAVE/4];
            wi = -Sinewave[j];
            if (shift) {
                wr >>= 1; wi >>= 1;}
            for (i=m; i<n; i+=istep) {
                j = i + l;
                tr = FIX_MPY(wr,iq[2*j]) - FIX_MPY(wi,iq[2*j+1]);
                ti = FIX_MPY(wr,iq[2*j+1]) + FIX_MPY(wi,iq[2*j]);
                qr = iq[2*i];
                qi = iq[2*i+1];
                if (shift) {
                    qr >>= 1; qi >>= 1;}
                iq[2*j] = qr - tr;
                iq[2*j+1] = qi - ti;
                iq[2*i] = qr + tr;
                iq[2*i+1] = qi + ti;
            }
        }
        --k;
        l = istep;
    }
    return 0;
}

double rectangle(int i, int length)
{
    return 1.0;
}

void fifth_order(int16_t *data, int length)
{
    int i;
    int a, b, c, d, e, f;
    a = data[0];
    b = data[2];
    c = data[4];
    d = data[6];
    e = data[8];
    f = data[10];
    /* a downsample should improve resolution, so don't fully shift */
    /* ease in instead of being stateful */
    data[0] = ((a+b)*10 + (c+d)*5 + d + f) >> 4;
    data[2] = ((b+c)*10 + (a+d)*5 + e + f) >> 4;
    data[4] = (a + (b+e)*5 + (c+d)*10 + f) >> 4;
    for (i=12; i<length; i+=4) {
        a = c;
        b = d;
        c = e;
        d = f;
        e = data[i-2];
        f = data[i];
        data[i/2] = (a + (b+e)*5 + (c+d)*10 + f) >> 4;
    }
}

void remove_dc(int16_t *data, int length)
{
    int i;
    int16_t ave;
    long sum = 0L;
    for (i=0; i < length; i+=2) {
        sum += data[i];
    }
    ave = (int16_t)(sum / (long)(length));
    if (ave == 0) {
        return;
    }
    for (i=0; i < length; i+=2) {
        data[i] -= ave;
    }
}

void downsample_iq(int16_t *data, int length)
{
    fifth_order(data, length);
    fifth_order(data+1, length-1);
}

long real_conj(int16_t real, int16_t imag)
{
    return ((long)real*(long)real + (long)imag*(long)imag);
}



/**
  * Create tuning_states to sample with a device not known yet
  */
void create_tuning_states(int frequency_start, int frequency_end, int dongle_bandwidth)
{
    fprintf(stderr, "Start frequency: %i Hz\nEnd frequency: %i Hz\n", frequency_start, frequency_end);
    
    // Frequency calculation
    struct tuning_state *tuning_state;
    int total_bandwidth     = frequency_end - frequency_start;
        total_tuning_states = ceil((double)total_bandwidth / (double)dongle_bandwidth);
        sample_rate         = total_bandwidth / total_tuning_states;
    
    // Print frequency info
    fprintf(stderr, "Total bandwidth: %i Hz\n", total_bandwidth);
    fprintf(stderr, "Dongle bandwidth: %i Hz\n", dongle_bandwidth);
    fprintf(stderr, "Total bandwidth / Dongle bandwidth = %i tuning states\n", total_tuning_states);
    fprintf(stderr, "Calculated sample rate: %i Hz\n", sample_rate);
    
    // Create buffer shizzle
    int bin_e, buffer_length, i, j;
    int downsample  = 1;
    int max_size    = 1000;
    double bin_size;
    
    // Calculate the minimum buffer size starting with 1250000 as bin_size
    for (bin_e = 1; bin_e <= 21; bin_e++) {
        // bin_size = Hz
        bin_size = (double)total_bandwidth / (double)((1<<bin_e) * downsample);
        //         2500000  /  (1, 2, 4, 8, 16, 32, 64)  *  1
        
        // Max size is downsample e.g. 1000 (1k)
        if (bin_size <= (double)max_size) {
            break;
        }
    }
    
    // Print bin_size info
    fprintf(stderr, "bin_e: %i, bin_size: %0.2fHz\n", bin_e, bin_size);
    
    buffer_length = 2 * (1<<bin_e) * downsample;
    
    if (buffer_length < DEFAULT_BUFFER_LENGTH) {
        buffer_length = DEFAULT_BUFFER_LENGTH;
    }
    
    // Print buffer_length info
    fprintf(stderr, "Buffer size: %i bytes (%0.2fms)\n", buffer_length, 1000 * 0.5 * (float)buffer_length / (float)total_bandwidth);
    
    for (i = 0; i < total_tuning_states; i++) {
        tuning_state                    = &tuning_states[i];
        tuning_state->frequency_start   = frequency_start + (i * sample_rate);
        tuning_state->frequency_end     = frequency_start + (i * sample_rate) + sample_rate;
        tuning_state->frequency_center  = frequency_start + (i * sample_rate) + (sample_rate / 2);
        tuning_state->rate              = sample_rate;
        tuning_state->bin_e             = bin_e;
        tuning_state->buffer_length     = buffer_length;
        tuning_state->buf8              = (uint8_t*)malloc(buffer_length * sizeof(uint8_t));
        tuning_state->sample_count      = 0;
        tuning_state->avg               = (long*)malloc((1<<bin_e) * sizeof(long));
//        tuning_state->gain              = gain; // this is set with initialisation
        
        // Checks
        if (!tuning_state->buf8) {
            fprintf(stderr, "Error: malloc tuning_state->buf8.\n");
            exit(1);
        }
        if (!tuning_state->avg) {
            fprintf(stderr, "Error: malloc.\n");
            exit(1);
        }
        
        for (j = 0; j < (1<<bin_e); j++) {
            tuning_state->avg[j] = 0L;
        }
    }
}

void initialize_devices()
{
    device_count = rtlsdr_get_device_count();
    static rtlsdr_dev_t *device;
    int iterations_count = ceil((double)total_tuning_states / (double)device_count);
    int response_code;
    
    // Reduce devices if possible when total iterations stays the same
    // Example: 7 tuning_states, 3 devices attached. For this are 3 iterations needed
    // Reducing the devices to 2, there are still 3 iterations needed.
    // Even known the devices are used in threads, personal experience tells me there will be
    // a more equal result when using an equal number of devices while sampling. The lesser, the better.
    while (device_count > 1 && iterations_count == ceil((double)total_tuning_states / (double)(device_count - 1))) {
        device_count--;
    }
    
    for (int device_index = 0; device_index < device_count; device_index++) {
        response_code = rtlsdr_open(&device, (uint32_t)device_index);
        
        if (response_code < 0) {
            fprintf(stderr, "Failed to open rtlsdr device #%d.\n", device_index);
            exit(1);
        }
        
        // Tuner gain (global)
        if (-1 == gain) {
            verbose_auto_gain(device);
        } else {
            verbose_gain_set(device, nearest_gain(device, nearest_gain(device, gain)));
        }
        
        // Reset buffer of device
        verbose_reset_buffer(device);
        
        // Sample rate (global)
        rtlsdr_set_sample_rate(device, sample_rate);
        
        // Add device to global devices array
        devices[device_index] = device;
    }
    
    fprintf(stderr, "%i device(s) initialized.\n", device_count);
}

void retune(rtlsdr_dev_t *device, int frequency)
{
    fprintf(stderr, "Retuning the RTLSDR device.\n");
    
    uint8_t dump[BUFFER_DUMP];
    int n_read;
    rtlsdr_set_center_freq(device, (uint32_t)frequency);
    /* wait for settling and flush buffer */
    usleep(5000);
    rtlsdr_read_sync(device, &dump, BUFFER_DUMP, &n_read);
    if (n_read != BUFFER_DUMP) {
        fprintf(stderr, "Error: bad retune.\n");
    }
}

static void *scan(void *arg)
{
    fprintf(stderr, "This is the scan fn speaking.\n");
    
    int tuning_state_index = arg;
    int n_read, i, offset, index;
    struct tuning_state *tuning_state = &tuning_states[tuning_state_index];
    int count = 0;
    int downsample = 1;
    int bin_length = 1 << tuning_state->bin_e;
    int32_t w;
    
    retune(tuning_state->device, tuning_state->frequency_center);
    
    while (1) {
        rtlsdr_read_sync(tuning_state->device, tuning_state->buf8, tuning_state->buffer_length, &n_read);
        
        // Loop trough buffer and put it in fft_buf array casted with int16_t
        for (i = 0; i < tuning_state->buffer_length; i++) {
            fft_buf[i] = (int16_t)tuning_state->buf8[i] - 127;
        }
        
        remove_dc(fft_buf, tuning_state->buffer_length / downsample);
        remove_dc(fft_buf + 1, (tuning_state->buffer_length / downsample) - 1);
        
        // Window function and FFT
        for (offset = 0; offset < (tuning_state->buffer_length / downsample); offset += (2 * bin_length)) {
            for (i = 0; i < bin_length; i++) {
                index = offset + (i * 2);
                w =  (int32_t)fft_buf[index];
                w *= (int32_t)(window_coefs[i]);
                fft_buf[index] = (int16_t)w;
                w =  (int32_t)fft_buf[index + 1];
                w *= (int32_t)(window_coefs[i]);
                fft_buf[index + 1] = (int16_t)w;
            }
        
            fix_fft(fft_buf + offset, tuning_state->bin_e);
        
            if (peak_hold) {
                for (i = 0; i < bin_length; i++) {
                    tuning_state->avg[i] = MAX(real_conj(fft_buf[offset + (i * 2)], fft_buf[offset + (i * 2) + 1]), tuning_state->avg[i]);
                }
            } else {
                for (i = 0; i < bin_length; i++) {
                    tuning_state->avg[i] += real_conj(fft_buf[offset + (i * 2)], fft_buf[offset + (i * 2) + 1]);
                }
            }
        
            tuning_state->sample_count += downsample;
        }
        
        // Output that shit
        int len, bin_count;
        long tmp;
        double dbm;
        len = 1 << tuning_state->bin_e;
        
        
        
        
        // Get Date, hours, minutes and seconds
        //strftime(stderr, 50, "%Y-%m-%d, %H:%M:%S, ", localtime(time(NULL)));
        if (count < 10) {
            fprintf(stderr, "2017-11-23, 23:57:0%i", count);
        } else {
            fprintf(stderr, "2017-11-23, 23:57:%i", count);
        }
        
        // Just info
        fprintf(stderr, ", %i, %i, %.2f, %i, ", tuning_state->frequency_start, tuning_state->frequency_end, (double)tuning_state->rate / (double)(len * downsample), tuning_state->sample_count * total_tuning_states);
        
        /* fix FFT stuff quirks */
        if (tuning_state->bin_e > 0) {
            /* nuke DC component (not effective for all windows) */
            tuning_state->avg[0] = tuning_state->avg[1];
            /* FFT is translated by 180 degrees */
            for (i = 0; i < (len / 2); i++) {
                tmp = tuning_state->avg[i];
                tuning_state->avg[i] = tuning_state->avg[i+len/2];
                tuning_state->avg[i + (len / 2)] = tmp;
            }
        }
        
        for (i=0; i <= (len - 1); i++) {
            dbm  = (double)tuning_state->avg[i];
            dbm /= (double)tuning_state->rate;
            
            if (0 == peak_hold) {
                dbm /= (double) tuning_state->sample_count;
            }
            
            dbm  = 10 * log10(dbm);
            fprintf(stderr, "%.2f, ", dbm);
        }
        
        // Clean up that puppy
        for (i=0; i<len; i++) {
            tuning_state->avg[i] = 0L;
        }
        
        tuning_state->sample_count = 0;
        
        // Print some info
//        fprintf(stderr, "%i\n", count);
//        fprintf(stderr, "buf8: %i, buffer_length: %i, n_read: %i\n", tuning_state->buf8, tuning_state->buffer_length, n_read);
        
        fprintf(stderr, "\n");
        
        count++;
        usleep(250000);
    }
}

//int parse_frequency(char *arg)
//{
//    char *frequency;
//    int upper, lower,;
//
//    start = arg;
//    stop = strchr(start, ':') + 1;
//    stop[-1] = '\0';
//    step = strchr(stop, ':') + 1;
//    step[-1] = '\0';
//    lower = (int)atofs(start);
//    upper = (int)atofs(stop);
//    stop[-1] = ':';
//    step[-1] = ':';
//    downsample = 1;
//    downsample_passes = 0;
//}

int main(int argc, char **argv)
{
    fprintf(stderr, "Main fn.\n");
    int interval            = 100;      // Milliseconds
    char *freqency_argument = "";       // 123000K:124M
    double dongle_bandwidth = 2600000;  // Hz per frequency hop
    struct timeval tp;
    long int epoch_time_ms, next_tick;
    int i, opt, length;
    double (*window_fn)(int, int) = rectangle; // 1.0

    while ((opt = getopt(argc, argv, "f:g:i:p")) != -1) {
        switch (opt) {
            case 'f': // start freqency:end frequency
                freqency_argument = strdup(optarg);
                break;
            case 'g': // Gain input is 0 to 50
                gain = (int)(atof(optarg) * 10);
                break;
            case 'i': // Interval in milliseconds (actually 'capture time')
                interval = (int)round(atoft(optarg));
                break;
            case 'p': // Use this flag to disable peak hold
                peak_hold = 0;
                break;
            case 'b': // This is the maximum bandwidth. It could be less.
                dongle_bandwidth = (double)(atof(optarg));
                break;
        }
    }
        
    if ("" == freqency_argument) {
        fprintf(stderr, "No frequency range provided.\n");
        exit(1);
    }
    
//    int frequency_start = 123
//    frequency_end   = parse_frequency(freqency_argument, 'start');
    
    // Create tuning_states freestanding from devices and store them in tuning_states array
//    create_tuning_states(380000000, 382500000, dongle_bandwidth);
    create_tuning_states(382500000, 385000000, dongle_bandwidth);
        
    initialize_devices();
    
    // Reserve some memory for buffer
    fft_buf = malloc(tuning_states[0].buffer_length * sizeof(int16_t));
    
    // Reserve some memory and add some precalculated shizzle
    length = 1 << tuning_states[0].bin_e;
    window_coefs = malloc(length * sizeof(int));
    
    sine_table(tuning_states[0].bin_e);
    
    for (i = 0; i < length; i++) {
        window_coefs[i] = (int)(256 * window_fn(i, length));
    }
    
    // Don't use frequency hopping from here...
    if (device_count != total_tuning_states) {
        fprintf(stderr, "%i device(s) found, %i device(s) needed (frequency hopping not supported yet).\n", device_count, total_tuning_states);
        exit(1);
    }
    
    // Create a thread for each tuning state / device
    struct tuning_state *tuning_state;
    pthread_t threads[device_count];
    for (i = 0; i < device_count; i++) {
        tuning_state = &tuning_states[i];
        tuning_state->device = devices[i]; // Add the device because we can only use one argument
        
//        fprintf(stderr, "Freq start %i, Freq end %i\n", tuning_state->frequency_start, tuning_state->frequency_end);
        
//        pthread_create(&threads[i], NULL, scan, (void *)(&tuning_state));
        pthread_create(&threads[i], NULL, scan, (void *)(i));
    }
    
    // Join all threads
    for (i = 0; i < device_count; i++) {
        pthread_join(threads[i], NULL);
    }
    
    
//    // Sample data
//    while (!do_exit) {
//        epoch_time_ms = (tp.tv_sec * 1000) + (tp.tv_usec / 1000);
//        next_tick = epoch_time_ms + interval;
//
//        while (epoch_time_ms < next_tick) {
//            // Create a job for each tuning state and add them to the thread pool
//            for (i = 0; i < total_tuning_states; i++) {
//                pool.AddJob(scan);
//            }
//
//            // Join all threads
//            // Don't do this, just use an extra thread to lock the bins of each tuner_state and concat and average all
//
//            usleep(interval * 1000);
//            gettimeofday(&tp, NULL);
//            epoch_time_ms = (tp.tv_sec * 1000) + (tp.tv_usec / 1000);
//        }
//
//        // Concat sample data
//        // ...
//
//        // Output data
//
//
//    }
    
    // On exit
    fprintf(stderr, "\nUser cancel, exiting...\n");
    
    // Free up reserved memory
    free(fft_buf);
    
    // Close all the devices
    for (i = 0; i < device_count; i++) {
        rtlsdr_close(devices[i]);
    }
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

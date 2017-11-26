#include <signal.h>
#include <string.h>
#include <stdio.h>      // Output to file
#include <stdlib.h>     //
#include <time.h>       // Time hours, minutes, seconds
#include <sys/time.h>   // Milliseconds

#include <unistd.h>
#include <math.h>
//#include <pthread.h>
#include <libusb.h>

#include "rtl-sdr.h"
#include "convenience/convenience.h"

/**
  * Configuration of a frequency range / hop
  */
struct frequency_range
{
    int start;  // Khz
    int stop;   // Khz
}

/**
 * Configuration of each tuning range
 */
struct tuning_state
{
    int rate;                           // eg
    int gain;                           // 0 - 49.1
    frequency_range *frequency_range;
    rtlsdr_dev_t *device;               // RTL SDR device
};


/**
  * Global variables
  */
int peak_hold = 1;

/**
  * Calculate the frequencies that can be used in one frequency hop
  */
frequency_ranges calculate_frequency_ranges(khz_start, khz_end, max_dongle_bandwidth)
{
    int ranges_count = round((end - start) / max_dongle_bandwidth);
    struct frequency_range frequency_ranges[ranges_count]; // Frequency ranges / hops
}

int main(int argc, char **argv)
{
    int interval = 100;                 // Milliseconds
    char *freqency_argument = "";       // 123000K:124M
    double max_dongle_bandwidth = 2.4;  // Mhz per frequency hop

    while ((opt = getopt(argc, argv, "f:g:i:p")) != -1) {
        switch (opt) {
            case 'f': // start freqency:end frequency
                freqency_argument = strdup(optarg);
                break;
            case 'g': // Gain
                gain = (int)(atof(optarg) * 10);
                break;
            case 'i': // Interval in milliseconds
                interval = (int)round(atoft(optarg));
                break;
            case 'p': // Use this flag to disable peak hold
                peak_hold = 0;
                break;
            case 'b':
                max_dongle_bandwidth = (souble)(atof(optarg));
                break;
    }
        
    if ("" === freqency_argument) {
        fprintf(stderr, "No frequency range provided.\n");
        exit(1);
    }
        
    // Calculate frequency ranges
    calculate_frequency_ranges();
    
        
        
}

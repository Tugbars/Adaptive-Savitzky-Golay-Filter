#ifndef ADAPTIVE_SAVGOL_H
#define ADAPTIVE_SAVGOL_H

#include <stdio.h>
#include <stdlib.h>

#define MAX_WINDOW 31 // Define the maximum window size
#define MIN_WINDOW 5  // Define the minimum window size

#define SMOOTHNESS_THRESHOLD 1.1
#define CORRELATION_THRESHOLD 0.99

// Adaptive Savitzky-Golay filter function
int adaptive_savgol_filter(MqsRawDataPoint_t* noisySignal, MqsRawDataPoint_t* smoothedSignal, int dataSize, int polyorder, double crit_val);
double calculate_smoothness(const MqsRawDataPoint_t* data, int n);
double calculate_correlation(const MqsRawDataPoint_t* x, const MqsRawDataPoint_t* y, int n);

// Function to print the data
void printData(const MqsRawDataPoint_t data[], size_t dataSize);

#endif // ADAPTIVE_SAVGOL_H

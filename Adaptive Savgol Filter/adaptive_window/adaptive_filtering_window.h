#ifndef ADAPTIVE_SAVGOL_H
#define ADAPTIVE_SAVGOL_H

#include "../mqs_def.h"
#include "../filter_state_machine/adaptive_filtering_config.h"

#define SMOOTHNESS_THRESHOLD 1.1
#define CORRELATION_THRESHOLD 0.99

// Adaptive Savitzky-Golay filter function
int adaptive_savgol_filter(MqsRawDataPoint_t* noisySignal, MqsRawDataPoint_t* smoothedSignal, int dataSize, int polyorder, double crit_val);
double calculate_smoothness(const MqsRawDataPoint_t* data, int n);
double calculate_correlation(const MqsRawDataPoint_t* x, const MqsRawDataPoint_t* y, int n);

#endif // ADAPTIVE_SAVGOL_H

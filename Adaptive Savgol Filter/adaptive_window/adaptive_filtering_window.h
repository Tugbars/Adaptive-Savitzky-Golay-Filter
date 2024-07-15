#ifndef ADAPTIVE_SAVGOL_H
#define ADAPTIVE_SAVGOL_H

#include "../mqs_def.h"
#include "../filter_state_machine/adaptive_filtering_config.h"

#define SMOOTHNESS_THRESHOLD 1.1
#define CORRELATION_THRESHOLD 0.99

// Adaptive Savitzky-Golay filter function
int adaptive_savgol_filter(MqsRawDataPoint_t* noisySignal, MqsRawDataPoint_t* smoothedSignal, uint16_t dataSize, uint8_t polyorder, double crit_val, uint16_t peakIndex, uint16_t interval_size);
double calculate_smoothness(const MqsRawDataPoint_t* data, uint16_t peakIndex, uint16_t interval_size);
double calculate_correlation(const MqsRawDataPoint_t* x, const MqsRawDataPoint_t* y, uint16_t peakIndex, uint16_t interval_size);
#endif // ADAPTIVE_SAVGOL_H

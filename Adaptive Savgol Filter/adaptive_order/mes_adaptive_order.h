#ifndef DEN_ORD_REG_H
#define DEN_ORD_REG_H

#include "../mqs_def.h"

typedef struct {
    int optimal_order;
    int optimal_window;
} OptimalOrderWindow;

typedef struct {
    int order;
    int window;
    double smoothness;
    double correlation;
} CacheEntry;

// Main denoising function using Optimum Order SG filter
double calculate_sigma(MqsRawDataPoint_t* noisy_sig, int len, char type);
void find_peak_range(MqsRawDataPoint_t* noisy_sig, int len, int interval_size, int* start, int* end);
void apply_optimal_filter_with_cache(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, int len, int start, int end, OptimalOrderWindow* results, double* best_smoothness, double* best_correlation, int* best_order, int* best_window);
void apply_optimal_filter(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, int len, int best_order, int best_window);
void evaluate_optimal_order_for_all_indexes(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, int start, int end, int len, double sigma, OptimalOrderWindow* optimal_order_windows, uint16_t peakIndex, int interval_size);
double find_primary_peak(MqsRawDataPoint_t* noisy_sig, int len, uint16_t* peakIndex);
#endif // DEN_ORD_REG_H

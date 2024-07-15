#ifndef DEN_ORD_REG_H
#define DEN_ORD_REG_H

#include "../mqs_def.h"

typedef struct {
    uint8_t optimal_order;
    uint8_t optimal_window;
} OptimalOrderWindow;

typedef struct {
    uint8_t order;
    uint8_t window;
    double smoothness;
    double correlation;
} CacheEntry;

// Main denoising function using Optimum Order SG filter
double find_primary_peak(MqsRawDataPoint_t* noisy_sig, uint16_t len, uint16_t* peakIndex);
double calculate_sigma(MqsRawDataPoint_t* noisy_sig, uint16_t len, char type);
void find_peak_range(MqsRawDataPoint_t* noisy_sig, uint16_t len, uint16_t interval_size, uint16_t* start, uint16_t* end);
void apply_optimal_filter_with_cache(
    MqsRawDataPoint_t* noisy_sig,
    MqsRawDataPoint_t* smoothed_sig,
    uint16_t len,
    uint16_t start,
    uint16_t end,
    OptimalOrderWindow* results,
    double* best_smoothness,
    double* best_correlation,
    uint8_t* best_order,
    uint8_t* best_window);

void evaluate_optimal_order_for_all_indexes(
    MqsRawDataPoint_t* noisy_sig,
    MqsRawDataPoint_t* smoothed_sig,
    uint16_t start,
    uint16_t end,
    uint16_t len,
    double sigma,
    OptimalOrderWindow* optimal_order_windows,
    uint16_t peakIndex,
    uint16_t interval_size
);

void apply_optimal_filter(
    MqsRawDataPoint_t* noisy_sig,
    MqsRawDataPoint_t* smoothed_sig,
    uint16_t len,
    uint8_t best_order,
    uint8_t best_window);


#endif // DEN_ORD_REG_H

#ifndef DEN_ORD_REG_H
#define DEN_ORD_REG_H

#include "../mqs_def.h"

#define RANGE_SIZE 5
#define ORDER_RANGE 5

typedef struct {
    int optimal_order;
    int optimal_window;
} OptimalOrderWindow;

#define MAX_RANGE_SIZE 5

typedef struct {
    int order;
    int window;
    double smoothness;
    double correlation;
} CacheEntry;

#define MAX_UNIQUE_COMBINATIONS 5

// Main denoising function using Optimum Order SG filter
double calculate_sigma(MqsRawDataPoint_t* noisy_sig, int len, char type);
void find_peak_range(MqsRawDataPoint_t* noisy_sig, int len, int* start, int* end);
void run_denoising_process(MqsRawDataPoint_t* rawData, size_t dataSize, char type, double lambda, double sigma, int M, int pmax, MqsRawDataPoint_t* denoised);
void apply_optimal_filter_with_cache(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, int len, int start, int end, OptimalOrderWindow* results, double* best_smoothness, double* best_correlation, int* best_order, int* best_window);
void apply_optimal_filter(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, int len, int best_order, int best_window);
void evaluate_optimal_order_for_all_indexes(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, int start, int end, int len, double sigma, double lambda, int pmin, int pmax, double GUE_MSE[RANGE_SIZE][ORDER_RANGE], OptimalOrderWindow* optimal_order_windows);


#endif // DEN_ORD_REG_H

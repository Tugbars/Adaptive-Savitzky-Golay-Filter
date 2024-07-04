#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "../savitzky_golay_filter/mes_savgol.h"
#include "mes_adaptive_order.h"
#include "../adaptive_window/adaptive_filtering_window.h"

/*
%%%% Input
     % M: half-window length
     % pmax: Maximum Order
     % noisy_sig: Noisy signal (1XN)
     % type: Noise type ('G': Gaussian, 'L': Laplacian, 'U': Uniform)
     % lamda: Regularization parameter
     % sigma: noise standard deviation
%%%% Output
    % denoised: Denoised output
    % Order: Optimum order at each time instants
    % GUE_MSE:  Regularized GUE-MSE corresponding to the order values pmin to pmax, at each time instants
    */

    // Function to estimate sigma

/**
 * @brief Calculate the noise standard deviation (sigma) for a given signal.
 *
 * This function estimates the noise standard deviation (sigma) for a given noisy signal
 * using different noise models: Gaussian, Laplacian, or Uniform.
 *
 * @param noisy_sig The noisy signal array.
 * @param len The length of the signal array.
 * @param type The type of noise ('G' for Gaussian, 'L' for Laplacian, 'U' for Uniform).
 * @return The estimated noise standard deviation.
 */
double calculate_sigma(MqsRawDataPoint_t* noisy_sig, int len, char type) {
    double sigma = 0;
    double median_diff = 0;

    for (int i = 1; i < len; i++) {
        median_diff += fabs(noisy_sig[i].phaseAngle - noisy_sig[i - 1].phaseAngle);
    }

    if (type == 'G') {
        sigma = (median_diff / (len - 1)) / (0.6745 * sqrt(2));
    }
    else if (type == 'L') {
        sigma = (median_diff / (len - 1)) / (1.1461 * sqrt(2));
    }
    else if (type == 'U') {
        sigma = (median_diff / (len - 1)) / ((2 - sqrt(2)) * sqrt(3));
    }

    return sigma;
}

/**
 * @brief Calculate the Generalized Unbiased Estimate of Mean Squared Error (GUE-MSE) for a given signal.
 *
 * This function computes the GUE-MSE, which is a regularized measure of the error between the noisy signal
 * and the smoothed signal, taking into account the noise standard deviation and a regularization parameter.
 *
 * @param noisy_sig The noisy signal array.
 * @param smoothed_sig The smoothed signal array.
 * @param len The length of the signal arrays.
 * @param sigma The noise standard deviation.
 * @param lambda The regularization parameter.
 * @param GUE_MSE The output array to store the GUE-MSE values.
 */
void calculate_gue_mse(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, int len, double sigma, double lambda, double* GUE_MSE) {
    for (int i = 0; i < len; i++) {
        double term1 = smoothed_sig[i].phaseAngle * smoothed_sig[i].phaseAngle;
        double term2 = -2 * smoothed_sig[i].phaseAngle * noisy_sig[i].phaseAngle;
        double term3 = 2 * sigma * sigma;
        double term4 = lambda * sigma * sigma;
        double term5 = noisy_sig[i].phaseAngle * noisy_sig[i].phaseAngle;
        GUE_MSE[i] = (term1 + term2 + term3 + term4 + term5) - sigma * sigma;
    }
}

/**
 * @brief Find the index of the maximum value in an array.
 *
 * This function searches for the maximum value in the 'phaseAngle' field of the given array and returns its index.
 *
 * @param a The input array of MqsRawDataPoint_t.
 * @param size The size of the input array.
 * @param col The column index (unused in this implementation).
 * @param max_val Pointer to store the maximum value found.
 * @param max_index Pointer to store the index of the maximum value found.
 * @return The index of the maximum value found in the array.
 */
static inline int maxrow(const MqsRawDataPoint_t a[], int size, int col, float* max_val, int* max_index) {
    for (int i = 0; i < size; i++) {
        if (*max_val < a[i].phaseAngle) {
            *max_val = a[i].phaseAngle;
            *max_index = i;
        }
    }
    return *max_index;
}

/**
 * @brief Recursively find the primary peak in a signal.
 *
 * This function uses a recursive approach to find the primary peak in the given signal array.
 *
 * @param a The input array of MqsRawDataPoint_t.
 * @param size The size of the input array.
 * @param l The left boundary of the current search window.
 * @param r The right boundary of the current search window.
 * @param peakIndex Pointer to store the index of the found peak.
 * @return The value of the peak found.
 */
static double findPeakRec(const MqsRawDataPoint_t a[], int size, int l, int r, uint16_t* peakIndex) {
    if (l > r) return -1;

    int mid = (l + r) / 2;
    float max_val = 0.0f;
    int max_index = 0;

    int max_row_index = maxrow(a, size, mid, &max_val, &max_index);

    if (mid == 0 || mid == size - 1) {
        *peakIndex = max_row_index;
        return max_val;
    }

    if (max_val < a[mid - 1].phaseAngle)
        return findPeakRec(a, size, l, mid - 1, peakIndex);
    else if (max_val < a[mid + 1].phaseAngle)
        return findPeakRec(a, size, mid + 1, r, peakIndex);
    else {
        *peakIndex = max_row_index;
        return max_val;
    }
}

/**
 * @brief Find the primary peak in a noisy signal.
 *
 * This function identifies the primary peak in the noisy signal array.
 *
 * @param noisy_sig The noisy signal array.
 * @param len The length of the signal array.
 * @param peakIndex Pointer to store the index of the found peak.
 * @return The value of the primary peak found.
 */
static double find_primary_peak(MqsRawDataPoint_t* noisy_sig, int len, uint16_t* peakIndex) {
    return findPeakRec(noisy_sig, len, 0, len - 1, peakIndex);
}

/**
 * @brief Find the range around the primary peak in a noisy signal.
 *
 * This function determines the start and end indices of a range around the primary peak in the noisy signal array.
 *
 * @param noisy_sig The noisy signal array.
 * @param len The length of the signal array.
 * @param start Pointer to store the start index of the range.
 * @param end Pointer to store the end index of the range.
 */
void find_peak_range(MqsRawDataPoint_t* noisy_sig, int len, int* start, int* end) {
    uint16_t peakIndex;
    double peakVal = find_primary_peak(noisy_sig, len, &peakIndex);

    *start = (peakIndex - 2 >= 0) ? peakIndex - 2 : 0;
    *end = (peakIndex + 2 < len) ? peakIndex + 2 : len - 1;
}

/**
 * @brief Calculate GUE-MSE for a specific order and update the GUE-MSE matrix.
 *
 * This function calculates the GUE-MSE for a specific polynomial order and updates the GUE-MSE matrix for the given range.
 *
 * @param noisy_sig The noisy signal array.
 * @param smoothed_sig The smoothed signal array.
 * @param start The start index of the range.
 * @param end The end index of the range.
 * @param len The length of the signal array.
 * @param sigma The noise standard deviation.
 * @param lambda The regularization parameter.
 * @param temp_gue_mse Temporary array to store GUE-MSE values for the current order.
 * @param GUE_MSE The output matrix to store GUE-MSE values for all orders.
 * @param p The polynomial order.
 * @param optimal_window Pointer to store the optimal window size for the current order.
 */
void calculate_gue_mse_for_order(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, int start, int end, int len, double sigma, double lambda, double* temp_gue_mse, double GUE_MSE[RANGE_SIZE][ORDER_RANGE], int p, int* optimal_window) {
    int window_size = adaptive_savgol_filter(noisy_sig, smoothed_sig, len, p, 0.99);
    calculate_gue_mse(&noisy_sig[start], &smoothed_sig[start], end - start + 1, sigma, lambda, temp_gue_mse);

    for (int i = start; i <= end; i++) {
        GUE_MSE[i - start][p - 3] = temp_gue_mse[i - start];
#ifdef DEBUG_PRINT
        printf("GUE_MSE[%d][%d] = %f (window size = %d)\n", i - start, p - 3, GUE_MSE[i - start][p - 3], window_size);
#endif
    }

    *optimal_window = window_size;
}

/**
 * @brief Evaluate the optimal polynomial order for each index in the given range.
 *
 * This function evaluates the optimal polynomial order for each index in the given range by calculating the GUE-MSE
 * for each order and selecting the order with the mini mum GUE-MSE.
 *
 * @param noisy_sig The noisy signal array.
 * @param smoothed_sig The smoothed signal array.
 * @param start The start index of the range.
 * @param end The end index of the range.
 * @param len The length of the signal array.
 * @param sigma The noise standard deviation.
 * @param lambda The regularization parameter.
 * @param pmin The minimum polynomial order to evaluate.
 * @param pmax The maximum polynomial order to evaluate.
 * @param GUE_MSE The output matrix to store GUE-MSE values for all orders.
 * @param optimal_order_windows The output array to store the optimal order and window size for each index.
 */
void evaluate_optimal_order_for_all_indexes(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, int start, int end, int len, double sigma, double lambda, int pmin, int pmax, double GUE_MSE[RANGE_SIZE][ORDER_RANGE], OptimalOrderWindow* optimal_order_windows) {
    for (int i = start; i <= end; i++) {
#ifdef DEBUG_PRINT
        printf("Testing index %d in the noisy dataset\n", i);
#endif
        double min_gue_mse = INFINITY;
        int best_order = pmin;
        int best_window = 0;

        for (int p = pmin; p <= pmax; p++) {
#ifdef DEBUG_PRINT
            printf("Testing order p = %d\n", p);
#endif
            int optimal_window;
            calculate_gue_mse_for_order(noisy_sig, smoothed_sig, start, end, len, sigma, lambda, GUE_MSE[i - start], GUE_MSE, p, &optimal_window);
            if (GUE_MSE[i - start][p - 3] < min_gue_mse) {
                min_gue_mse = GUE_MSE[i - start][p - 3];
                best_order = p;
                best_window = optimal_window;
            }
        }

        int real_index = i - start;
        optimal_order_windows[real_index].optimal_order = best_order;
        optimal_order_windows[real_index].optimal_window = best_window;
#ifdef DEBUG_PRINT
        // printf("real order/window Index i %d\n", real_index);
        printf("Index %d: Optimal Order = %d, Optimal Window = %d\n", i, best_order, best_window);
#endif
    }
}

/**
 * @brief Apply the optimal filter using cached results.
 *
 * This function applies the optimal filter using cached results of previously calculated smoothness and correlation values.
 *
 * @param noisy_sig The noisy signal array.
 * @param smoothed_sig The smoothed signal array.
 * @param len The length of the signal array.
 * @param start The start index of the range.
 * @param end The end index of the range.
 * @param results The array of optimal order and window size for each index.
 * @param best_smoothness Pointer to store the best smoothness value found.
 * @param best_correlation Pointer to store the best correlation value found.
 * @param best_order Pointer to store the best polynomial order found.
 * @param best_window Pointer to store the best window size found.
 */
void apply_optimal_filter_with_cache(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, int len, int start, int end, OptimalOrderWindow* results, double* best_smoothness, double* best_correlation, int* best_order, int* best_window) {
    static CacheEntry cache[MAX_UNIQUE_COMBINATIONS];
    int cache_count = 0;

    double best_smoothness_found = INFINITY;
    double best_correlation_found = 0.0;

    // Extract unique combinations
    int unique_count = 0;

    for (int i = start; i <= end; ++i) {
        bool found = false;
        for (int j = 0; j < unique_count; ++j) {
            if (cache[j].order == results[i - start].optimal_order && cache[j].window == results[i - start].optimal_window) {
                found = true;
                break;
            }
        }
        if (!found) {
            cache[unique_count].order = results[i - start].optimal_order;
            cache[unique_count].window = results[i - start].optimal_window;
            unique_count++;
        }
    }

    // Evaluate each unique combination
    for (int i = 0; i < unique_count; ++i) {
        int current_order = cache[i].order;
        int current_window = cache[i].window;
        int half_window_size = (current_window - 1) / 2;

        mes_savgolFilter(noisy_sig, len, half_window_size, smoothed_sig, current_order, 0, 0);
        cache[i].smoothness = calculate_smoothness(smoothed_sig, len);
        cache[i].correlation = calculate_correlation(noisy_sig, smoothed_sig, len);

        //printf("Order %d: Window = %d, smoothness = %f, correlation = %f\n", current_order, current_window, cache[i].smoothness, cache[i].correlation);
    }

    // Find the best combination based on smoothness and correlation criteria
    for (int i = 0; i < unique_count; ++i) {
        double smoothness = cache[i].smoothness;
        double correlation = cache[i].correlation;

        printf("Evaluating: Order %d, Window %d, Smoothness = %f, Correlation = %f\n", cache[i].order, cache[i].window, smoothness, correlation);

        if (correlation >= CORRELATION_THRESHOLD) {
            if (smoothness < best_smoothness_found) {
                best_smoothness_found = smoothness;
                best_correlation_found = correlation;
                *best_order = cache[i].order;
                *best_window = cache[i].window;
            }
        }
    }

    *best_smoothness = best_smoothness_found;
    *best_correlation = best_correlation_found;

    printf("Best Smoothness: %f for Order %d and Window %d\n", best_smoothness_found, *best_order, *best_window);
    printf("Best Correlation: %f for Order %d and Window %d\n", best_correlation_found, *best_order, *best_window);
}

/**
 * @brief Apply the optimal filter to the noisy signal.
 *
 * This function applies the optimal filter to the noisy signal using the specified polynomial order and window size.
 *
 * @param noisy_sig The noisy signal array.
 * @param smoothed_sig The smoothed signal array.
 * @param len The length of the signal array.
 * @param best_order The optimal polynomial order.
 * @param best_window The optimal window size.
 */
void apply_optimal_filter(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, int len, int best_order, int best_window) {
    int half_window_size = (best_window - 1) / 2;
    mes_savgolFilter(noisy_sig, len, half_window_size, smoothed_sig, best_order, 0, 0);
}

// denoise function using smoothed_sig directly
void den_ord_reg(int M, int pmax, MqsRawDataPoint_t* noisy_sig, int len, char type, double lambda, double sigma, MqsRawDataPoint_t* smoothed_sig, double GUE_MSE[RANGE_SIZE][ORDER_RANGE], OptimalOrderWindow* optimal_order_windows, double* best_smoothness, double* best_correlation) {
    int pmin = 3;

    if (sigma == -1) {
        sigma = calculate_sigma(noisy_sig, len, type);
    }

    int start, end;
    find_peak_range(noisy_sig, len, &start, &end);

#ifdef DEBUG_PRINT
    printf("Calculating GUE-MSE for each order and interval around the peak...\n");
#endif

    evaluate_optimal_order_for_all_indexes(noisy_sig, smoothed_sig, start, end, len, sigma, lambda, pmin, pmax, GUE_MSE, optimal_order_windows);
#ifdef DEBUG_PRINT
    printf("ARE WE HERE--------------------------------------------------");
#endif
    int best_order, best_window;
    apply_optimal_filter_with_cache(noisy_sig, smoothed_sig, len, start, end, optimal_order_windows, best_smoothness, best_correlation, &best_order, &best_window);

    // Apply the filter with the best order and window
    apply_optimal_filter(noisy_sig, smoothed_sig, len, best_order, best_window);
}

void run_denoising_process(MqsRawDataPoint_t* rawData, size_t dataSize, char type, double lambda, double sigma, int M, int pmax, MqsRawDataPoint_t* smoothed_sig) {
    int start, end;
    find_peak_range(rawData, dataSize, &start, &end);
    int rangeSize = end - start + 1;

    double GUE_MSE[RANGE_SIZE][ORDER_RANGE];
    double best_smoothness;
    double best_correlation;
    OptimalOrderWindow optimal_order_windows[MAX_RANGE_SIZE];

    den_ord_reg(M, pmax, rawData, dataSize, type, lambda, sigma, smoothed_sig, GUE_MSE, optimal_order_windows, &best_smoothness, &best_correlation);
}

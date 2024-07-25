#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "../savitzky_golay_filter/mes_savgol.h"
#include "../adaptive_window/adaptive_filtering_window.h"
#include "../filter_state_machine/adaptive_filtering_config.h"

#undef DEBUG_WINDOW
#define MAX_ITERATIONS ((MAX_WINDOW - MIN_WINDOW) / 2 + 1)

typedef struct {
    double correlation;
    double smoothness;
    uint8_t window_size;
} Result;

/*!
 * @brief Calculate the Pearson correlation coefficient for two MqsRawDataPoint_t arrays.
 *
 * This function calculates the Pearson correlation coefficient between two arrays
 * of MqsRawDataPoint_t, using their phase angles.
 *
 * @param x The first array of MqsRawDataPoint_t.
 * @param y The second array of MqsRawDataPoint_t.
 * @param n The number of elements in the arrays.
 * @return The Pearson correlation coefficient.
 */
#define SWEEP_SIZE 501 // ctx.len can be used instead. 
double calculate_correlation(const MqsRawDataPoint_t* x, const MqsRawDataPoint_t* y, uint16_t peakIndex, uint16_t interval_size) {
    int half_range = interval_size / 2;
    int start = (peakIndex - half_range >= 0) ? peakIndex - half_range : 0;
    int end = (start + interval_size - 1 < SWEEP_SIZE) ? start + interval_size - 1 : SWEEP_SIZE - 1;
    
    // Adjust the start if the end exceeds the array length
    if (end >= SWEEP_SIZE) {
        end = SWEEP_SIZE - 1;
        start = (end - interval_size + 1 >= 0) ? end - interval_size + 1 : 0;
    }

    double sum_x = 0, sum_y = 0;
    for (int i = start; i <= end; ++i) {
        sum_x += x[i].phaseAngle;
        sum_y += y[i].phaseAngle;
    }

    int n = end - start + 1;
    double mean1 = sum_x / n;
    double mean2 = sum_y / n;

    double sq_sum_x = 0, sq_sum_y = 0;
    double num = 0;

    for (int i = start; i <= end; ++i) {
        double diff_x = x[i].phaseAngle - mean1;
        double diff_y = y[i].phaseAngle - mean2;
        sq_sum_x += diff_x * diff_x;
        sq_sum_y += diff_y * diff_y;
        num += diff_x * diff_y;
    }

    double std1 = sqrt(sq_sum_x / n);
    double std2 = sqrt(sq_sum_y / n);

    double den = n * std1 * std2;
    return num / den;
}

/*!
 * @brief Check if a number is odd.
 *
 * This function returns true if the input number is odd, otherwise false.
 *
 * @param num The number to check.
 * @return True if the number is odd, otherwise false.
 */
static bool is_odd(int num) {
    return num % 2 != 0;
}

/*!
 * @brief Calculate the smoothness of the dataset.
 *
 * This function calculates the smoothness of the dataset using the sum of squared second derivatives
 * of the phase angles of MqsRawDataPoint_t elements.
 *
 * @param data The array of MqsRawDataPoint_t elements.
 * @param n The number of elements in the array.
 * @return The calculated smoothness of the dataset.
 */
double calculate_smoothness(const MqsRawDataPoint_t* data, uint16_t peakIndex, uint16_t interval_size) {
    int half_range = interval_size / 2;
    int start = (peakIndex - half_range >= 0) ? peakIndex - half_range : 0;
    int end = (start + interval_size - 1 < SWEEP_SIZE) ? start + interval_size - 1 : SWEEP_SIZE - 1;
    
    // Adjust the start if the end exceeds the array length
    if (end >= SWEEP_SIZE) {
        end = SWEEP_SIZE - 1;
        start = (end - interval_size + 1 >= 0) ? end - interval_size + 1 : 0;
    }

    double smoothness = 0;
    for (int i = start + 1; i < end - 1; ++i) {
        double second_derivative = data[i + 1].phaseAngle - 2 * data[i].phaseAngle + data[i - 1].phaseAngle;
        smoothness += second_derivative * second_derivative;
    }
    return smoothness;
}

/*!
 * @brief Apply the Savitzky-Golay filter and print the results.
 *
 * This function applies the Savitzky-Golay filter to a noisy signal and prints the smoothed signal.
 *
 * @param noisySignal The input array of noisy MqsRawDataPoint_t elements.
 * @param smoothedSignal The output array to store the smoothed MqsRawDataPoint_t elements.
 * @param dataSize The number of elements in the input array.
 * @param halfWindowSize The half window size for the Savitzky-Golay filter.
 * @param polyorder The polynomial order for the Savitzky-Golay filter.
 */
static inline void apply_savgol_filter(MqsRawDataPoint_t* noisySignal, MqsRawDataPoint_t* smoothedSignal, uint16_t dataSize, uint8_t halfWindowSize, uint8_t polyorder) {
    mes_savgolFilter(noisySignal, dataSize, halfWindowSize, smoothedSignal, polyorder, 0, 0);
}

/*!
 * @brief Select the best available window size based on correlation and smoothness criteria.
 *
 * This function iterates through the results array to find the best window size
 * that meets the specified correlation and smoothness thresholds. It updates the
 * best correlation, best smoothness, and best window size based on the criteria.
 *
 * @param results The array of Result structures containing correlation, smoothness, and window size.
 * @param iteration The number of iterations (or the number of valid elements in the results array).
 * @param crit_val The critical value for the correlation coefficient.
 * @param best_correlation Pointer to store the best correlation coefficient found.
 * @param best_smoothness Pointer to store the best smoothness value found.
 * @param best_window_size Pointer to store the best window size found.
 */
static void select_best_available_window(const Result* results, uint16_t iteration, double crit_val, double* best_correlation, double* best_smoothness, uint8_t* best_window_size) {
    *best_correlation = results[0].correlation;
    *best_smoothness = results[0].smoothness;
    *best_window_size = results[0].window_size;

    for (int i = 1; i < iteration; ++i) {
        if (results[i].correlation > *best_correlation && results[i].smoothness < SMOOTHNESS_THRESHOLD) {
            *best_correlation = results[i].correlation;
            *best_smoothness = results[i].smoothness;
            *best_window_size = results[i].window_size;
        }
        else if (results[i].correlation == *best_correlation && results[i].smoothness < *best_smoothness && results[i].smoothness < SMOOTHNESS_THRESHOLD) {
            *best_smoothness = results[i].smoothness;
            *best_window_size = results[i].window_size;
        }
    }
}

/*!
 * @brief Evaluate the results of the current window size and update the best values if criteria are met.
 *
 * This function checks if the current correlation and smoothness values meet the
 * specified thresholds. If the criteria are met, it updates the best correlation,
 * best smoothness, and best window size values.
 *
 * @param correlation The current correlation value.
 * @param smoothed_smoothness The current smoothness value.
 * @param crit_val The critical value for the correlation coefficient.
 * @param best_correlation Pointer to the best correlation value found so far.
 * @param best_smoothness Pointer to the best smoothness value found so far.
 * @param windowSize The current window size being evaluated.
 * @param best_window_size Pointer to the best window size found so far.
 * @return True if the current window size is better than the previous best, otherwise False.
 */
static bool evaluate_results_and_update_best(double correlation, double smoothed_smoothness, double crit_val, double* best_correlation, double* best_smoothness, uint8_t windowSize, uint8_t* best_window_size) {
    if (smoothed_smoothness < SMOOTHNESS_THRESHOLD) {
        if (correlation > *best_correlation && correlation > crit_val) {
            *best_correlation = correlation;
            *best_smoothness = smoothed_smoothness;
            *best_window_size = windowSize;
            return true;
        }
        else if (correlation == *best_correlation && smoothed_smoothness < *best_smoothness) {
            *best_smoothness = smoothed_smoothness;
            *best_window_size = windowSize;
            return true;
        }
    }
    return false;
}

/*!
 * @brief Calculate the smoothness and correlation of the smoothed signal.
 *
 * This function calculates the smoothness of the smoothed signal using the sum
 * of squared second derivatives of the phase angles. It also calculates the
 * Pearson correlation coefficient between the noisy and smoothed signals.
 *
 * @param noisySignal The input array of noisy MqsRawDataPoint_t elements.
 * @param smoothedSignal The output array of smoothed MqsRawDataPoint_t elements.
 * @param peakIndex The index of the primary peak in the signal.
 * @param interval_size The size of the interval around the peak for calculations.
 * @param smoothed_smoothness Pointer to store the calculated smoothness of the smoothed signal.
 * @param correlation Pointer to store the calculated Pearson correlation coefficient.
 */
static inline void calculate_smoothness_and_correlation(const MqsRawDataPoint_t* noisySignal, const MqsRawDataPoint_t* smoothedSignal, uint16_t peakIndex, uint16_t interval_size, double* smoothed_smoothness, double* correlation) {
    *correlation = calculate_correlation(noisySignal, smoothedSignal, peakIndex, interval_size);
    *smoothed_smoothness = calculate_smoothness(smoothedSignal, peakIndex, interval_size);
}

/*!
 * @brief Find the best window size for the Savitzky-Golay filter.
 *
 * This function finds the best window size for the Savitzky-Golay filter based on the correlation
 * and smoothness criteria. It iterates through possible window sizes and selects the best one.
 *
 * @param noisySignal The input array of noisy MqsRawDataPoint_t elements.
 * @param smoothedSignal The output array to store the smoothed MqsRawDataPoint_t elements.
 * @param dataSize The number of elements in the input array.
 * @param polyorder The polynomial order for the Savitzky-Golay filter.
 * @param crit_val The critical value for the correlation coefficient.
 * @param best_correlation A pointer to store the best correlation coefficient found.
 * @return The best window size for the Savitzky-Golay filter.
 */
static int find_best_window_size(MqsRawDataPoint_t* noisySignal, MqsRawDataPoint_t* smoothedSignal, uint16_t dataSize, uint8_t polyorder, double crit_val, double* best_correlation, uint16_t peakIndex, uint16_t interval_size) {
    int windowSize = g_adaptive_filtering_config.min_window;
    int best_window_size = g_adaptive_filtering_config.min_window;

    Result results[MAX_ITERATIONS];

    *best_correlation = -1;
    double best_smoothness = INFINITY;  // Initialize to maximum possible value
    int iteration = 0;

    double noisy_smoothness = calculate_smoothness(noisySignal, peakIndex, interval_size);

    while (windowSize <= g_adaptive_filtering_config.max_window && iteration < MAX_ITERATIONS) {
        int halfWindowSize = (windowSize - 1) / 2;
#ifdef DEBUG_WINDOW
        printf("Testing window size: %d halfWindowSize %d\n", windowSize, halfWindowSize);
#endif
        apply_savgol_filter(noisySignal, smoothedSignal, dataSize, halfWindowSize, polyorder);

        double correlation, smoothed_smoothness;
        calculate_smoothness_and_correlation(noisySignal, smoothedSignal, peakIndex, interval_size, &smoothed_smoothness, &correlation);

        results[iteration].correlation = correlation;
        results[iteration].smoothness = smoothed_smoothness;
        results[iteration].window_size = windowSize;

#ifdef DEBUG_WINDOW
        printf("Smoothed signal smoothness: %f signal correlation %f\n", smoothed_smoothness, correlation);
#endif

        if (smoothed_smoothness < SMOOTHNESS_THRESHOLD) {
            if (correlation > *best_correlation && correlation > crit_val) {
                *best_correlation = correlation;
                best_smoothness = smoothed_smoothness;
                best_window_size = windowSize;
            } else if (correlation == *best_correlation && smoothed_smoothness < best_smoothness) {
                best_smoothness = smoothed_smoothness;
                best_window_size = windowSize;
            }
        }

        // Stop if the correlation is decreasing and the smoothness is below the threshold
        if (iteration > 0 && correlation < results[iteration - 1].correlation && smoothed_smoothness < SMOOTHNESS_THRESHOLD) {
            break;
        }

        windowSize += 2;
        iteration++;
    }

    // If no suitable window size is found within the thresholds, select the best available one
    if (*best_correlation < crit_val) {
        select_best_available_window(results, iteration, crit_val, best_correlation, &best_smoothness, &best_window_size);
    }

    printf("Selected window size: %d\n", best_window_size);

    return best_window_size;
}

/*!
 * @brief Apply the adaptive Savitzky-Golay filter to the input signal.
 *
 * This function applies the adaptive Savitzky-Golay filter to the input signal. It finds the best
 * window size based on correlation and smoothness criteria, and then applies the filter using the
 * best window size.
 *
 * @param noisySignal The input array of noisy MqsRawDataPoint_t elements.
 * @param smoothedSignal The output array to store the smoothed MqsRawDataPoint_t elements.
 * @param dataSize The number of elements in the input array.
 * @param polyorder The polynomial order for the Savitzky-Golay filter.
 * @param crit_val The critical value for the correlation coefficient.
 */
int adaptive_savgol_filter(MqsRawDataPoint_t* noisySignal, MqsRawDataPoint_t* smoothedSignal, uint16_t dataSize, uint8_t polyorder, double crit_val, uint16_t peakIndex, uint16_t interval_size) {
    if (!noisySignal || !smoothedSignal || g_adaptive_filtering_config.min_window < 5 || !is_odd(g_adaptive_filtering_config.min_window) || dataSize <= 0) {
        //fprintf(stderr, "error\n");
        return g_adaptive_filtering_config.min_window;
    }

    double best_correlation;
    int best_window_size = find_best_window_size(noisySignal, smoothedSignal, dataSize, polyorder, crit_val, &best_correlation, peakIndex, interval_size);

    int best_half_window_size = (best_window_size - 1) / 2;
    mes_savgolFilter(noisySignal, dataSize, best_half_window_size, smoothedSignal, polyorder, 0, 0);

    return best_window_size;
}

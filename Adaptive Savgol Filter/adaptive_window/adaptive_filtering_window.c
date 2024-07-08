#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "../savitzky_golay_filter/mes_savgol.h"
#include "../adaptive_window/adaptive_filtering_window.h"
#include "../filter_state_machine/adaptive_filtering_config.h"


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
double calculate_correlation(const MqsRawDataPoint_t* x, const MqsRawDataPoint_t* y, int n) {
    double sum_x = 0, sum_y = 0;
    for (int i = 0; i < n; ++i) {
        sum_x += x[i].phaseAngle;
        sum_y += y[i].phaseAngle;
    }

    double mean1 = sum_x / n;
    double mean2 = sum_y / n;

    double sq_sum_x = 0, sq_sum_y = 0;
    double num = 0;

    for (int i = 0; i < n; ++i) {
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
double calculate_smoothness(const MqsRawDataPoint_t* data, int n) {
    double smoothness = 0;
    for (int i = 1; i < n - 1; ++i) {
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
static void apply_savgol_filter(MqsRawDataPoint_t* noisySignal, MqsRawDataPoint_t* smoothedSignal, int dataSize, int halfWindowSize, int polyorder) {
    mes_savgolFilter(noisySignal, dataSize, halfWindowSize, smoothedSignal, polyorder, 0, 0);
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

#define MAX_ITERATIONS ((MAX_WINDOW - MIN_WINDOW) / 2 + 1)

static int find_best_window_size(MqsRawDataPoint_t* noisySignal, MqsRawDataPoint_t* smoothedSignal, int dataSize, int polyorder, double crit_val, double* best_correlation) {
    int windowSize = g_adaptive_filtering_config.min_window;
    int best_window_size = g_adaptive_filtering_config.min_window;

    double correlations[MAX_ITERATIONS] = { 0 };
    int window_sizes[MAX_ITERATIONS] = { 0 };

    *best_correlation = -1;
    int iteration = 0;

    double noisy_smoothness = calculate_smoothness(noisySignal, dataSize);
#ifdef DEBUG_PRINT
    printf("Noisy signal smoothness: %f\n", noisy_smoothness);
#endif

    while (windowSize <= g_adaptive_filtering_config.max_window && iteration < MAX_ITERATIONS) {
        int halfWindowSize = (windowSize - 1) / 2;

#ifdef DEBUG_PRINT
        printf("Testing window size: %d halfWindowSize %d\n", windowSize, halfWindowSize);
#endif
        apply_savgol_filter(noisySignal, smoothedSignal, dataSize, halfWindowSize, polyorder);

        double correlation = calculate_correlation(noisySignal, smoothedSignal, dataSize);
#ifdef DEBUG_PRINT
        printf("correlation %f\n", correlation);
#endif
        double smoothed_smoothness = calculate_smoothness(smoothedSignal, dataSize);
#ifdef DEBUG_PRINT
        printf("Smoothed signal smoothness: %f\n", smoothed_smoothness);
#endif
        correlations[iteration] = correlation;
        window_sizes[iteration] = windowSize;

        if (correlation > *best_correlation && fabs(correlation) > crit_val && smoothed_smoothness < SMOOTHNESS_THRESHOLD) {
            *best_correlation = correlation;
            best_window_size = windowSize;
        }

        windowSize += 2;
        iteration++;
    }

    // If no suitable window size is found within the thresholds, select the best available one
    if (*best_correlation < crit_val) {
        *best_correlation = correlations[0];
        best_window_size = window_sizes[0];

        for (int i = 1; i < iteration; ++i) {
            if (correlations[i] > *best_correlation) {
                *best_correlation = correlations[i];
                best_window_size = window_sizes[i];
            }
        }
    }

#ifdef DEBUG_PRINT
    printf("Selected window size: %d\n", best_window_size);
#endif

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
int adaptive_savgol_filter(MqsRawDataPoint_t* noisySignal, MqsRawDataPoint_t* smoothedSignal, int dataSize, int polyorder, double crit_val) {
    if (!noisySignal || !smoothedSignal || g_adaptive_filtering_config.min_window < 5 || !is_odd(g_adaptive_filtering_config.min_window) || dataSize <= 0) {
        fprintf(stderr, "error\n");
        return g_adaptive_filtering_config.min_window;
    }

    double best_correlation;
    int best_window_size = find_best_window_size(noisySignal, smoothedSignal, dataSize, polyorder, crit_val, &best_correlation);

#ifdef DEBUG_PRINT
    printf("Selected window size: %d\n", best_window_size);
#endif
    int best_half_window_size = (best_window_size - 1) / 2;
    mes_savgolFilter(noisySignal, dataSize, best_half_window_size, smoothedSignal, polyorder, 0, 0);

    return best_window_size;
}

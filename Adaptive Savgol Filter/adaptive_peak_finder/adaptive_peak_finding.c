/*!
 * Adaptive Peak Finding and Analysis Implementation.
 *
 * Author: Tugbars Heptaskin
 * Date: 06/18/2024
 * Company: Aminic Aps
 *
 * This implementation provides algorithms for detecting and analyzing peaks in noisy signals.
 * The approach is based on adaptive peak detection, which dynamically adjusts detection thresholds
 * and criteria based on local signal properties, such as mean and standard deviation.
 *
 * The algorithms in this implementation utilize local mean and local standard deviation to find peaks that adhere
 * to the overall trend of the dataset in a given window. The implementation considers both the full width at half
 * maximum (FWHM) and the width at 10% of the peak height, ensuring robust peak detection and analysis.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "adaptive_peak_finding.h"


 /*!
  * @brief Determines if a peak is still climbing at the end of a dataset.
  *
  * This function assesses whether the identified peak in a dataset is still rising
  * as it reaches the end of the dataset. This is important in peak finding algorithms,
  * particularly when analyzing segments of data where the peak might extend beyond
  * the current dataset's boundary.
  *
  * The function iterates from the peak index to the end of the dataset, calculating
  * the derivative (rate of change) at each point. It checks if this derivative is
  * less than or equal to a specified noise tolerance. If the condition fails more than once,
  * it indicates that the peak is no longer climbing.
  *
  * This check helps to determine if the current dataset's peak is part of a larger peak
  * that might be fully realized in subsequent datasets. If the peak is still climbing
  * at the end of the current dataset, there may be a need to analyze the next dataset
  * to find the true peak.
  *
  * @param b The array of data points (MqsRawDataPoint_t) containing the peak.
  * @param sizeB The size of the array.
  * @param peakIndex The index of the peak within the array.
  * @param noiseTolerance The tolerance level for the derivative to be considered noise.
  * @return True if the peak is still climbing; false otherwise.
  */
static bool isPeakClimbing(MqsRawDataPoint_t b[], int sizeB, int peakIndex, float noiseTolerance)
{
    if (peakIndex <= 0 || peakIndex >= sizeB - 1)
    {
        return false;
    }

    int failCount = 0; // Counter for the number of times condition is not met

    for (int i = peakIndex; i < sizeB - 1; i++)
    {
        float derivativeAfter = b[i + 1].phaseAngle - b[i].phaseAngle;

        // Check if the derivative after is less than or equal to the noise tolerance
        if (derivativeAfter <= noiseTolerance)
        {
            failCount++;
            if (failCount >= 2) // Check if it's the second time
            {
                return false; // Peak is not climbing if condition failed twice
            }
        }
    }

    // Return true only if failCount is less than 2
    return failCount < 2;
}

/*!
 * @brief Finds the index of the maximum value in a column of a 1D array.
 *
 * @param a The array of data points (MqsRawDataPoint_t) to search through.
 * @param size The number of elements in the array.
 * @param col The column in the array to search for the maximum value.
 * @param max_val A pointer to store the maximum value found.
 * @param max_index A pointer to store the index of the maximum value.
 * @return The index of the maximum value found in the specified column.
 */
static inline int maxrow(const MqsRawDataPoint_t a[], int size, int col, float* max_val, int* max_index)
{
    for (int i = 0; i < size; i++)
    {
        if (*max_val < a[i].phaseAngle)
        {
            *max_val = a[i].phaseAngle;
            *max_index = i;
        }
    }
    return *max_index;
}


/*!
 * @brief Recursively finds a peak in a dataset using a divide-and-conquer approach.
 *
 * This function implements a recursive peak finding algorithm. It divides the dataset
 * into two halves at each recursive step and determines the direction (left or right)
 * to continue the search based on the comparison of adjacent values. This divide-and-conquer
 * approach significantly reduces the time complexity compared to a linear search, improving
 * performance, especially in large datasets.
 *
 * The function also supports ignoring specific indices in the dataset, which can be useful
 * in cases where certain data points have low FWHM.
 *
 * @param a The array of data points (MqsRawDataPoint_t) to search through for a peak.
 * @param size The size of the array.
 * @param l The starting index of the current search window.
 * @param r The ending index of the current search window.
 * @param peakIndex A pointer to store the index of the found peak.
 * @param ignoreIndices An array of indices to be ignored during the search.
 * @param numIgnoreIndices The number of indices to ignore.
 * @return The value of the peak found, or -1 if no peak is found.
 */
static double findPeakRec(const MqsRawDataPoint_t a[], int size, int l, int r, uint16_t* peakIndex)
{
    if (l > r)
        return -1;

    int mid = (l + r) / 2;
    float max_val = 0.0f;
    int max_index = 0;

    int max_row_index = maxrow(a, size, mid, &max_val, &max_index);

    if (mid == 0 || mid == size - 1)
    {
        *peakIndex = max_row_index;
        return max_val;
    }

    if (max_val < a[mid - 1].phaseAngle)
        return findPeakRec(a, size, l, mid - 1, peakIndex);
    else if (max_val < a[mid + 1].phaseAngle)
        return findPeakRec(a, size, mid + 1, r, peakIndex);
    else
    {
        *peakIndex = max_row_index;
        return max_val;
    }
}

// Function prototypes
/*!
 * @brief Calculate local mean and standard deviation for a given index in the signal.
 *
 * This function computes the local mean and standard deviation of the signal within a specified
 * window size centered at a given index.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param window_size The size of the sliding window for local statistics.
 * @param index The index at which to compute the local statistics.
 * @param local_mean Pointer to store the computed local mean.
 * @param local_std Pointer to store the computed local standard deviation.
 */
static inline void calculate_local_stats(const MqsRawDataPoint_t* signal, int length, int window_size, int index, double* local_mean, double* local_std) {
    int start = fmax(0, index - window_size / 2);
    int end = fmin(length - 1, index + window_size / 2);
    double sum = 0, sum_sq = 0;
    int count = 0;

    for (int i = start; i <= end; ++i) {
        double value = signal[i].phaseAngle;
        sum += value;
        sum_sq += value * value;
        count++;
    }

    *local_mean = sum / count;
    *local_std = sqrt(sum_sq / count - (*local_mean) * (*local_mean));
} //boundary problem 

/*!
 * @brief Detects peaks in the signal using an adaptive gradient method.
 *
 * This function identifies peaks in the given signal by analyzing the gradient within a specified range
 * and using a sliding window to compute local statistics.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param start The starting index for peak detection.
 * @param end The ending index for peak detection.
 * @param window_size The size of the sliding window for local statistics.
 * @param peaks Array to store the detected peak indices.
 * @param num_peaks Pointer to store the number of detected peaks.
 */

//restrict is not the way to go here. 
static void adaptive_gradient_find_peaks(const MqsRawDataPoint_t* signal, int length, int start, int end, int window_size, int*  peaks, int* num_peaks) {
    // Define array sizes based on PEAK_DETECTION_WINDOW_SIZE
    double dY[PEAK_DETECTION_WINDOW_SIZE - 2] = { 0 };
    int S[PEAK_DETECTION_WINDOW_SIZE - 2] = { 0 };
    int ddS[PEAK_DETECTION_WINDOW_SIZE - 4] = { 0 };

    int range_length = end - start + 1;
    double local_mean, local_std;

    if (range_length > PEAK_DETECTION_WINDOW_SIZE) {
#ifdef DEBUG_PRINT
        printf("Range length exceeds maximum allowed size.\n");
#endif
        return;
    }

    for (int i = start + 1; i < end - 1; ++i) { //start ve end passlendiğinden boundary problem olmayacak.
        dY[i - start - 1] = signal[i + 1].phaseAngle - signal[i - 1].phaseAngle;
    }

    for (int i = 0; i < range_length - 2; ++i) {
        S[i] = (dY[i] > 0) - (dY[i] < 0);  // Using subtraction to set -1, 0, or 1
    }

    for (int i = 1; i < range_length - 3; ++i) {
        ddS[i - 1] = S[i] - S[i - 1];
    }

    *num_peaks = 0;
    for (int i = 1; i < range_length - 3; ++i) {
        if (ddS[i - 1] == -2) {
            int peak_index = start + i;  // Adjusting index to account for convolution offset
            calculate_local_stats(signal, length, window_size, peak_index, &local_mean, &local_std);
            double adaptive_threshold = local_mean + 0.9 * local_std;  // Adjusted multiplier
            if (signal[peak_index].phaseAngle > adaptive_threshold) {
                peaks[*num_peaks] = peak_index;  // Adjusted peak index stored here
                (*num_peaks)++;
            }
        }
    }
}


/*!
 * @brief Calculate the width of each detected peak.
 *
 * This function computes the width of each detected peak in the signal based on both the full width
 * at half maximum (FWHM) criterion and the width at 10% of the peak height.
 *
 * Additionally, the function calculates the width at 10% of the peak height to account for noise and
 * to get a better understanding of the broader peak characteristics. This is especially useful in
 * noisy signals where the 10% width helps in assessing the influence of noise and the general shape of
 * the dataset on the peak width.
 *
 * The function determines these widths by finding the points where the signal crosses the 50% and 10%
 * thresholds on either side of the peak.
 *
 * @param signal The input signal array.
 * @param peaks Array of detected peak indices.
 * @param num_peaks The number of detected peaks.
 * @param widths Array to store the computed widths of the peaks.
 * @param length The length of the input signal array.
 */
static void calculate_peak_widths(const MqsRawDataPoint_t* signal, const int* peaks, int num_peaks, double* widths, int length) {
    for (int i = 0; i < num_peaks; i++) {
        int peak = peaks[i];
        double peak_height = signal[peak].phaseAngle;
        double half_max = peak_height / 2.0;

        int start = peak, stop = peak;
        while (start > 0 && signal[start].phaseAngle > peak_height * 0.1) start--;  // 10% height as the threshold
        while (stop < length && signal[stop].phaseAngle > peak_height * 0.1) stop++;

        int left_half_max = peak, right_half_max = peak;
        while (left_half_max > 0 && signal[left_half_max].phaseAngle > half_max) left_half_max--;
        while (right_half_max < length && signal[right_half_max].phaseAngle > half_max) right_half_max++;

        widths[i] = (right_half_max - left_half_max) * 0.5 + (stop - start) * 0.5;
    }
} //no boundary problem 

/*!
 * @brief Calculate the prominence of a given peak.
 *
 * This function computes the prominence of a peak in the signal, which is the height of the peak
 * relative to the lowest contour line surrounding it.
 *
 * Additionally, the function adjusts the prominence by a shape factor. The shape factor is based
 * on the deviation of the peak from the local mean and standard deviation within a specified window.
 * This adjustment helps to account for the influence of noise and the general shape of the dataset
 * on the peak prominence.
 *
 * The function determines the prominence by finding the lowest point between the peak and the nearest
 * higher peaks or the signal boundaries on both sides.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param peakIndex The index of the peak for which to calculate prominence.
 * @return The computed prominence of the peak.
 */
static double find_prominence(const MqsRawDataPoint_t* signal, int length, int peakIndex) { //no boundary problem 
    int leftBoundary = 0;
    int rightBoundary = length - 1;
    double peak_val = signal[peakIndex].phaseAngle;

    // Find the nearest higher peak or end on the left
    for (int i = peakIndex - 1; i >= 0; i--) {
        if (signal[i].phaseAngle > peak_val) {
            leftBoundary = i;
            break;
        }
    }

    // Find the nearest higher peak or end on the right
    for (int i = peakIndex + 1; i < length; i++) {
        if (signal[i].phaseAngle > peak_val) {
            rightBoundary = i;
            break;
        }
    }

    // Find the minimum value within the boundaries
    double minValue = signal[rightBoundary].phaseAngle;
    for (int i = leftBoundary; i <= rightBoundary; i++) {
        if (signal[i].phaseAngle < minValue) {
            minValue = signal[i].phaseAngle;
        }
    }

    double local_mean, local_std;
    calculate_local_stats(signal, length, WINDOW_SIZE, peakIndex, &local_mean, &local_std);
    double shape_factor = exp(-(signal[peakIndex].phaseAngle - local_mean) * (signal[peakIndex].phaseAngle - local_mean) / (2 * local_std * local_std));
    double adjusted_prominence = (peak_val - minValue) * shape_factor;

    return adjusted_prominence;
}


/*!
 * @brief Calculate the prominences of all detected peaks.
 *
 * This function computes the prominence of each detected peak in the signal.
 *
 * @param signal The input signal array.
 * @param peaks Array of detected peak indices.
 * @param num_peaks The number of detected peaks.
 * @param prominences Array to store the computed prominences of the peaks.
 * @param length The length of the input signal array.
 */
static void calculate_peak_prominences(const MqsRawDataPoint_t* signal, const int* peaks, int num_peaks, double* prominences, int length) {
    for (int i = 0; i < num_peaks; i++) {
        prominences[i] = find_prominence(signal, length, peaks[i]);
    }
}

/*!
 * @brief Find the primary peak in the signal.
 *
 * This function identifies the primary peak in the signal by recursively searching for the maximum
 * value.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @return The index of the primary peak.
 */
static uint16_t find_primary_peak(const MqsRawDataPoint_t* signal, int length) {
    uint16_t primary_peak = 0;
    float peak_value = findPeakRec(signal, length, 0, length - 1, &primary_peak);
#ifdef DEBUG_PRINT
    printf("Primary peak at index %d with value %f\n", primary_peak, peak_value);
#endif
    return primary_peak;
}

/*!
 * @brief Detect peaks within a specified range in the signal.
 *
 * This function identifies peaks within a specified range in the signal using an adaptive gradient method.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param range_start The starting index of the range for peak detection.
 * @param range_end The ending index of the range for peak detection.
 * @param window_size The size of the sliding window for local statistics.
 * @param peaks Array to store the detected peak indices.
 * @return The number of detected peaks.
 */
static int detect_peaks_in_range(const MqsRawDataPoint_t* signal, int length, int range_start, int range_end, int window_size, int* peaks) {
    int num_peaks = 0;
    adaptive_gradient_find_peaks(signal, length, range_start, range_end, window_size, peaks, &num_peaks);

#ifdef DEBUG_PRINT
    printf("Detected Peaks:\n");
    for (int i = 0; i < num_peaks; i++) {
        printf("Peak at index %d, value = %f\n", peaks[i], signal[peaks[i]].phaseAngle);
    }
#endif
    return num_peaks;
}

/*!
 * @brief Detects the peak with the widest width in the signal.
 *
 * This function identifies the peak with the widest width in the signal using the adaptive gradient method
 * and verifies its prominence.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param peak_range The range for detecting peaks.
 * @param widest_peak_index Pointer to store the index of the peak with the widest width.
 * @return True if a valid peak with sufficient width and prominence is found; false otherwise.
 */
static bool find_widest_peak(const MqsRawDataPoint_t* signal, int length, int* peaks, int num_peaks, double* widths, uint16_t* widest_peak_index) {
    if (num_peaks == 0) {
        return false;
    }

    calculate_peak_widths(signal, peaks, num_peaks, widths, length);

    double* prominences = (double*)malloc(num_peaks * sizeof(double));
    if (!prominences) {
        fprintf(stderr, "Memory allocation failed\n");
        return false;
    }

    calculate_peak_prominences(signal, peaks, num_peaks, prominences, length);

#ifdef DEBUG_PRINT
    printf("Peak Widths and Prominences:\n");
    for (int i = 0; i < num_peaks; i++) {
        printf("Peak at index %d has width %f and prominence %f\n", peaks[i], widths[i], prominences[i]);
    }
#endif

    double max_width = 0.0;
    int max_width_index = -1;
    for (int i = 0; i < num_peaks; i++) {
        if (prominences[i] > PROMINENCE_THRESHOLD && widths[i] > max_width) {
            max_width = widths[i];
            max_width_index = i;
        }
    }

    if (max_width_index != -1) {
        if (prominences[max_width_index] < PROMINENCE_THRESHOLD) {
#ifdef DEBUG_PRINT
            printf("Low Prominence\n");
#endif
            free(prominences);
            return false;
        }
        *widest_peak_index = peaks[max_width_index];
        free(prominences);
        return true;
    }

    free(prominences);
    return false;
}

/*!
 * @brief Detects the peak with the widest width in the signal.
 *
 * This function identifies the peak with the widest width in the signal using the adaptive gradient method
 * and verifies its prominence.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param peak_range The range for detecting peaks.
 * @param widest_peak_index Pointer to store the index of the peak with the widest width.
 * @return True if a valid peak with sufficient width and prominence is found; false otherwise.
 */
static bool detect_peak_with_width(const MqsRawDataPoint_t* signal, int length, int peak_range, uint16_t* widest_peak_index) {
    uint16_t primary_peak = find_primary_peak(signal, length);

    int range_start = primary_peak - peak_range / 2;
    int range_end = primary_peak + peak_range / 2;
    if (range_start < 0) range_start = 0;
    if (range_end >= length) range_end = length - 1;

    int peaks[PEAK_DETECTION_WINDOW_SIZE];
    for (int i = 0; i < peak_range; i++) peaks[i] = -1;  // Initialize peaks array

    int num_peaks = detect_peaks_in_range(signal, length, range_start, range_end, WINDOW_SIZE, peaks);

    double widths[PEAK_DETECTION_WINDOW_SIZE];
    return find_widest_peak(signal, length, peaks, num_peaks, widths, widest_peak_index);
}

/*!
 * @brief Process the input data to find and verify the peak with the widest width.
 *
 * This function processes the input data to identify and verify the peak with the widest width,
 * printing its index and magnitude if found.
 *
 * @param a The input signal array.
 * @param size The size of the input signal array.
 * @param peakIndex Pointer to store the index of the peak with the widest width.
 * @return True if a valid peak with sufficient width and prominence is found; false otherwise.
 */
bool processPeak(MqsRawDataPoint_t a[], int size, uint16_t* peakIndex, bool* isEdgeCase) {
    int peak_range = PEAK_DETECTION_WINDOW_SIZE;  // Example range

    for (int i = 0; i < size; i++) {
        printf(" %f,", a[i].phaseAngle);
    }

    bool result = detect_peak_with_width(a, size, peak_range, peakIndex);

    // Check if peak is near the end and potentially still climaxing
    if (*peakIndex >= size - PEAK_THRESHOLD)
    {
        *isEdgeCase = isPeakClimbing(a, size, *peakIndex, NOISE_TOLERANCE);
    }

    if (result) {
        double peak_magnitude = a[*peakIndex].phaseAngle;

#ifdef DEBUG_PRINT
        printf("Peak found at index: %d\n", *peakIndex);
        printf("Peak magnitude: %f\n", peak_magnitude); 
#endif
    }
    else {
#ifdef DEBUG_PRINT
        printf("No peak with width detected.\n");
#endif
    }

    return result;
}

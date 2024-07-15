#ifndef PEAK_DETECTION_H
#define PEAK_DETECTION_H

#include "../mqs_def.h"
#include <stdbool.h>


#define WINDOW_SIZE 51  // Size of the sliding window for local statistics

/*!
 * @brief Defines the noise tolerance level for validating edge case climbing peaks.
 *
 * This constant represents the threshold for noise tolerance used in determining whether a peak
 * is still climbing at the end of a dataset. It is used in the context of peak analysis to
 * distinguish between genuine rising peaks and minor fluctuations that could be attributed
 * to noise. A lower value indicates stricter criteria for a peak to be considered as climbing.
 */
#define NOISE_TOLERANCE 0.9f 

 /*!
  * @brief Defines the threshold for identifying edge case peaks in a dataset.
  *
  * This constant sets the threshold for defining edge case peaks. An edge case peak is
  * identified when the peak is near the end of the dataset and has not reached its maximum
  * (climax) within the dataset's interval. This threshold value determines how close to the
  * end of the dataset a peak must be to be considered an edge case. It is used to decide
  * whether to check if a peak is still climbing or if it may continue in a subsequent dataset.
  */	
#define PEAK_THRESHOLD  30 

#define PROMINENCE_THRESHOLD 9.0

#define PEAK_DETECTION_WINDOW_SIZE 51

bool processPeak(MqsRawDataPoint_t a[], int size, uint16_t* peakIndex, bool* isEdgeCase);

#endif // PEAK_DETECTION_H

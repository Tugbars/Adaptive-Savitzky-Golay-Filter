#ifndef ADAPTIVE_FILTERING_H
#define ADAPTIVE_FILTERING_H

// Preprocessor constants
#define NOISE_TYPE 'G' //GAUSSIAN

#define SIGMA -1
#define MAX_SIGNAL_LENGTH 501 //MES SWEEP LENGTH

#include "../mqs_def.h"
#include "stdbool.h"
#include "adaptive_filtering_config.h"

// Global variables for peak index and result
extern uint16_t g_peakIndex;
extern bool g_peakResult;

// Enumeration for the states in the denoising process
typedef enum {
    DEN_STATE_WAITING,
    DEN_STATE_INIT,
    DEN_STATE_FIND_PEAK_RANGE,
    DEN_STATE_EVAL_OPTIMAL_ORDER,
    DEN_STATE_APPLY_FILTER_WITH_CACHE,
    DEN_STATE_APPLY_FILTER,
    DEN_STATE_PROCESS_PEAK,
    DEN_STATE_COUNT // This should always be the last entry
} DenState_t;

// Context structure to hold state and other necessary data
typedef struct {
    DenState_t current_state;
    double sigma;
    double best_smoothness;
    double best_correlation;
    int start;
    int end;
    int best_order;
    int best_window;
    int len;
    int interval_size;
    void (*callback)(void);
} DenoiseContext;

// Function to encapsulate the initialization and state machine start-up
void startDenoisingProcess(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, size_t len, void (*callback)(void));

void populate_noisy_sig(MqsRawDataPoint_t* noisy_sig, const double* dataset, size_t dataSize);
#endif // ADAPTIVE_FILTERING_H

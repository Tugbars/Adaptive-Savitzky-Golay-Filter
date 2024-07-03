#ifndef ADAPTIVE_FILTERING_H
#define ADAPTIVE_FILTERING_H

// Preprocessor constants
#define NOISE_TYPE 'G' //GAUSSIAN
#define LAMBDA 0.5
#define SIGMA -1
#define M 5 //MINIMAL WINDOW
#define PMAX 5 //MAX ORDER THAT NEEDS TO BE TESTED. THIS WAS ORIGINALLY 5 AFAIK.
#define PMIN 3 //MAX ORDER THAT NEEDS TO BE TESTED. THIS WAS ORIGINALLY 5 AFAIK.

#define MAX_SIGNAL_LENGTH 501 //MES SWEEP LENGTH

// Enumeration for the states in the denoising process
typedef enum {
    DEN_STATE_INIT,
    DEN_STATE_FIND_PEAK_RANGE,
    DEN_STATE_EVAL_OPTIMAL_ORDER,
    DEN_STATE_APPLY_FILTER,
    DEN_STATE_PROCESS_PEAK, 
    DEN_STATE_DONE,
    DEN_STATE_COUNT // Keep this last to count the states
} DenState_t;

// Context structure to hold state and other necessary data
typedef struct {
    DenState_t current_state;
    double sigma;
    int start;
    int end;
    int best_order;
    int best_window;
    double best_smoothness;
    double best_correlation;
    int len; // Length of the signal
} DenoiseContext;

// Function to encapsulate the initialization and state machine start-up
void startDenoisingProcess(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, size_t len);

void populate_noisy_sig(MqsRawDataPoint_t* noisy_sig, const double* dataset, size_t dataSize);
#endif // ADAPTIVE_FILTERING_H

/**
 * @file mes_adaptive_filtering.c
 * @brief Adaptive Filtering and Peak Processing for Signal Denoising
 *
 * This file contains the implementation of adaptive filtering and peak processing for signal denoising.
 * The denoising process involves multiple states including initialization, peak range finding,
 * evaluating optimal filter order, applying the filter, and processing the peak.
 *
 * Author: Tugbars Heptaskin
 * Date: 07/03/2024
 * Company: Aminic Aps
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "../savitzky_golay_filter/mes_savgol.h"
#include "../adaptive_window/adaptive_filtering_window.h"
#include "../adaptive_order/mes_adaptive_order.h"
#include "mes_adaptive_filtering.h"
#include "adaptive_filtering_config.h"
#include "../adaptive_peak_finder/adaptive_peak_finding.h"


 // Global array for optimal order windows
static OptimalOrderWindow optimal_order_windows[ORDER_EVAL_INTERVAL]; // ORDER_EVAL_INTERVAL
static DenoiseContext ctx;

typedef struct {
    bool isFilteringRequested;
} Status_t;

static Status_t currentStatus = { false };

// Global variables for peak index and result
uint16_t g_peakIndex = 0;
bool g_peakResult = false;

typedef void (*OnEntry_t)(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
typedef void (*OnExit_t)(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);

/**
 * @struct StateFuncs_t
 * @brief Structure to hold function pointers for state entry and exit functions.
 */
typedef struct {
    OnEntry_t onEntry; /**< Function pointer for state entry function */
    OnExit_t onExit;   /**< Function pointer for state exit function */
    bool isComplete;   /**< Flag to indicate if the state processing is complete */
} StateFuncs_t;

// Function prototypes for state entry and exit functions
static void doNothing(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onEntryInit(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onExitInit(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onEntryFindPeakRange(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onExitFindPeakRange(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onEntryEvalOptimalOrder(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onExitEvalOptimalOrder(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onEntryApplyFilterWithCache(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onExitApplyFilterWithCache(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onEntryApplyFilter(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onExitApplyFilter(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onEntryProcessPeak(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
static void onExitProcessPeak(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);

/**
 * @brief State function pointers for each state in the denoising process.
 */
static StateFuncs_t STATE_FUNCS[DEN_STATE_COUNT] = {
    {doNothing, doNothing, true},                                          // DEN_STATE_WAITING
    {onEntryInit, onExitInit, false},                                      // DEN_STATE_INIT
    {onEntryFindPeakRange, onExitFindPeakRange, false},                    // DEN_STATE_FIND_PEAK_RANGE
    {onEntryEvalOptimalOrder, onExitEvalOptimalOrder, false},              // DEN_STATE_EVAL_OPTIMAL_ORDER
    {onEntryApplyFilterWithCache, onExitApplyFilterWithCache, false},      // DEN_STATE_APPLY_FILTER_WITH_CACHE
    {onEntryApplyFilter, onExitApplyFilter, false},                        // DEN_STATE_APPLY_FILTER
    {onEntryProcessPeak, onExitProcessPeak, false}                         // DEN_STATE_PROCESS_PEAK
};

static const char* STATE_NAMES[DEN_STATE_COUNT] = {
    "DEN_STATE_WAITING",
    "DEN_STATE_INIT",
    "DEN_STATE_FIND_PEAK_RANGE",
    "DEN_STATE_EVAL_OPTIMAL_ORDER",
    "DEN_STATE_APPLY_FILTER_WITH_CACHE",
    "DEN_STATE_APPLY_FILTER",
    "DEN_STATE_PROCESS_PEAK"
};

/**
 * @brief Determines the next state in the denoising state machine.
 *
 * @return The next state in the denoising process.
 */
static DenState_t NextState(void) {
    if (ctx.current_state != DEN_STATE_WAITING && !STATE_FUNCS[ctx.current_state].isComplete) {
        return ctx.current_state;  // Stay in the current state if it is not complete
    }

    switch (ctx.current_state) {
    case DEN_STATE_WAITING:
        if (currentStatus.isFilteringRequested)
            return DEN_STATE_INIT;
        return DEN_STATE_WAITING;
    case DEN_STATE_INIT:
        return DEN_STATE_FIND_PEAK_RANGE;
    case DEN_STATE_FIND_PEAK_RANGE:
        return DEN_STATE_EVAL_OPTIMAL_ORDER;
    case DEN_STATE_EVAL_OPTIMAL_ORDER:
        return DEN_STATE_APPLY_FILTER_WITH_CACHE;
    case DEN_STATE_APPLY_FILTER_WITH_CACHE:
        return DEN_STATE_APPLY_FILTER;
    case DEN_STATE_APPLY_FILTER:
        return DEN_STATE_PROCESS_PEAK;
    case DEN_STATE_PROCESS_PEAK:
        return DEN_STATE_WAITING;
    default:
        return DEN_STATE_WAITING; // Default to WAITING to end the state machine
    }
}
void resetFilteringRequest() {
    currentStatus.isFilteringRequested = false;
}

/**
 * @brief Processes state changes in the denoising state machine.
 *
 * This function handles the transition between states in the denoising process.
 * It ensures that the appropriate entry and exit functions are called for each state transition.
 */
static void ProcessStateChange(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    bool isChanged;
    do {
        const DenState_t next = NextState();
        isChanged = (next != ctx.current_state);

        if (isChanged) {
            printf("Transitioning from %s to %s\n", STATE_NAMES[ctx.current_state], STATE_NAMES[next]);

            if (STATE_FUNCS[ctx.current_state].onExit) {
                STATE_FUNCS[ctx.current_state].onExit(noisy_sig, smoothed_sig);
            }

            ctx.current_state = next;
            STATE_FUNCS[ctx.current_state].isComplete = false;  // Reset the completion flag for the new state

            if (STATE_FUNCS[ctx.current_state].onEntry) {
                STATE_FUNCS[ctx.current_state].onEntry(noisy_sig, smoothed_sig);
            }
        }

        // Check the completion flag to ensure the state processing is complete before transitioning
    } while (isChanged && STATE_FUNCS[ctx.current_state].isComplete);

    // Call the callback function when the state machine reaches WAITING state
    if (ctx.current_state == DEN_STATE_WAITING && ctx.callback) {
        void (*doneCallback)(void) = ctx.callback; // Copy the callback pointer
        ctx.callback = NULL; // Clear the callback to prevent re-entry
        doneCallback();

        // Reset the completion flags for all states except WAITING
        for (int i = 0; i < DEN_STATE_COUNT; i++) {
            if (i != DEN_STATE_WAITING) {
                STATE_FUNCS[i].isComplete = false;
            }
        }

        resetFilteringRequest();
    }
}

/**
 * @brief Starts the denoising process.
 *
 * This function initializes the DenoiseContext structure and begins the state machine process.
 */
void startDenoisingProcess(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, uint16_t len, void (*callback)(void)) {
    ctx = (DenoiseContext){
        .current_state = DEN_STATE_WAITING,
        .sigma = -1,
        .best_smoothness = 0.0,
        .best_correlation = 0.0,
        .start = 0,
        .end = 0,
        .best_order = 0,
        .best_window = 0,
        .len = len,
        .interval_size = ORDER_EVAL_INTERVAL, // Default interval size, can be set as needed
        .callback = callback
    };
    currentStatus.isFilteringRequested = true;
    STATE_FUNCS[DEN_STATE_WAITING].isComplete = true;
    ProcessStateChange(noisy_sig, smoothed_sig);
}

static void doNothing(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    // This function intentionally left blank
}

/**
 * @brief Entry function for the INIT state.
 */
static void onEntryInit(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    printf("Sigma calculation\n");
    if (ctx.sigma == -1) {
        ctx.sigma = calculate_sigma(noisy_sig, ctx.len, NOISE_TYPE);
    }
    STATE_FUNCS[ctx.current_state].isComplete = true;
    ProcessStateChange(noisy_sig, smoothed_sig);
}

void onExitInit(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    // Any exit actions for INIT state if needed
}

/**
 * @brief Entry function for the FIND_PEAK_RANGE state.
 */
static void onEntryFindPeakRange(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    find_peak_range(noisy_sig, ctx.len, ctx.interval_size, &ctx.start, &ctx.end);

    if (ctx.start < 0 || ctx.end >= ctx.len) {
        printf("Error: Peak range out of bounds. Start: %d, End: %d, Length: %d\n", ctx.start, ctx.end, ctx.len);
    }
    else {
        STATE_FUNCS[ctx.current_state].isComplete = true;
        ProcessStateChange(noisy_sig, smoothed_sig);
    }
}

static void onExitFindPeakRange(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    // Any exit actions for FIND_PEAK_RANGE state if needed
}

/**
 * @brief Entry function for the EVAL_OPTIMAL_ORDER state.
 */
static void onEntryEvalOptimalOrder(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    uint16_t peakIndex;
    double peakVal = find_primary_peak(noisy_sig, ctx.len, &peakIndex);

    evaluate_optimal_order_for_all_indexes(
        noisy_sig, smoothed_sig, ctx.start, ctx.end,
        ctx.len, ctx.sigma, optimal_order_windows, peakIndex, SMOOTH_CORR_INTERVAL
    );

    printf("Optimal Order Windows: interval_size %d\n", ctx.interval_size);
    for (int i = 0; i < ctx.interval_size; i++) {
        printf("Index %d: Optimal Order = %d, Optimal Window = %d\n",
            i, optimal_order_windows[i].optimal_order, optimal_order_windows[i].optimal_window);
    }

    STATE_FUNCS[ctx.current_state].isComplete = true;
    ProcessStateChange(noisy_sig, smoothed_sig);
}

/**
 * @brief Exit function for the EVAL_OPTIMAL_ORDER state.
 */
static void onExitEvalOptimalOrder(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    // Any exit actions for EVAL_OPTIMAL_ORDER state if needed
}

/**
 * @brief Entry function for the APPLY_FILTER_WITH_CACHE state.
 */
static void onEntryApplyFilterWithCache(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    printf("optimal_filter EVAL:\n");
    apply_optimal_filter_with_cache(
        noisy_sig, smoothed_sig, ctx.len, ctx.start, ctx.end,
        optimal_order_windows, &ctx.best_smoothness, &ctx.best_correlation,
        &ctx.best_order, &ctx.best_window
    );

    STATE_FUNCS[ctx.current_state].isComplete = true;
    ProcessStateChange(noisy_sig, smoothed_sig);
}

static void onExitApplyFilterWithCache(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    // Any exit actions for APPLY_FILTER_WITH_CACHE state if needed
}

/**
 * @brief Entry function for the APPLY_FILTER state.
 */
static void onEntryApplyFilter(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    apply_optimal_filter(noisy_sig, smoothed_sig, ctx.len, ctx.best_order, ctx.best_window);

    STATE_FUNCS[ctx.current_state].isComplete = true;
    ProcessStateChange(noisy_sig, smoothed_sig);
}

static void onExitApplyFilter(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    // Any exit actions for APPLY_FILTER state if needed
}

/**
 * @brief Entry function for the PROCESS_PEAK state.
 */
static void onEntryProcessPeak(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    bool isEdgeCase;
    g_peakResult = processPeak(smoothed_sig, ctx.len, &g_peakIndex, &isEdgeCase);

    if (g_peakResult) {
        printf("Peak processing completed successfully. Peak index: %d, Edge case: %d\n", g_peakIndex, isEdgeCase);
    }
    else {
        printf("No peak found during processing.\n");
    }

    STATE_FUNCS[ctx.current_state].isComplete = true;
    ProcessStateChange(noisy_sig, smoothed_sig);
}

static void onExitProcessPeak(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    printf("Exiting ProcessPeak State\n");
    // Any exit actions for PROCESS_PEAK state if needed
}

/**
 * @brief Populates the noisy signal array with data from the dataset.
 */
void populate_noisy_sig(MqsRawDataPoint_t* noisy_sig, const double* dataset, size_t dataSize) {
    for (size_t i = 0; i < dataSize; ++i) {
        noisy_sig[i].phaseAngle = dataset[i];
        noisy_sig[i].impedance = 0.0;  // Set the impedance to a default value
    }
}


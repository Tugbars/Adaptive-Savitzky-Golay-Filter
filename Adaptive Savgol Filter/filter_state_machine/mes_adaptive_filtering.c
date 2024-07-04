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

DenoiseContext ctx;

typedef void (*OnEntry_t)(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
typedef void (*OnExit_t)(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);

/**
 * @struct StateFuncs_t
 * @brief Structure to hold function pointers for state entry and exit functions.
 */
typedef struct {
    OnEntry_t onEntry; /**< Function pointer for state entry function */
    OnExit_t onExit;   /**< Function pointer for state exit function */
} StateFuncs_t;


// Function prototypes for state entry and exit functions
void onEntryInit(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
void onExitInit(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
void onEntryFindPeakRange(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
void onExitFindPeakRange(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
void onEntryEvalOptimalOrder(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
void onExitEvalOptimalOrder(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
void onEntryApplyFilter(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
void onExitApplyFilter(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
void onEntryProcessPeak(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);
void onExitProcessPeak(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig);

/**
 * @brief State function pointers for each state in the denoising process.
 */
static const StateFuncs_t STATE_FUNCS[DEN_STATE_COUNT] = {
    {onEntryInit, onExitInit},
    {onEntryFindPeakRange, onExitFindPeakRange},
    {onEntryEvalOptimalOrder, onExitEvalOptimalOrder},
    {onEntryApplyFilter, onExitApplyFilter},
    {onEntryProcessPeak, onExitProcessPeak},
    {NULL, NULL} // DEN_STATE_DONE does not require entry/exit functions
};

/**
 * @brief Processes state changes in the denoising state machine.
 *
 * This function handles the transition between states in the denoising process.
 * It ensures that the appropriate entry and exit functions are called for each state transition.
 */
static void ProcessStateChange(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    bool isChanged;
    do {
        const DenState_t next = NextState(&ctx);
        isChanged = (next != ctx.current_state);

        if (isChanged) {
            if (STATE_FUNCS[ctx.current_state].onExit) {
                STATE_FUNCS[ctx.current_state].onExit(&ctx, noisy_sig, smoothed_sig);
            }

            ctx.current_state = next;

            if (STATE_FUNCS[ctx.current_state].onEntry) {
                STATE_FUNCS[ctx.current_state].onEntry(&ctx, noisy_sig, smoothed_sig);
            }
        }
    } while (isChanged);
}

/**
 * @brief Determines the next state in the denoising state machine.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 * @return The next state in the denoising process.
 */
static DenState_t NextState(DenoiseContext* ctx) {
    switch (ctx->current_state) {
    case DEN_STATE_INIT:
        return DEN_STATE_FIND_PEAK_RANGE;
    case DEN_STATE_FIND_PEAK_RANGE:
        return DEN_STATE_EVAL_OPTIMAL_ORDER;
    case DEN_STATE_EVAL_OPTIMAL_ORDER:
        return DEN_STATE_APPLY_FILTER;
    case DEN_STATE_APPLY_FILTER:
        return DEN_STATE_PROCESS_PEAK; // Transition to peak processing
    case DEN_STATE_PROCESS_PEAK:
        return DEN_STATE_DONE;
    case DEN_STATE_DONE:
        // Handle the final state if needed
        break;
    default:
        return DEN_STATE_DONE; // Default to DONE to end the state machine
    }
    return DEN_STATE_DONE;
}

/**
 * @brief Starts the denoising process.
 *
 * This function initializes the DenoiseContext structure and begins the state machine process.
 */
void startDenoisingProcess(MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig, size_t len) {
    ctx = (DenoiseContext){
        .current_state = DEN_STATE_INIT,
        .sigma = -1,
        .best_smoothness = 0.0,
        .best_correlation = 0.0,
        .start = 0,
        .end = 0,
        .best_order = 0,
        .best_window = 0,
        .len = len // use the provided length
    };

    ProcessStateChange(noisy_sig, smoothed_sig);
}

/**
 * @brief Entry function for the INIT state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onEntryInit(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    if (ctx->sigma == -1) {
        ctx->sigma = calculate_sigma(noisy_sig, ctx->len, NOISE_TYPE);
    }
    ProcessStateChange(noisy_sig, smoothed_sig);
}

void onExitInit(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    // Any exit actions for INIT state if needed
}

/**
 * @brief Entry function for the FIND_PEAK_RANGE state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onEntryFindPeakRange(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    find_peak_range(noisy_sig, ctx->len, &ctx->start, &ctx->end);
    ProcessStateChange(noisy_sig, smoothed_sig);
}

void onExitFindPeakRange(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    // Any exit actions for FIND_PEAK_RANGE state if needed
}

/**
 * @brief Entry function for the EVAL_OPTIMAL_ORDER state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onEntryEvalOptimalOrder(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    double GUE_MSE[RANGE_SIZE][ORDER_RANGE];
    OptimalOrderWindow optimal_order_windows[MAX_RANGE_SIZE];

    evaluate_optimal_order_for_all_indexes(
        noisy_sig, smoothed_sig, ctx->start, ctx->end,
        ctx->len, ctx->sigma, LAMBDA, g_adaptive_filtering_config.pmin, g_adaptive_filtering_config.pmax,
        GUE_MSE, optimal_order_windows
    );

    apply_optimal_filter_with_cache(
        noisy_sig, smoothed_sig, ctx->len, ctx->start, ctx->end,
        optimal_order_windows, &ctx->best_smoothness, &ctx->best_correlation,
        &ctx->best_order, &ctx->best_window
    );

    ProcessStateChange(noisy_sig, smoothed_sig);
}

/**
 * @brief Exit function for the EVAL_OPTIMAL_ORDER state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onExitEvalOptimalOrder(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    // Any exit actions for EVAL_OPTIMAL_ORDER state if needed
}

/**
 * @brief Entry function for the APPLY_FILTER state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onEntryApplyFilter(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    apply_optimal_filter(noisy_sig, smoothed_sig, ctx->len, ctx->best_order, ctx->best_window);

    // The denoised array is not needed as per your instructions
    ProcessStateChange(noisy_sig, smoothed_sig);
}

void onExitApplyFilter(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    // Any exit actions for APPLY_FILTER state if needed
}

/**
 * @brief Entry function for the PROCESS_PEAK state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onEntryProcessPeak(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    //printf("Entering ProcessPeak State\n");

    uint16_t peakIndex;
    bool isEdgeCase;
    bool result = processPeak(smoothed_sig, ctx->len, &peakIndex, &isEdgeCase);

    if (result) {
        printf("Peak processing completed successfully. Peak index: %d, Edge case: %d\n", peakIndex, isEdgeCase);
    }
    else {
        printf("No peak found during processing.\n");
    }

    ProcessStateChange(noisy_sig, smoothed_sig);
}

void onExitProcessPeak(DenoiseContext* ctx, MqsRawDataPoint_t* noisy_sig, MqsRawDataPoint_t* smoothed_sig) {
    printf("Exiting ProcessPeak State\n");
    // Any exit actions for PROCESS_PEAK state if needed
}

/**
 * @brief Populates the noisy signal array with data from the dataset.
 *
 * @param dataset Pointer to the dataset array.
 * @param dataSize Size of the dataset array.
 */
void populate_noisy_sig(MqsRawDataPoint_t* noisy_sig, const double* dataset, size_t dataSize) {
    for (size_t i = 0; i < dataSize; ++i) {
        noisy_sig[i].phaseAngle = dataset[i];
        noisy_sig[i].impedance = 0.0;  // Set the impedance to a default value
    }
}

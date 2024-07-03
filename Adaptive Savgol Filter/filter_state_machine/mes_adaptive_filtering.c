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
#include "savitzky_golay_filter/mes_savgol.h"
#include "adaptive_order/mes_adaptive_order.h"
#include "adaptive_window/adaptive_filtering_window.h"

#include "mes_adaptive_filtering.h"

 // Global variables
MqsRawDataPoint_t noisy_sig[MAX_SIGNAL_LENGTH];
MqsRawDataPoint_t smoothed_sig[MAX_SIGNAL_LENGTH];
DenoiseContext ctx;

typedef void (*OnEntry_t)(DenoiseContext* ctx);
typedef void (*OnExit_t)(DenoiseContext* ctx);

/**
 * @struct StateFuncs_t
 * @brief Structure to hold function pointers for state entry and exit functions.
 */
typedef struct {
	OnEntry_t onEntry; /**< Function pointer for state entry function */
	OnExit_t onExit;   /**< Function pointer for state exit function */
} StateFuncs_t;


// Function prototypes for state entry and exit functions
void onEntryInit(DenoiseContext* ctx);
void onExitInit(DenoiseContext* ctx);
void onEntryFindPeakRange(DenoiseContext* ctx);
void onExitFindPeakRange(DenoiseContext* ctx);
void onEntryEvalOptimalOrder(DenoiseContext* ctx);
void onExitEvalOptimalOrder(DenoiseContext* ctx);
void onEntryApplyFilter(DenoiseContext* ctx);
void onExitApplyFilter(DenoiseContext* ctx);
void onEntryProcessPeak(DenoiseContext* ctx); // New entry function for peak processing
void onExitProcessPeak(DenoiseContext* ctx);  // New exit function for peak processing

/**
 * @brief State function pointers for each state in the denoising process.
 */
static const StateFuncs_t STATE_FUNCS[DEN_STATE_COUNT] = {
	{onEntryInit, onExitInit},
	{onEntryFindPeakRange, onExitFindPeakRange},
	{onEntryEvalOptimalOrder, onExitEvalOptimalOrder},
	{onEntryApplyFilter, onExitApplyFilter},
	{onEntryProcessPeak, onExitProcessPeak}, // New state function entry
	{NULL, NULL} // DEN_STATE_DONE does not require entry/exit functions
};

/**
 * @brief Processes state changes in the denoising state machine.
 *
 * This function handles the transition between states in the denoising process.
 * It ensures that the appropriate entry and exit functions are called for each state transition.
 */
static void ProcessStateChange(void) {
	bool isChanged;
	do {
		const DenState_t next = NextState(&ctx);
		isChanged = (next != ctx.current_state);

		if (isChanged) {
			if (STATE_FUNCS[ctx.current_state].onExit) {
				STATE_FUNCS[ctx.current_state].onExit(&ctx);
			}

			ctx.current_state = next;

			if (STATE_FUNCS[ctx.current_state].onEntry) {
				STATE_FUNCS[ctx.current_state].onEntry(&ctx);
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
void startDenoisingProcess() {
	ctx = (DenoiseContext){
		.current_state = DEN_STATE_INIT,
		.sigma = -1,
		.best_smoothness = 0.0,
		.best_correlation = 0.0,
		.start = 0,
		.end = 0,
		.best_order = 0,
		.best_window = 0,
		.len = 360 // or any other appropriate length
	};

	ProcessStateChange(&ctx);
}

/**
 * @brief Entry function for the INIT state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onEntryInit(DenoiseContext* ctx) {
	if (SIGMA == -1) {
		ctx->sigma = calculate_sigma(noisy_sig, ctx->len, NOISE_TYPE);
	}
	ProcessStateChange(ctx);
}

void onExitInit(DenoiseContext* ctx) {
	// Any exit actions for INIT state if needed
}

/**
 * @brief Entry function for the FIND_PEAK_RANGE state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onEntryFindPeakRange(DenoiseContext* ctx) {
	find_peak_range(noisy_sig, ctx->len, &ctx->start, &ctx->end);
	ProcessStateChange(ctx);
}

void onExitFindPeakRange(DenoiseContext* ctx) {
	// Any exit actions for FIND_PEAK_RANGE state if needed
}

/**
 * @brief Entry function for the EVAL_OPTIMAL_ORDER state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onEntryEvalOptimalOrder(DenoiseContext* ctx) {
	double GUE_MSE[RANGE_SIZE][ORDER_RANGE];
	OptimalOrderWindow optimal_order_windows[MAX_RANGE_SIZE];

	evaluate_optimal_order_for_all_indexes(
		noisy_sig, smoothed_sig, ctx->start, ctx->end,
		ctx->len, ctx->sigma, LAMBDA, PMIN, PMAX, //dogru mu diye bak. bu veriler. 
		GUE_MSE, optimal_order_windows
	);

	apply_optimal_filter_with_cache(
		noisy_sig, smoothed_sig, ctx->len, ctx->start, ctx->end,
		optimal_order_windows, &ctx->best_smoothness, &ctx->best_correlation,
		&ctx->best_order, &ctx->best_window
	);

	ProcessStateChange(ctx);
}

/**
 * @brief Exit function for the EVAL_OPTIMAL_ORDER state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onExitEvalOptimalOrder(DenoiseContext* ctx) {
	// Any exit actions for EVAL_OPTIMAL_ORDER state if needed
}

/**
 * @brief Exit function for the APPLY_FILTER state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onEntryApplyFilter(DenoiseContext* ctx) {
	apply_optimal_filter(noisy_sig, smoothed_sig, ctx->len, ctx->best_order, ctx->best_window);

	// The denoised array is not needed as per your instructions
	ProcessStateChange(ctx);
}

void onExitApplyFilter(DenoiseContext* ctx) {
	// Any exit actions for APPLY_FILTER state if needed
}

/**
 * @brief Entry function for the PROCESS_PEAK state.
 *
 * @param ctx Pointer to the DenoiseContext structure.
 */
void onEntryProcessPeak(DenoiseContext* ctx) {
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

	ProcessStateChange();
}

void onExitProcessPeak(DenoiseContext* ctx) {
	printf("Exiting ProcessPeak State\n");
	// Any exit actions for PROCESS_PEAK state if needed
}

/**
 * @brief Populates the noisy signal array with data from the dataset.
 *
 * @param dataset Pointer to the dataset array.
 * @param dataSize Size of the dataset array.
 */
void populate_noisy_sig(const double* dataset, size_t dataSize) {
	for (size_t i = 0; i < dataSize; ++i) {
		noisy_sig[i].phaseAngle = dataset[i];
		noisy_sig[i].impedance = 0.0;  // Set the impedance to a default value
	}
}

//how can I 
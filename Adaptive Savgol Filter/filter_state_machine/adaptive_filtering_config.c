#include "adaptive_filtering_config.h"

AdaptiveFilteringConfig g_adaptive_filtering_config;

// Ensure window sizes are odd and orders are >= 2
static int validate_and_adjust_config(int* min_window, int* max_window, int* pmin, int* pmax) {
    int valid = 1;

    if (*min_window % 2 == 0) {
        //printf("Warning: min_window %d is not an odd number. Adjusting to %d.\n", *min_window, *min_window + 1);
        *min_window += 1;
        valid = 0;
    }
    if (*max_window % 2 == 0) {
        //printf("Warning: max_window %d is not an odd number. Adjusting to %d.\n", *max_window, *max_window + 1);
        *max_window += 1;
        valid = 0;
    }
    if (*pmin < 2) {
        //printf("Warning: pmin %d is less than 2. Adjusting to 2.\n", *pmin);
        *pmin = 2;
        valid = 0;
    }
    if (*pmax < 2) {
        //printf("Warning: pmax %d is less than 2. Adjusting to 2.\n", *pmax);
        *pmax = 2;
        valid = 0;
    }

    return valid;
}

void initialize_adaptive_filtering_config(int min_window, int max_window, int pmin, int pmax) {
    if (!validate_and_adjust_config(&min_window, &max_window, &pmin, &pmax)) {
        //printf("Error: Invalid configuration provided.\n");
        return;
    }

    g_adaptive_filtering_config.min_window = min_window;
    g_adaptive_filtering_config.max_window = max_window;
    g_adaptive_filtering_config.pmin = pmin;
    g_adaptive_filtering_config.pmax = pmax;
}

void update_adaptive_filtering_config(int min_window, int max_window, int pmin, int pmax) {
    if (!validate_and_adjust_config(&min_window, &max_window, &pmin, &pmax)) {
        printf("Error: Invalid configuration provided.\n");
        return;
    }

    g_adaptive_filtering_config.min_window = min_window;
    g_adaptive_filtering_config.max_window = max_window;
    g_adaptive_filtering_config.pmin = pmin;
    g_adaptive_filtering_config.pmax = pmax;
}
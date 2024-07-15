#include "adaptive_filtering_config.h"

AdaptiveFilteringConfig g_adaptive_filtering_config;

// Ensure window sizes are odd and orders are >= 2
static int validate_and_adjust_config(uint8_t* min_window, uint8_t* max_window, uint8_t* pmin, uint8_t* pmax) {
    int valid = 1;

    if (*min_window % 2 == 0) {
        *min_window += 1;
        valid = 0;
    }
    if (*max_window % 2 == 0) {
        *max_window += 1;
        valid = 0;
    }
    if (*pmin < 2) {
        *pmin = 2;
        valid = 0;
    }
    if (*pmax < 2) {
        *pmax = 2;
        valid = 0;
    }

    return valid;
}

void initialize_adaptive_filtering_config(uint8_t min_window, uint8_t max_window, uint8_t pmin, uint8_t pmax) {
    if (!validate_and_adjust_config(&min_window, &max_window, &pmin, &pmax)) {
        return;
    }

    g_adaptive_filtering_config.min_window = min_window;
    g_adaptive_filtering_config.max_window = max_window;
    g_adaptive_filtering_config.pmin = pmin;
    g_adaptive_filtering_config.pmax = pmax;
}

void update_adaptive_filtering_config(uint8_t min_window, uint8_t max_window, uint8_t pmin, uint8_t pmax) {
    if (!validate_and_adjust_config(&min_window, &max_window, &pmin, &pmax)) {
        return;
    }

    g_adaptive_filtering_config.min_window = min_window;
    g_adaptive_filtering_config.max_window = max_window;
    g_adaptive_filtering_config.pmin = pmin;
    g_adaptive_filtering_config.pmax = pmax;
}

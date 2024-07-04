#include "adaptive_filtering_config.h"

AdaptiveFilteringConfig g_adaptive_filtering_config;

void initialize_adaptive_filtering_config(int min_window, int max_window, int pmin, int pmax) {
    g_adaptive_filtering_config.min_window = min_window;
    g_adaptive_filtering_config.max_window = max_window;
    g_adaptive_filtering_config.pmin = pmin;
    g_adaptive_filtering_config.pmax = pmax;
}

void update_adaptive_filtering_config(int min_window, int max_window, int pmin, int pmax) {
    g_adaptive_filtering_config.min_window = min_window;
    g_adaptive_filtering_config.max_window = max_window;
    g_adaptive_filtering_config.pmin = pmin;
    g_adaptive_filtering_config.pmax = pmax;
}

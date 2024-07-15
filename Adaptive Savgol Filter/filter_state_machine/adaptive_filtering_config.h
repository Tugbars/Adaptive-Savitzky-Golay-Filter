#ifndef ADAPTIVE_FILTERING_CONFIG_H
#define ADAPTIVE_FILTERING_CONFIG_H
#include <stdint.h>
typedef struct {
    uint8_t min_window;
    uint8_t max_window;
    uint8_t pmin;
    uint8_t pmax;
} AdaptiveFilteringConfig;

#define MIN_WINDOW 13
#define MAX_WINDOW 23
#define MAX_ITERATIONS ((MAX_WINDOW - MIN_WINDOW) / 2 + 1)

#define LAMBDA 0.5
#define MIN_ORDER 3
#define MAX_ORDER 5

#define ORDER_EVAL_INTERVAL 1
#define SMOOTH_CORR_INTERVAL 81

extern AdaptiveFilteringConfig g_adaptive_filtering_config;

void initialize_adaptive_filtering_config(uint8_t min_window, uint8_t max_window, uint8_t pmin, uint8_t pmax);
void update_adaptive_filtering_config(uint8_t min_window, uint8_t max_window, uint8_t pmin, uint8_t pmax);
#endif // ADAPTIVE_FILTERING_CONFIG_H

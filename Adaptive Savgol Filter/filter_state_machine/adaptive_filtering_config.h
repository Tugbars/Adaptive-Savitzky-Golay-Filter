#ifndef ADAPTIVE_FILTERING_CONFIG_H
#define ADAPTIVE_FILTERING_CONFIG_H

typedef struct {
    int min_window;
    int max_window;
    int pmin;
    int pmax;
} AdaptiveFilteringConfig;

#define MIN_WINDOW 13
#define MAX_WINDOW 23
#define MAX_ITERATIONS ((MAX_WINDOW - MIN_WINDOW) / 2 + 1)

#define MIN_ORDER 3
#define MAX_ORDER 5

extern AdaptiveFilteringConfig g_adaptive_filtering_config;

void initialize_adaptive_filtering_config(int min_window, int max_window, int pmin, int pmax);
void update_adaptive_filtering_config(int min_window, int max_window, int pmin, int pmax);

#endif // ADAPTIVE_FILTERING_CONFIG_H
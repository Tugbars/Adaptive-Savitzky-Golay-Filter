#ifndef ADAPTIVE_FILTERING_CONFIG_H
#define ADAPTIVE_FILTERING_CONFIG_H

typedef struct {
    int min_window;
    int max_window;
    int pmin;
    int pmax;
} AdaptiveFilteringConfig;

extern AdaptiveFilteringConfig g_adaptive_filtering_config;

void initialize_adaptive_filtering_config(int min_window, int max_window, int pmin, int pmax);
void update_adaptive_filtering_config(int min_window, int max_window, int pmin, int pmax);

#endif // ADAPTIVE_FILTERING_CONFIG_H
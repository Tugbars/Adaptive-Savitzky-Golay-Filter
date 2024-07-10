# Savitzky-Golay Adaptive Filtering Implementation

This project provides an implementation of the Savitzky-Golay filter that dynamically determines and applies the most optimal window size and polynomial order for denoising signals. The algorithm adapts to the shape and noise level of the peaks in the data based on a user-modifiable threshold. The implementation is designed for high-performance and precision in analyzing noisy signals, making it particularly suitable for finding peaks in noisy impedance curves. This code is developed to smoothen the impedance curves enough that small noise spikes can't influence finding the peak with the largest width in a given signal.

## Key Features:

- **Optimal Window Selection:** Utilizes smoothness and Pearson correlation coefficient to determine the optimal window size for filtering.
- **Optimal Order Selection:** Uses Generalized Unbiased Estimate of Mean Squared Error (GUE-MSE) to find the optimal polynomial order.
- **Peak Detection:** Specifically developed to find one peak in a narrow range of a noisy impedance analyzer output but can easily be modified to find multiple peaks in full range of the dataset. The algorithm first identifies a primary peak and then evaluates potential peaks around it, considering the local mean and standard deviation. Peak detection algorithm also calculates FWHM(%10 + %50) and prominence of the peak with the consideration of the local shape factor around the peak. 
- **State Machine Architecture:** Utilizes a state machine for better scalability and to allow for future feedback-based recursive search for the optimal window and order.
- **Industry Standard Accuracy:** The output of the code matches with the MATLAB's output. 

## User-Modifiable Thresholds:
The system includes hardcoded thresholds that the user can modify for their needs:

- **Smoothness Threshold:** Defined by `SMOOTHNESS_THRESHOLD` (default `1.1`). A lower value indicates less noise in the signal.
- **Correlation Threshold:** Defined by `CORRELATION_THRESHOLD` (default `0.99`). A value of `1.0` indicates a perfectly aligned smoothed signal compared to the actual noisy signal.
- **Min and Max Window Sizes: These values set the range for the window size in the Savitzky-Golay filter and are defined during the initialization of the adaptive filtering configuration. Example setup:
  ```initialize_adaptive_filtering_config(MIN_WINDOW, MAX_WINDOW, MIN_ORDER, MAX_ORDER);```

- Min and Max Polynomial Orders: Similarly, the range for the polynomial order in the Savitzky-Golay filter is defined during initialization. Example setup:
```initialize_adaptive_filtering_config(5, 31, MIN_ORDER, MAX_ORDER);```

- The algorithm is designed to return the most optimal parameters if the smoothness and correlation thresholds are not met within the defined min-max window and order interval.

## Example Usage:
The implementation is particularly suitable for high-performance and precision analysis of noisy signals, especially in the context of finding peaks in noisy impedance curves. It was specifically developed to find a peak in expectation of having only one peak but can be adapted for process all the peaks in a given dataset. 
```
**Set your callback function, which will be called when the denoising process is completed.**
static void denoisingFinishedCallback(void) {
    //add your callback feature here. 
}

int main() {
    const double dataset[] = { /* insert your data here */ };
    const double dataset2[] = { /* insert your data here */ };

    const size_t dataSize = sizeof(dataset) / sizeof(dataset[0]);
    populate_noisy_sig(noisy_sig, dataset, dataSize); 

    **set the min and max order/window here.**
    initialize_adaptive_filtering_config(MIN_WINDOW, MAX_WINDOW, 3, 5); 

    **start the denoising state machine.**
    startDenoisingProcess(rawData->data, filterData->data, 360, denoisingFinishedCallback);

    return 0; 
}
```

## Output example: 
![image](https://github.com/Tugbars/Adaptive-Savitzky-Golay-Filter/assets/23309063/211370bb-08bb-4286-9fce-f36c64a29dbf)

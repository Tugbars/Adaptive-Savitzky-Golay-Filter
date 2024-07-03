# Savitzky-Golay Adaptive Filtering Implementation
This project provides an implementation of the Savitzky-Golay filter that dynamically determines and applies the most optimal window size and polynomial order for denoising signals. The algorithm adapts to the shape and noise level of the peaks in the data based on a user-modifiable threshold. The implementation is designed for high-performance and precision in analyzing noisy signals, making it particularly suitable for finding peaks in noisy impedance curves.

## Key Features:

- Optimal Window Selection: Utilizes smoothness and Pearson correlation coefficient to determine the optimal window size for filtering.
- Optimal Order Selection: Uses Generalized Unbiased Estimate of Mean Squared Error (GUE-MSE) to find the optimal polynomial order.
- Peak Detection: Specifically developed to find one peak in a narrow range of a noisy impedance analyzer output but can easily be modified to find multiple peaks in a larger range. The algorithm first identifies a primary peak and then evaluates potential peaks around it, considering the local mean and standard deviation. The program can be modified to check all the datapoints in the noisy signal. For my application, this was not necessary. 
- State Machine Architecture: Utilizes a state machine for better scalability and to allow for future feedback-based recursive search for the optimal window and order.

## User-Modifiable Thresholds:
The system includes hardcoded thresholds that the user can modify for their needs:

-Smoothness Threshold: Defined by SMOOTHNESS_THRESHOLD (default 1.1). A lower value indicates less noise in the signal.
-Correlation Threshold: Defined by CORRELATION_THRESHOLD (default 0.99). A value of 1.0 indicates a perfectly aligned smoothed signal compared to the actual noisy signal.
-Min and Max Window Sizes: Defined by MIN_WINDOW (default 5) and MAX_WINDOW (default 31). These values set the range for the window size in the Savitzky-Golay filter.
-Min and Max Polynomial Orders: These values can be initialized by the user to set the range for the polynomial order in the Savitzky-Golay filter.

# Example Usage:
The implementation is particularly suitable for high-performance and precision analysis of noisy signals, especially in the context of finding peaks in noisy impedance curves. It was specifically developed to find a peak in expectation of having only one peak but can be adapted for multiple peaks.

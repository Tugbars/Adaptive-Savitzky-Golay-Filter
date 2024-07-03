# Savitzky-Golay Adaptive Filtering Implementation
This project provides an implementation of the Savitzky-Golay filter that dynamically determines and applies the most optimal window size and polynomial order for denoising signals. The algorithm adapts to the shape and noise level of the peaks in the data based on a user-modifiable threshold. The implementation is designed for high-performance and precision in analyzing noisy signals, making it particularly suitable for finding peaks in noisy impedance curves.

## Key Features:

- Optimal Window Selection: Utilizes smoothness and Pearson correlation coefficient to determine the optimal window size for filtering.
- Optimal Order Selection: Uses Generalized Unbiased Estimate of Mean Squared Error (GUE-MSE) to find the optimal polynomial order.
- Peak Detection: Specifically developed to find one peak in a narrow range of a noisy impedance analyzer output but can easily be modified to find multiple peaks in a larger range.
- State Machine Architecture: Utilizes a state machine for better scalability and to allow for future feedback-based recursive search for the optimal window and order.

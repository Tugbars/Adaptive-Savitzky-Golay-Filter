function adaptiveSavgolTesting()
    % Define default parameters for smoothing
    defaultWindowSize1 = 13;
    defaultPolyOrder1 = 3;
    defaultWindowSize2 = 13;
    defaultPolyOrder2 = 3;
    
    test4045 = [11.272, 11.254, 11.465, 11.269, 11.31, 11.388, 11.385, 11.431, 11.333, 11.437, 11.431, 11.527, 11.483, 11.449, 11.544, 11.39, 11.469, 11.526, 11.498, 11.522, 11.709, 11.503, 11.564, 11.428, 11.714, 11.707, 11.619, 11.751, 11.626, 11.681, 11.838, 11.658, 11.859, 11.916, 11.814, 11.833, 12.046, 11.966, 12.031, 12.079, 11.958, 12.114, 12.041, 12.186, 12.048, 12.258, 12.312, 12.126, 12.159, 12.393, 12.221, 12.45, 12.439, 12.282, 12.373, 12.573, 12.647, 12.545, 12.467, 12.629, 12.686, 12.668, 12.748, 12.71, 12.852, 13.02, 12.848, 13.144, 13.225, 13.211, 13.496, 13.311, 13.33, 13.634, 13.189, 13.623, 13.671, 13.618, 13.645, 13.779, 14.006, 14.13, 14.071, 14.277, 14.223, 14.457, 14.378, 14.698, 14.599, 14.84, 15.143, 15.106, 15.343, 15.506, 15.665, 15.889, 15.878, 16.055, 16.153, 15.966, 16.637, 16.783, 16.746, 17.193, 16.877, 17.656, 17.522, 17.842, 18.086, 18.336, 18.863, 18.977, 19.534, 19.308, 19.626, 19.956, 20.221, 20.673, 20.59, 21.229, 21.767, 22.225, 22.477, 22.695, 22.828, 23.586, 23.776, 24.39, 25.316, 24.639, 25.767, 26.469, 26.976, 27.651, 27.807, 28.089, 28.869, 29.964, 30.367, 30.159, 31.133, 32.034, 33.131, 32.775, 34.372, 34.516, 35.603, 36.214, 37.742, 38.868, 38.702, 39.811, 40.818, 41.422, 41.521, 42.57, 42.819, 42.871, 42.944, 43.851, 44.086, 44.272, 44.466, 44.274, 44.473, 44.348, 43.932, 43.817, 43.48, 42.943, 42.491, 41.793, 41.071, 39.491, 39.231, 38.365, 37.833, 36.583, 35.787, 34.949, 33.006, 32.827, 32.266, 31.012, 30.436, 29.737, 28.097, 28.76, 27.068, 26.195, 25.262, 24.677, 24.211, 23.574, 22.868, 22.781, 22.258, 21.475, 21.247, 21.982, 20.771, 20.383, 20.349, 19.866, 19.433, 18.573, 18.723, 18.325, 18.084, 18.226, 17.492, 17.505, 16.762, 16.907, 16.606, 16.265, 16.234, 15.983, 16.147, 15.811, 15.667, 15.509, 15.325, 15.031, 14.884, 14.881, 14.836, 14.814, 14.706, 14.158, 14.399, 14.123, 14.084, 14.173, 13.963, 13.981, 14.218, 13.898, 13.869, 13.701, 13.397, 13.528, 13.321, 13.071, 13.393, 13.164, 12.876, 13.021, 12.989, 12.869, 13.004, 12.833, 12.795, 12.661, 12.761, 12.547, 12.775, 12.388, 12.425, 12.564, 12.408, 12.301, 12.469, 12.173, 12.323, 12.248, 12.281, 12.208, 11.887, 12.149, 12.073, 12.053, 11.88, 12.066, 11.958, 12.007, 11.868, 11.921, 11.898, 11.804, 11.7, 11.81, 11.758, 11.717, 11.715, 11.611, 11.719, 11.679, 11.619, 11.58, 11.576, 11.589, 11.491, 11.659, 11.506, 11.431, 11.535, 11.349, 11.464, 11.343, 11.492, 11.407, 11.479, 11.269, 11.355, 11.323, 11.341, 11.238, 11.32, 11.333, 11.262, 11.31, 11.221, 11.302, 11.135, 11.139, 11.217, 11.343, 11.225, 11.089, 11.079, 11.127, 11.082, 11.141, 11.186, 11.184, 11.231, 11.025, 11.058, 11.076, 11.087, 11.047, 11.02, 10.996, 10.906, 11.144, 11.005, 10.911, 10.993, 10.858, 11.086, 10.954, 10.906, 11.026, 11.005, 10.934, 10.922, 10.914, 10.955, 11.057, 10.967, 10.811, 10.833, 10.747, 10.821, 10.946, 10.844, 10.838, 10.848, 10.847];
    test4046 = [11.26, 11.13, 11.276, 11.136, 11.194, 11.22, 11.18, 11.371, 11.269, 11.389, 11.377, 11.448, 11.248, 11.569, 11.442, 11.332, 11.525, 11.532, 11.427, 11.48, 11.431, 11.435, 11.515, 11.561, 11.489, 11.575, 11.693, 11.484, 11.619, 11.572, 11.708, 11.676, 11.528, 11.828, 11.598, 11.697, 11.823, 11.839, 11.889, 11.839, 11.927, 11.99, 11.826, 11.964, 11.813, 12.089, 12.076, 11.986, 12.192, 12.118, 11.997, 12.144, 12.311, 12.315, 12.293, 12.441, 12.35, 12.354, 12.399, 12.419, 12.455, 12.453, 12.566, 12.538, 12.725, 12.588, 12.89, 12.914, 12.913, 12.94, 13.023, 12.842, 13.141, 13.083, 13.321, 13.352, 13.372, 13.488, 13.555, 13.674, 13.66, 13.651, 13.997, 13.868, 13.969, 14.105, 14.191, 14.531, 14.278, 14.432, 14.523, 14.744, 14.441, 14.802, 15.066, 14.889, 15.088, 15.493, 15.494, 15.774, 16.081, 15.847, 15.77, 16.463, 16.188, 16.436, 16.388, 16.708, 16.842, 16.95, 17.413, 17.781, 17.929, 18.59, 18.478, 18.806, 18.738, 19.276, 19.744, 19.83, 20.157, 20.391, 20.367, 20.866, 21.19, 21.535, 21.994, 22.578, 22.869, 22.826, 23.654, 23.704, 24.724, 24.866, 25.609, 25.97, 26.53, 26.882, 28.051, 28.477, 29.039, 29.041, 30.509, 30.461, 31.413, 31.605, 32.46, 34.372, 34.084, 35.04, 35.014, 35.321, 36.64, 36.967, 37.323, 38.137, 39.721, 39.566, 40.457, 41.6, 42.555, 42.642, 43.222, 43.819, 43.611, 44.176, 44.248, 44.378, 44.385, 44.624, 44.413, 44.228, 43.884, 43.298, 43.342, 42.24, 41.29, 41.042, 40.263, 39.176, 38.162, 38.101, 36.628, 37.511, 35.735, 34.998, 33.44, 32.395, 31.241, 30.892, 29.863, 29.021, 28.65, 27.799, 26.718, 27.519, 25.865, 25.085, 23.494, 23.679, 22.736, 22.747, 22.706, 21.48, 21.577, 21.855, 20.304, 20.415, 19.763, 19.755, 19.377, 18.424, 18.938, 18.784, 18.147, 18.072, 17.786, 17.284, 17.527, 17.015, 16.833, 16.676, 16.693, 16.288, 16.071, 15.955, 15.876, 15.482, 15.065, 15.203, 15.227, 14.815, 15.08, 14.961, 14.7, 14.835, 14.526, 14.142, 14.214, 14.261, 14.06, 14.094, 13.824, 13.916, 13.689, 13.736, 13.817, 13.489, 13.69, 13.645, 13.332, 13.344, 13.16, 13.13, 13.111, 12.94, 12.997, 12.829, 12.906, 12.569, 12.824, 12.533, 12.456, 12.47, 12.578, 12.363, 12.328, 12.399, 12.285, 12.3, 12.28, 12.232, 12.366, 12.303, 12.174, 12.035, 11.957, 12.123, 11.931, 12.031, 11.943, 12.024, 11.989, 11.93, 11.702, 11.943, 11.827, 11.818, 11.877, 11.696, 11.784, 11.726, 11.617, 11.542, 11.503, 11.58, 11.639, 11.688, 11.514, 11.541, 11.267, 11.388, 11.45, 11.537, 11.489, 11.32, 11.366, 11.376, 11.269, 11.24, 11.417, 11.314, 11.296, 11.306, 11.231, 11.381, 11.173, 11.321, 11.21, 11.185, 11.298, 11.121, 11.287, 11.227, 11.112, 11.199, 11.208, 11.224, 11.21, 11.168, 11.168, 11.266, 11.075, 11.15, 10.992, 11.005, 11.081, 10.916, 10.984, 11.074, 10.954, 11.052, 11.105, 10.999, 10.953, 11.02, 10.945, 11.056, 11.065, 10.913, 10.958, 11.022, 11.038, 10.969, 10.887, 10.904, 10.936, 10.899, 10.988, 10.752];

    % Plot smoothed data with separate parameters for each array
    plotSmoothedData(test4045, test4046, defaultWindowSize1, defaultPolyOrder1, defaultWindowSize2, defaultPolyOrder2);
end

function plotSmoothedData(arr1, arr2, windowSize1, polyOrder1, windowSize2, polyOrder2)
    % Validate windowSize and polyOrder for both arrays
    if mod(windowSize1, 2) == 0
        error('Window size for array 1 must be odd.');
    end
    if windowSize1 < 1 || polyOrder1 >= windowSize1
        error('Invalid window size or polynomial order for array 1.');
    end
    
    if mod(windowSize2, 2) == 0
        error('Window size for array 2 must be odd.');
    end
    if windowSize2 < 1 || polyOrder2 >= windowSize2
        error('Invalid window size or polynomial order for array 2.');
    end
    
    % Apply Savitzky-Golay filter to both arrays with their respective parameters
    smoothArr1 = sgolayfilt(arr1, polyOrder1, windowSize1);
    smoothArr2 = sgolayfilt(arr2, polyOrder2, windowSize2);

    % Find the highest peak in the smoothed data
    [~, maxIdx1] = max(smoothArr1);
    [~, maxIdx2] = max(smoothArr2);
    
    % Plot the smoothed and original data
    figure;
    hold on;
    plot(arr1, 'b:', 'LineWidth', 1); % Original data in dotted blue
    plot(arr2, 'r:', 'LineWidth', 1); % Original data in dotted red
    plot(smoothArr1, 'b', 'LineWidth', 2); % Smoothed data in solid blue
    plot(smoothArr2, 'r', 'LineWidth', 2); % Smoothed data in solid red
    
    % Highlight the highest peaks in the smoothed data
    plot(maxIdx1, smoothArr1(maxIdx1), 'bp', 'MarkerSize', 12, 'MarkerFaceColor', 'b');
    plot(maxIdx2, smoothArr2(maxIdx2), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    
    % Add legend
    legend('Original Array 1', 'Original Array 2', 'Smoothed Array 1', 'Peak 1', 'Smoothed Array 2', 'Peak 2');
    
    % Get the legend handle and position
    lgd = legend();
    legendPos = lgd.Position;
    
    % Calculate annotation positions based on legend position
    annotationXPos = legendPos(1) + legendPos(3) + 0.05;
    annotationYPos1 = legendPos(2) + legendPos(4) - 0.05;
    annotationYPos2 = legendPos(2) + legendPos(4) - 0.10;
    
    % Annotate the Savitzky-Golay filter parameters for each array
    annotation('textbox', [annotationXPos, annotationYPos1, 0.1, 0.1], 'String', ...
        sprintf('Array 1:\nWindow Size = %d\nPoly Order = %d', windowSize1, polyOrder1), 'Color', 'b', 'EdgeColor', 'none', 'FontSize', 10);
    annotation('textbox', [annotationXPos, annotationYPos2, 0.1, 0.1], 'String', ...
        sprintf('Array 2:\nWindow Size = %d\nPoly Order = %d', windowSize2, polyOrder2), 'Color', 'r', 'EdgeColor', 'none', 'FontSize', 10);
    
    title('Data Comparison: Original and Smoothed with Separate Parameters');
    xlabel('Data Points');
    ylabel('Amplitude');
    hold off;
end
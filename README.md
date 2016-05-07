lab-canny-filter
================

Canny edge filter implementation. Extracts sharp edges from any input image.

# Required environment

 * Scilab
 * SIVP Scilab toolbox

# How to use

 1. Open Scilab
 2. Source the SIVP toolbox (from Scilab menu bar)
 3. Source the `canny_filter.sci` file
 4. Run the `doCannyFilter()` function in the prompt
 5. Pick the image file to be filtered

**Warning: do not pick a large image. This being an experiment, it relies on the implementation of custom low-level mathematical operations, such as the convolution operation. To speed things up, you may use Scilab native conv2() function instead.**

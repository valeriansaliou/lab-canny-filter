// =====================
// === CANNY FILTER ====
// == Valerian Saliou ==
// =====================

// -----
// UTILS
// -----

function imageMatrix=loadImage(path, isRGB)
    if isRGB == 0 then
        imageMatrix=double(imread(path));
    else
        imageMatrix=double(rgb2gray(imread(path)));
    end
endfunction

function displayImage(imageMatrix)
    imshow(uint8(imageMatrix));
endfunction

function image=writeImage(imageMatrix, fileName)
    image=imwrite(imageMatrix, fileName);
endfunction


// ------
// FILTER
// ------

function doCannyFilter()
    path_file_name=uiputfile(["*.png"; "*.jpg"; "*.jpeg"; "*.bmp"]);
    original_image=loadImage(path_file_name, 1);

    // Apply filters
    gaussian_image = applyGaussianFilter(original_image);
    [sobel_image, angles] = applySobelFilter(gaussian_image);
    maxima_image = applyNonMaximaFilter(sobel_image, angles);
    hysteresis_image = applyHysteresisFilter(maxima_image, sobel_image, angles);

    disp("All processed");

    // Concatenate images
    image = cat(2, original_image, gaussian_image, sobel_image, maxima_image, hysteresis_image);

    displayImage(image);
endfunction

function matrix=gaussianMatrix(n, _kernel)
    x = n;
    y = n;

    [h1, h2] = meshgrid(-(x-1)/2:(x-1)/2, -(y-1)/2:(y-1)/2);
    hg = exp(- (h1.^2+h2.^2) / (2*_kernel^2));

    matrix = hg ./ sum(hg(:));
endfunction

function image=applyGaussianFilter(image)
    // Generate Gaussian matrix
    mask = gaussianMatrix(5, 0.86);

    disp("==> BEGIN Gaussian mask <==");
    disp("Processing...");
    disp(mask);
    image = applyMaskFilter(image, mask);
    disp("==> END Gaussian mask <==");
endfunction

function filtered_image=applyMaskFilter(image, mask)
    // Proceed: convolution(mask, image)

    image_size = size(image)

    // Uncomment this to bypass custom convolutional filter
    //filtered_image = conv2(image, mask);
    //filtered_image = resize_matrix(filtered_image, image_size(1), image_size(2));

    filtered_image = zeros(image_size(1), image_size(2));

    // Mask parameters
    mask_size = size(mask);

    mask_height = mask_size(1);
    mask_width = mask_size(2);

    mask_center_y = ceil(mask_height / 2);
    mask_center_x = ceil(mask_width / 2);

    image_size = size(image)
    image_height = image_size(1)
    image_width = image_size(2)

    // Apply mask filter
    for y = 1:image_height
        for x = 1:image_width
            // Apply mask to pixel(x,y)
            if (size(image(y, x)) > 1)
                masked_pixel = [0, 0, 0];
            else
                masked_pixel = 0;
            end

            for my = 1:mask_height
                // Y-inset value?
                y_add = 0;

                if (my < mask_center_y)
                    y_add = -1 * (mask_center_y - my);
                elseif (my > mask_center_y)
                    y_add = my - mask_center_y;
                end

                for mx = 1:mask_width
                    // X-inset value?
                    x_add = 0;

                    if (mx < mask_center_x)
                        x_add = -1 * (mask_center_x - mx);
                    elseif (mx > mask_center_x)
                        x_add = mx - mask_center_x;
                    end

                    y_new = y + y_add;
                    x_new = x + x_add;

                    // Pixel exists in image? (not out-of-borders)
                    if (y_new <= image_height & y_new > 0)
                        if (x_new <= image_width & x_new > 0)
                          if (size(image(y, x)) > 1)
                                for cord = 1:3
                                    masked_pixel(cord) = masked_pixel(cord) + (mask(my, mx) * image(y_new, x_new, cord));
                                end
                            else
                                masked_pixel = masked_pixel + (mask(my, mx) * image(y_new, x_new));
                            end
                        end
                    end
                end
            end

            if (size(image(y, x)) > 1)
              filtered_image(y, x, 1) = round(masked_pixel(1));
              filtered_image(y, x, 2) = round(masked_pixel(2));
              filtered_image(y, x, 3) = round(masked_pixel(3));
            else
              filtered_image(y, x) = round(masked_pixel);
            end
        end
    end
endfunction

function merged_image=mergeGradientMatrixes(image, gradient_x, gradient_y)
    gradient_size = size(image);
    merged_image = zeros(gradient_size(1), gradient_size(2));

    for y = 1:gradient_size(1)
        for x = 1:gradient_size(2)
            merged_image(y, x) = sqrt(gradient_x(y, x)^2 + gradient_y(y, x)^2);
        end
    end
endfunction

function angles=processAngleMatrix(image, gradient_x, gradient_y)
    gradient_size = size(image);
    angles = zeros(gradient_size(1), gradient_size(2));

    for y = 1:gradient_size(1)
        for x = 1:gradient_size(2)
            J_x = gradient_x(y, x);
            J_y = gradient_y(y, x);

            if (J_x ~= 0 & J_y ~= 0)
                angle_value = atand((-1 * J_x) / J_y);

                if (angle_value < 0)
                    angle_value = angle_value + 180;
                end

                // Cap angle to directional values (0, 45, 90, 135)
                // eg:  horizontal      (0)
                //      vertical        (90)
                //      diagonal right  (45)
                //      diagonal left   (135)
                if (angle_value <= 22.5)
                    angle_value = 0;
                elseif (angle_value >= 157.5)
                    angle_value = 0;
                elseif (angle_value < 67.5)
                    angle_value = 45;
                elseif (angle_value < 112.5)
                    angle_value = 90;
                elseif (angle_value < 157.5)
                    angle_value = 135;
                end

                angles(y, x) = angle_value;
            end
        end
    end
endfunction

function [pixel_1_coords, pixel_2_coords]=pixelCoordsGradient(pixel, angle, y, x)
    pixel_1_coords = [-1, -1];
    pixel_2_coords = [-1, -1];

    if (angle == 0)
        // X-Gradient
        pixel_1_coord_x = (x - 1);
        pixel_2_coord_x = (x + 1);

        if (pixel_1_coord_x >= 1)
            pixel_1_coords = [y, pixel_1_coord_x];
        end

        if (pixel_2_coord_x <= image_width)
            pixel_2_coords = [y, pixel_2_coord_x];
        end
    elseif (angle == 90)
        // Y-Gradient
        pixel_1_coord_y = (y - 1);
        pixel_2_coord_y = (y + 1);

        if (pixel_1_coord_y >= 1)
            pixel_1_coords = [pixel_1_coord_y, x];
        end

        if (pixel_2_coord_y <= image_height)
            pixel_2_coords = [pixel_2_coord_y, x];
        end
    elseif (angle == 45)
        // Diagonal-Right-Gradient
        pixel_1_coord_x = (x + 1);
        pixel_1_coord_y = (y - 1);

        pixel_2_coord_x = (x - 1);
        pixel_2_coord_y = (y + 1);

        if (pixel_1_coord_x <= image_width & pixel_1_coord_y >= 1)
            pixel_1_coords = [pixel_1_coord_y, pixel_1_coord_x];
        end

        if (pixel_2_coord_x >= 1 & pixel_2_coord_y <= image_height)
            pixel_2_coords = [pixel_2_coord_y, pixel_2_coord_x];
        end
    elseif (angle == 135)
        // Diagonal-Left-Gradient
        pixel_1_coord_x = (x - 1);
        pixel_1_coord_y = (y - 1);

        pixel_2_coord_x = (x + 1);
        pixel_2_coord_y = (y + 1);

        if (pixel_1_coord_x >= 1 & pixel_1_coord_y >= 1)
            pixel_1_coords = [pixel_1_coord_y, pixel_1_coord_x];
        end

        if (pixel_2_coord_x <= image_width & pixel_2_coord_y <= image_height)
            pixel_2_coords = [pixel_2_coord_y, pixel_2_coord_x];
        end
    end
endfunction

function perpendicular=getPerpendicularAngle(angle)
    perpendicular = angle + 90;

    // Circle back to inner trigo. circle (edge cases)
    // eg:  180deg  >  0deg  (+ others; default)
    //      225deg  >  45deg
    if (perpendicular >= 180)
        if (perpendicular == 225)
            perpendicular = 45;
        else
            perpendicular = 0;
        end
    end
endfunction

function maximum=processMaximumMatrix(image, angles)
    image_size = size(image);

    maximum = zeros(image_size(1), image_size(2));

    image_height = image_size(1);
    image_width = image_size(2);

    for y = 1:image_height
        for x = 1:image_width
            pixel = image(y, x);
            angle = angles(y, x);

            // Acquire surrounding gradient pixel coordinates
            [pixel_1_coords, pixel_2_coords] = pixelCoordsGradient(pixel, angle, y, x);

            // Check bounds & acquire pixel values
            pixel_1 = 0;
            pixel_2 = 0;

            if (pixel_1_coords(1) ~= -1 & pixel_1_coords(2) ~= -1)
               pixel_1 = image(pixel_1_coords(1), pixel_1_coords(2));
            end

            if (pixel_2_coords(1) ~= -1 & pixel_2_coords(2) ~= -1)
               pixel_2 = image(pixel_2_coords(1), pixel_2_coords(2));
            end

            maximum_pixel = max(pixel_1, pixel, pixel_2);

            // Maximum? Can insert.
            if (maximum_pixel == pixel)
                maximum(y, x) = pixel;
            end
        end
    end
endfunction

function hysteresis=processHysteresisMatrix(image, gradients, angles)
    // Configuration
    T_line = 255;
    T_h = processGradientModulus(gradients);
    T_i = (0.5 * T_h);

    image_size = size(image);

    hysteresis = zeros(image_size(1), image_size(2));

    image_height = image_size(1);
    image_width = image_size(2);

    for y = 1:image_height
        for x = 1:image_width
            pixel = image(y, x);

            if (pixel > T_h)
                // Upper-limit pixel
                hysteresis(y, x) = T_line;
            elseif (pixel > T_i)
                // Lower limit pixel, check if connected
                // Get perpendicular angle direction
                perpendicular = getPerpendicularAngle(angles(y, x));

                // Acquire surrounding gradient pixel coordinates
                [pixel_1_coords, pixel_2_coords] = pixelCoordsGradient(pixel, perpendicular, y, x);

                // Check bounds & acquire pixel values
                pixel_1 = 0;
                pixel_2 = 0;

                if (pixel_1_coords(1) ~= -1 & pixel_1_coords(2) ~= -1)
                   pixel_1 = image(pixel_1_coords(1), pixel_1_coords(2));
                end

                if (pixel_2_coords(1) ~= -1 & pixel_2_coords(2) ~= -1)
                   pixel_2 = image(pixel_2_coords(1), pixel_2_coords(2));
                end

                // Current pixel connected outline points?
                if (pixel_1 > T_h | pixel_2 > T_h)
                    // Consider this one as outline, too
                    hysteresis(y, x) = T_line;
                end
            end
        end
    end
endfunction

function [image, angles]=applySobelFilter(image)
    // Sobel filter
    mask_gx = [-1, 0, 1]';
    mask_gy = [-1, 0, 1];

    disp("==> BEGIN Sobel-Gx mask <==");
    disp("Processing...");
    disp(mask_gx);
    gradient_x = applyMaskFilter(image, mask_gx);
    disp(size(gradient_x))
    disp("==> END Sobel-Gx mask <==");


    disp("==> BEGIN Sobel-Gy mask <==");
    disp("Processing...");
    disp(mask_gy);
    gradient_y = applyMaskFilter(image, mask_gy);
    disp(size(gradient_y));
    disp("==> END Sobel-Gy mask <==");

    disp("==> BEGIN Gradient angles <==");
    disp("Processing...");
    angles = processAngleMatrix(image, gradient_x, gradient_y);
    disp(size(image));
    disp("==> END Gradient angles <==");

    image = mergeGradientMatrixes(image, gradient_x, gradient_y);
endfunction

function image=applyNonMaximaFilter(image, angles)
    disp("==> BEGIN Non maxima <==");
    disp("Processing...");
    image = processMaximumMatrix(image, angles);
    disp("==> END Non maximas <==");
endfunction

function image=applyHysteresisFilter(image, gradients, angles)
    disp("==> BEGIN Hysteresis <==");
    disp("Processing...");
    image = processHysteresisMatrix(image, gradients, angles);
    disp("==> END Hysteresis <==");
endfunction

function histogram=processHistogram(repartition_matrix)
    // 256 possible values, initiate all values to zero
    histogram = zeros(1, 256);

    repartition_matrix_size = size(repartition_matrix)

    for y = 1:repartition_matrix_size(1)
        for x = 1:repartition_matrix_size(2)
            index = floor(repartition_matrix(y, x)) + 1

            histogram(index) = histogram(index) + 1
        end
    end
endfunction

function percentile=processPercentile(repartition_matrix, percent)
    // First step: flatten matrix (convert to a straight vector)
    flat_matrix = repartition_matrix(:);

    // Second step: sort matrix
    flat_matrix = gsort(flat_matrix, 'g', 'i')

    // Third step: extract nth percentile
    index_percentile = floor((percent / 100) * size(flat_matrix));

    percentile = flat_matrix(index_percentile(1));
endfunction

function modulus=processGradientModulus(gradients)
    percentile = processPercentile(gradients, 84);

    modulus = percentile(1);
endfunction

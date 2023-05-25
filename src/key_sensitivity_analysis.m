img1 = imread('encrypted_image_1.png');
img2 = imread('encrypted_image_2.png');

[rows, cols, ~] = size(img1);

match_count = 0;

for row = 1:rows
    for col = 1:cols
        if img1(row, col) == img2(row, col)
            match_count = match_count + 1;
        end
    end
end

percent_match = match_count/(rows * cols)*100;
fprintf('Percentage of matching pixels: %.4f%%\n', percent_match);
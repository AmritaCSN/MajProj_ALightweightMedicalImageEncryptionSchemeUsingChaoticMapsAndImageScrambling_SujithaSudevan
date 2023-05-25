img = imread('encrypted_image_1.png');

histogram = imhist(img);
pdf = histogram / numel(img);

H = -sum(pdf .* log2(pdf));

disp(['Shannon''s entropy of encrypted image: ' num2str(H)]);
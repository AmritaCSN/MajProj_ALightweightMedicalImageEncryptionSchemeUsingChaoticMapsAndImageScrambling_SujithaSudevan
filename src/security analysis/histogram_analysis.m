plain_image = imread('plain_image_1.png');
encrypted_image = imread('encrypted_image_1.png');
decrypted_image = imread('decrypted_image_1.png');

subplot(1,3,1);
imhist(plain_image);
title('Original Image Histogram');

subplot(1,3,2);
imhist(encrypted_image);
title('Encrypted Image Histogram');

subplot(1,3,3);
imhist(decrypted_image);
title('Decrypted Image Histogram');
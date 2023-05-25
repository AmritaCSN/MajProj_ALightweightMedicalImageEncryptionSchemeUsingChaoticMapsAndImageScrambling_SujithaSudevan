clc;
clear;

%read the image
img = imread('xray_chest_1024.png');

%size of the image
[height, width, depth] = size(img);

%selecting random pixel location
x = randi(width);
y = randi(height);

out = img;

%assigning random value to pixel
out(y, x, :) = randi([0 255], [1 1 depth]);

% Display the modified image
imwrite(out,'xray_chest_1024_changed.png')

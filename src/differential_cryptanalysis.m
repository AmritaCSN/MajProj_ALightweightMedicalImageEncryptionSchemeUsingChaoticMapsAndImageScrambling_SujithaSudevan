%reading 2 encrypted images that are encrypted from 2 plain images that differ just by 1 bit
E1 = imread('encrypted_image_1.png');
E2 = imread('encrypted_image_2.png');

E1 = double(E1);
E2 = double(E2);

[rows,cols] = size(E1);
D = (E1~=E2);

NPCR = sum(sum(D))/(rows*cols)*100;
fprintf('NPCR: %f%%.',NPCR);
fprintf('\n');
UACI = sum(sum(abs(E1-E2)))/(255*rows*cols)*100;
fprintf('UACI: %f%%.',UACI);
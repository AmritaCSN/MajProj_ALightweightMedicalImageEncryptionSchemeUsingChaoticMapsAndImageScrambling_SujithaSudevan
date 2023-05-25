I = imread('encrypted_image_1.png');

C=0;
I = uint8(I);
ang0 = graycomatrix(I, 'offset', [0,1], 'Symmetric', false, 'NumLevels', 256);
ang45 = graycomatrix(I, 'offset', [-1,1], 'Symmetric', false, 'NumLevels', 256);
ang90 = graycomatrix(I, 'offset', [-1,0], 'Symmetric', false, 'NumLevels', 256);
ang135 = graycomatrix(I, 'offset', [-1,-1], 'Symmetric', false, 'NumLevels', 256);

s0 = sum(sum(ang0));
s45 = sum(sum(ang45));
s90 = sum(sum(ang90));
s135 = sum(sum(ang135));

ang0 = ang0/s0;
ang45 = ang45/s45;
ang90 = ang90/s90;
ang135 = ang135/s135;

mC = (ang0+ang45+ang90+ang135)/4;
[m,n] = size(mC);
for i=1:m
    for j=1:n
        C = C+mC(i,j)*(i-j)^2;
    end
end

fprintf('Contrast Value for Encrypted Image: %.2f%%\n',C);
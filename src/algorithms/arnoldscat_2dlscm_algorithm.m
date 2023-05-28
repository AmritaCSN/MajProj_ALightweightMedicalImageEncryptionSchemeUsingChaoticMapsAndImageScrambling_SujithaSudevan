clc;
clear;

tic;                %start timer
mem1 = memory;      %memory before

% reading image file
img = imread('256 x 256 (1).png');

%getting size of image
[rows,cols] = size(img);

%encryption key
x = 0.18; y = 0.91; cp = 0.142; n = 256; pi = 3.14;

%chaotic matrix generation
p = x; q = y;
Cx = zeros(rows,cols);
Cy = zeros(rows,cols);
for i = 1:rows
    for j = 1:cols
        p = round(sin(pi*(4*cp*p*(1-p)+(1-cp)*sin(pi*q))), 2);
        q = round(sin(pi*(4*cp*q*(1-q)+(1-cp)*sin(pi*p))), 2);
        Cx(i,j) = abs(p);
        Cy(i,j) = abs(q);
    end
end
Cx = mod(floor(Cx.*2^32),256);
Cy = mod(floor(Cy.*2^32),256);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%encryption%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%arnolds cat map transformations
enc = ArnoldCatMapEncryption(img,1);
enc = uint8(enc);

%permutation
O = zeros(rows,cols);
Z = zeros(1,cols);
[~,idx1] = sort(Cx,1);

for i = 1:rows
    for j = 1:cols
        Z(1,j) = Cx(idx1(i,j),j);
    end
    [~,idx2] = sort(Z,2);
    for k = 1:cols
        O(idx1(i,k),k) = enc(idx1(i,idx2(1,k)),idx2(1,k));
    end
end
permutation = O;

%diffusion
T = permutation;
D = zeros(size(T),'like',T);
D_ = zeros(size(T),'like',T);
T = double(T);
D = double(D);
D_ = double(D_);
%column diffusion
D(1,:) = mod(T(1,:) + T(rows,:) + T(rows-1,:) + Cx(1,:),256);
D(2,:) = mod(T(2,:) + D(1,:) + T(rows,:) + Cx(2,:),256);
for i = 3:rows
    D(i,:) = mod(T(i,:) + D(i-1,:) + D(i-2,:) + Cx(i,:),256);
end
%row diffusion
D_(:,1) = mod(D(:,1) + D(:,cols) + D(:,cols-1) + Cx(:,1),256);
D_(:,2) = mod(D(:,2) + D_(:,1) + D(:,cols) + Cx(:,2),256);
for i = 3:cols
    D_(:,i) = mod(D(:,i) + D_(:,i-1) + D_(:,i-2) + Cx(:,i),256);
end
diffusion = double(D_);
diffusion = uint8(diffusion);

%2nd round permutation
O1 = zeros(rows,cols);
Z1 = zeros(1,cols);
[~,idx11] = sort(Cx,1);
for i = 1:rows
    for j = 1:cols
        Z1(1,j) = Cx(idx11(i,j),j);
    end
    [~,idx21] = sort(Z1,2);
    for k = 1:cols
        O1(idx11(i,k),k) = diffusion(idx11(i,idx21(1,k)),idx21(1,k));
    end
end
permutation1 = O1;

%2nd round diffusion
T1 = permutation1;
D1 = zeros(size(T1),'like',T1);
D1_ = zeros(size(T1),'like',T1);
T1 = double(T1);
D1 = double(D1);
D1_ = double(D1_);
%column diffusion
D1(1,:) = mod(T1(1,:) + T1(rows,:) + T1(rows-1,:) + Cx(1,:),256);
D1(2,:) = mod(T1(2,:) + D1(1,:) + T1(rows,:) + Cx(2,:),256);
for i = 3:rows
    D1(i,:) = mod(T1(i,:) + D1(i-1,:) + D1(i-2,:) + Cx(i,:),256);
end
%row diffusion
D1_(:,1) = mod(D1(:,1) + D1(:,cols) + D1(:,cols-1) + Cx(:,1),256);
D1_(:,2) = mod(D1(:,2) + D1_(:,1) + D1(:,cols) + Cx(:,2),256);
for i = 3:cols
    D1_(:,i) = mod(D1(:,i) + D1_(:,i-1) + D1_(:,i-2) + Cx(:,i),256);
end
enc_img = double(D1_);
enc_img = uint8(enc_img);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%decryption%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%getting the size
[rows_,cols_] = size(enc_img);

%decryption key
x = 0.18; y = 0.91; cp = 0.142; n = 256; pi = 3.14;

%chaotic matrix generation
p = x; q = y;
Cx = zeros(rows,cols);
Cy = zeros(rows,cols);
for i = 1:rows
    for j = 1:cols
        p = round(sin(pi*(4*cp*p*(1-p)+(1-cp)*sin(pi*q))), 2);
        q = round(sin(pi*(4*cp*q*(1-q)+(1-cp)*sin(pi*p))), 2);
        Cx(i,j) = abs(p);
        Cy(i,j) = abs(q);
    end
end
Cx = mod(floor(Cx.*2^32),256);
Cy = mod(floor(Cy.*2^32),256);

%diffusion
A = enc_img;
P = zeros(size(A),'like',A);
P_ = zeros(size(A),'like',A);
A = double(A);
P = double(P);
P_ = double(P_);
%row diffusion
for i = cols:-1:3
    P(:,i) = mod(A(:,i)-A(:,i-1)-A(:,i-2)-Cx(:,i),256);
end
P(:,2) = mod(A(:,2)-P(:,cols)-A(:,1)-Cx(:,2),256);
P(:,1) = mod(A(:,1)-P(:,cols)-P(:,cols-1)-Cx(:,1),256);
%column diffusion
for i = rows_:-1:3
    P_(i,:) = mod(P(i,:)-P(i-1,:)-P(i-2,:)-Cx(i,:),256);
end
P_(2,:) = mod(P(2,:)-P_(rows_,:)-P(1,:)-Cx(2,:),256);
P_(1,:) = mod(P(1,:)-P_(rows_,:)-P_(rows_-1,:)-Cx(1,:),256);
permutation_ = double(P_);
permutation_ = uint8(permutation_);

%permutation
O_ = zeros(rows_,cols_);
Z_ = zeros(1,cols_);
[~,idx1_] = sort(Cx,1);
for i = 1:rows_
    for j = 1:cols_
        Z_(1,j) = Cx(idx1_(i,j),j);
    end
    [~,idx2_] = sort(Z_,2);
    for k = 1:cols_
        O_(idx1_(i,idx2_(1,k)),idx2_(1,k)) = permutation_(idx1_(i,k),k);
    end
end
circ_col_shift_ = O_;
circ_col_shift_ = uint8(circ_col_shift_);

%2nd round diffusion
A1 = circ_col_shift_;
P1 = zeros(size(A1),'like',A1);
P1_ = zeros(size(A1),'like',A1);
A1 = double(A1);
P1 = double(P1);
P1_ = double(P1_);
%row diffusion
for i = cols:-1:3
    P1(:,i) = mod(A1(:,i)-A1(:,i-1)-A1(:,i-2)-Cx(:,i),256);
end
P1(:,2) = mod(A1(:,2)-P1(:,cols_)-A1(:,1)-Cx(:,2),256);
P1(:,1) = mod(A1(:,1)-P1(:,cols_)-P1(:,cols_-1)-Cx(:,1),256);
%column diffusion
for i = rows:-1:3
    P1_(i,:) = mod(P1(i,:)-P1(i-1,:)-P1(i-2,:)-Cx(i,:),256);
end
P1_(2,:) = mod(P1(2,:)-P1_(rows_,:)-P1(1,:)-Cx(2,:),256);
P1_(1,:) = mod(P1(1,:)-P1_(rows_,:)-P1_(rows_-1,:)-Cx(1,:),256);
permutation1_ = double(P1_);
permutation1_ = uint8(permutation1_);

%2nd round permutation
O1_ = zeros(rows_,cols_);
Z1_ = zeros(1,cols_);
[~,idx11_] = sort(Cx,1);
for i = 1:rows_
    for j = 1:cols_
        Z1_(1,j) = Cx(idx11_(i,j),j);
    end
    [~,idx21_] = sort(Z1_,2);
    for k = 1:cols_
        O1_(idx11_(i,idx21_(1,k)),idx21_(1,k)) = permutation1_(idx11_(i,k),k);
    end
end

circ_col_shift1_ = O1_;
circ_col_shift1_ = uint8(circ_col_shift1_);

%arnolds cat map transformation
dec_img = ArnoldCatMapDecryption(circ_col_shift1_,1);
dec_img = uint8(dec_img);
dec_img = rot90(dec_img);
dec_img = rot90(dec_img);
dec_img = fliplr(dec_img);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imwrite(uint8(img),'plain_image_2.png');
imwrite(uint8(enc_img),'encrypted_image_2.png');
imwrite(uint8(dec_img),'decrypted_image_2.png');

%plain image
subplot(1, 3, 1);
imshow(img);
title('Original Image');

%encrypted image
subplot(1,3,2);
imshow(uint8(enc_img));
title({'Encrypted Image','using Arnolds Cat & 2D-LSCM Algorithm'});

%decrypted image
subplot(1,3,3);
imshow(uint8(dec_img));
title({'Decrypted Image','using Arnolds Cat & 2D-LSCM Algorithm'});

t = toc;            %stop timer
mem2 = memory;      %memory after

mem_used = (mem2.MemUsedMATLAB - mem1.MemUsedMATLAB)/(1024^2);

disp(['End-to-End Execution Time: ', num2str(t), ' seconds']);
disp(['End-to-End Memory Utilized: ', num2str(mem_used), ' MB']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [enc] = ArnoldCatMapEncryption(image, iterations)
    [m, n] = size(image);
    enc = zeros(m, n);
    for k = 1:iterations
        for i = 1:m
            for j = 1:n
                enc(mod(i+j-1,m)+1, mod(i+2*j-2,n)+1) = image(i, j);
            end
        end
        image = enc;
    end
end

function [dec_img] = ArnoldCatMapDecryption(enc, iterations)
    [m, n] = size(enc);
    dec_img = zeros(m, n);
    for k = 1:iterations
        for i = 1:m
            for j = 1:n
                dec_img(mod(j-2*i+2*m,m)+1, mod(j-i+n,n)+1) = enc(i, j);
            end
        end
        enc = dec_img;
    end
end
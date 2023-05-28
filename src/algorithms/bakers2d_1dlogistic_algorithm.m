clc;
clear;

tic;                %start timer
mem1 = memory;      %memory before

% reading image file
img = imread('256 x 256 (1).png');

%getting size of image
[rows,cols] = size(img);

%encryption rounds
T1 = 4; T2 = 4; R = 1;
T = 1;

%seed values
x1 = 0.3245;
x2 = 0.6467;
y1 = 0.3125;
y2 = 0.3813;
z1 = 0.1045;
z2 = 0.2509;
mu1 = 0.2506;
mu2 = 0.4526;
p1 = 4;
p2 = 4;
x0 = 0.4; % initial value of x
y0 = 0.4; % initial value of y
p = 0.4; % system parameter

X1k = zeros(rows,cols);
X2k = zeros(rows,cols);
Y1k = zeros(rows,cols);
Y2k = zeros(rows,cols);

%key generation
for i = 1:rows
    for j = 1:cols
        if x0 <= p
            X1k(i, j) = x0/p;
            Y1k(i, j) = p * y0;
        else
            X1k(i, j) = (x0 - p)/(1 - p);
            Y1k(i, j) = (1 - p) * y0 + 1 - p;
        end
        if y0 <= p
            X2k(i, j) = y0/p;
            Y2k(i, j) = p * x0;
        else
            X2k(i, j) = (y0 - p)/(1 - p);
            Y2k(i, j) = (1 - p) * x0 + 1 - p;
        end
        x0 = X1k(i, j);
        y0 = Y1k(i, j);
    end
end

range_min = 3.5699456;
range_max = 4;

%obtain sequence X1(n*T1) from X1k
X1nT1 = X1k(1:T1:(end/(2^(log2(rows)-2))));

%apply mathematical operation to bring values in range
X1nT1(X1nT1 < range_min) = range_min;
X1nT1(X1nT1 > range_max) = range_max;

%add modified values to μ1(n)
mu1n = x1 + (X1nT1);
X2nT2 = X2k(1:T2:(end/(2^(log2(rows)-2))));
X2nT2(X2nT2 < range_min) = range_min;
X2nT2(X2nT2 > range_max) = range_max;
mu2n = x2 + (X2nT2);

range_min = 0;
range_max = 1;

Y1nT1 = Y1k(1:T1:(end/(2^(log2(rows)-2))));
Y1nT1(Y1nT1 < range_min) = range_min;
Y1nT1(Y1nT1 > range_max) = range_max;
a1n = y1 + (Y1nT1);

Y2nT2 = Y2k(1:T2:(end/(2^(log2(rows)-2))));
Y2nT2(Y2nT2 < range_min) = range_min;
Y2nT2(Y2nT2 > range_max) = range_max;
a2n = y2 + (Y2nT2);

MxN = rows*cols;
mi = randi(numel(mu1n));
ai = randi(numel(a1n));
mu1 = mu1n(mi);
a1 = a1n(ai);

%initial values for Z1 sequence
z1 = 0.1;
z1_seq = zeros(1, MxN);
for k = 1:MxN
    if mod(k-1, T1) == 0
        mu_k = mu1;
    end
    if mod(k-1, T2) == 0
        z1 = a1*z1*(1-z1);
    end
    z1 = mu_k*z1*(1-z1);
    z1_seq(k) = z1;
end

%initial values for Z2 sequence
z2 = 0.2;
z2_seq = zeros(1, MxN);
for k = 1:MxN
    if mod(k-1, T1) == 0
        mu_k = mu1;
    end
    if mod(k-1, T2) == 0
        z2 = a1*z2*(1-z2);
    end
    z2 = mu_k*z2*(1-z2);
    z2_seq(k) = z2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%encryption%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = zeros(rows,cols);
[sorted_z1_seq, idx] = sort(z1_seq);

Fs = zeros(1, MxN);
for i = 1:rows
    for j = 1:cols
        E((i-1)*cols+j) = img(i,j);
    end
end
for k = 1:MxN
    Fk = E(idx(k));
    Fs(k) = Fk;
end

Fs = uint8(Fs);
E = reshape(E,[rows,cols]);

%generating uniformly distributed sequence {Gk}
G = mod(floor(z2_seq*10^15),256);
G = uint8(G);

%initializing substitution matrix
H = zeros(size(Fs));

%performing XOR to generate substitution matrix
for t = 1:numel(Fs)
    i = floor((t-1)/cols) + 1;
    j = t - (i-1)*cols;
    H = bitxor(G,Fs);
end
H1 = bitxor(G,Fs);

H = reshape(H1,[rows,cols]);
enc_img = uint8(H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%decryption%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%seed values
x1 = 0.3245;
x2 = 0.6467;
y1 = 0.3125;
y2 = 0.3813;
z1 = 0.1045;
z2 = 0.2509;
mu1 = 0.2506;
mu2 = 0.4526;
p1 = 4;
p2 = 4;
x0 = 0.4; % initial value of x
y0 = 0.4; % initial value of y
p = 0.4; % system parameter

X1k = zeros(rows,cols);
X2k = zeros(rows,cols);
Y1k = zeros(rows,cols);
Y2k = zeros(rows,cols);

%key generation
for i = 1:rows
    for j = 1:cols
        if x0 <= p
            X1k(i, j) = x0/p;
            Y1k(i, j) = p * y0;
        else
            X1k(i, j) = (x0 - p)/(1 - p);
            Y1k(i, j) = (1 - p) * y0 + 1 - p;
        end
        if y0 <= p
            X2k(i, j) = y0/p;
            Y2k(i, j) = p * x0;
        else
            X2k(i, j) = (y0 - p)/(1 - p);
            Y2k(i, j) = (1 - p) * x0 + 1 - p;
        end
        x0 = X1k(i, j);
        y0 = Y1k(i, j);
    end
end

range_min = 3.5699456;
range_max = 4;

%obtain sequence X1(n*T1) from X1k
X1nT1 = X1k(1:T1:(end/(2^(log2(rows)-2))));

%apply mathematical operation to bring values in range
X1nT1(X1nT1 < range_min) = range_min;
X1nT1(X1nT1 > range_max) = range_max;

%add modified values to μ1(n)
mu1n = x1 + (X1nT1);
X2nT2 = X2k(1:T2:(end/(2^(log2(rows)-2))));
X2nT2(X2nT2 < range_min) = range_min;
X2nT2(X2nT2 > range_max) = range_max;
mu2n = x2 + (X2nT2);

range_min = 0;
range_max = 1;

Y1nT1 = Y1k(1:T1:(end/(2^(log2(rows)-2))));
Y1nT1(Y1nT1 < range_min) = range_min;
Y1nT1(Y1nT1 > range_max) = range_max;
a1n = y1 + (Y1nT1);

Y2nT2 = Y2k(1:T2:(end/(2^(log2(rows)-2))));
Y2nT2(Y2nT2 < range_min) = range_min;
Y2nT2(Y2nT2 > range_max) = range_max;
a2n = y2 + (Y2nT2);

MxN = rows*cols;
mi = randi(numel(mu1n));
ai = randi(numel(a1n));
mu1 = mu1n(mi);
a1 = a1n(ai);

%initial values for Z1 sequence
z1 = 0.1;
z1_seq = zeros(1, MxN);
for k = 1:MxN
    if mod(k-1, T1) == 0
        mu_k = mu1;
    end
    if mod(k-1, T2) == 0
        z1 = a1*z1*(1-z1);
    end
    z1 = mu_k*z1*(1-z1);
    z1_seq(k) = z1;
end

%initial values for Z2 sequence
z2 = 0.2;
z2_seq = zeros(1, MxN);
for k = 1:MxN
    if mod(k-1, T1) == 0
        mu_k = mu1;
    end
    if mod(k-1, T2) == 0
        z2 = a1*z2*(1-z2);
    end
    z2 = mu_k*z2*(1-z2);
    z2_seq(k) = z2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%decryption%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_ = reshape(enc_img,[1,rows*cols]);

G_ = mod(floor(z2_seq*10^15),256);
G_ = uint8(G_);

Fs_ = bitxor(H_,G_);

E_ = zeros(rows,cols);
for k_ = 1:MxN
    Fk_ = Fs_(k_);
    E_(idx(k_)) = Fk_;
end

dec_img = transpose(E_);
dec_img = uint8(dec_img);

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
title({'Encrypted Image','using 2D-Bakers+1D-Logistic Algorithm'});

%decrypted image
subplot(1,3,3);
imshow(uint8(dec_img));
title({'Decrypted Image','using 2D-Bakers+1D-Logistic Algorithm'});

t = toc;            %stop timer
mem2 = memory;      %memory after

mem_used = (mem2.MemUsedMATLAB - mem1.MemUsedMATLAB)/(1024^2);

disp(['End-to-End Execution Time: ', num2str(t), ' seconds']);
disp(['End-to-End Memory Utilized: ', num2str(mem_used), ' MB']);

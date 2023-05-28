clc;
clear;

tic;                %start timer
mem1 = memory;      %memory before

%READING IMAGE FILE (PLAIN IMAGE)
img = imread('1024 x 1024 (1).png');

%GET SIZE OF IMAGE (ROWS AND COLUMNS)
[rows,cols] = size(img);

%INITIAL IMAGE SETUP
%dividing the image into blocks of 4x4
block_size = [16 16];
blocks = mat2cell(img, block_size(1)*ones(1,rows/16), block_size(2)*ones(1,cols/16));
[rows_,cols_] = size(blocks);

%PHASE 1: KEY GENERATION

%PHASE 1A: KEY GENERATION FOR BAKER'S MAP
x = 0.16; y = 0.91; n = 256; cp = 0.142;
xkey = zeros(n,1); ykey = zeros(n,1);
for i = 1:n
    xkey(i) = x;
    ykey(i) = y;
    condition = x<=0.5;
    x = condition.*(2*x) + (~condition).*(2*x-1);
    y = condition.*(cp*y) + (~condition).*(cp*y+0.5);
    x = x + rand()*0.001;
    y = y + rand()*0.001;
end
key1 = int64(mod(xkey*(10^3),rows));
key2 = int64(mod(ykey*(10^6),cols));
key1 = transpose(key1);
key2 = transpose(key2);

%KEY TRANSFORMATION
key1t = reshape(key1,[16,16]);
key1t = uint8(key1t);

%PHASE 1B: KEY GENERATION FOR IMAGE SCRAMBLING
key1is = randi([1,rows],1,rows);

%PHASE 1C: CHAOTIC MATRIX GENERATION FOR 2D-LSCM MAP
p = x; q = y;
Cx = zeros(rows,cols);
Cy = zeros(rows,cols);
for i = 1:rows
    for j = 1:cols
        p = round(sin(pi*(4*cp*p*(1-p)+(1-cp)*sin(pi*q))),2);
        q = round(sin(pi*(4*cp*q*(1-q)+(1-cp)*sin(pi*p))),2);
        Cx(i,j) = abs(p);
        Cy(i,j) = abs(q);
    end
end

Cx = mod(floor(Cx.*2^32),256);
Cy = mod(floor(Cy.*2^32),256);

%%%%%%%%%%%%%%%%%%%%%%%%ENCRYPTION STARTS HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PHASE 2: BAKER'S MAP ENCRYPTION

%PHASE 2A: KEY-BINDING
for i = 1:numel(blocks)
    key_binded{i} = bitxor(blocks{i}, key1t);
end
key_binded = cell2mat(reshape(key_binded,[rows_,cols_]));

%PHASE 2B: BAKER'S TRANSFORMATION (STRETCHING & STACKING)
bakers = key_binded;
for i = 1:4
    s1 = bakers(:,1:rows/2,:);
    s2 = bakers(:,rows/2+1:end,:);
    s1_stretched = reshape(s1,[floor(rows/2),cols]);
    s2_stretched = reshape(s2,[floor(rows/2),cols]);
    bakers = cat(1,s1_stretched,s2_stretched);
end

%PHASE 2C: SHIFTING PIXELS ROW-WISE
shifted_image = [];
for i = 1:rows
    shifted_image = cat(1,shifted_image,circshift(bakers(i,:),-i));
end

%PHASE 3: IMAGE SCRAMBLING

%PHASE 3A: ROW SWAPPING ACCORDING TO KEY
row_swap = shifted_image;
for i = 1:numel(key1is)
    row_swap([1,key1is(i)],:) = row_swap([key1is(i),1],:);
end

%PHASE 3B: COLUMN SWAPPING ACCORDING TO KEY
col_swap = row_swap;
for i = 1:numel(key1is)
    col_swap(:,[1,key1is(i)]) = col_swap(:,[key1is(i),1]);
end

%PHASE 3C: CIRCULAR ROW SHIFT ACCORDING TO KEY
circ_row_shift = zeros(size(col_swap), 'like', col_swap);
for i = 1:size(col_swap,1)
    circ_row_shift(i,:) = circshift(col_swap(i,:), [0, key1is(i)]);
end

%PHASE 3D: CIRCULAR COLUMN SHIFT ACCORDING TO KEY
circ_col_shift = zeros(size(circ_row_shift), 'like', circ_row_shift);
for i = 1:size(circ_row_shift,2)
    circ_col_shift(:,i) = circshift(circ_row_shift(:,i), [key1is(i), 0]);
end

%PHASE 4: 2D-LSCM

%PHASE 4A: 2D-LSCM PERMUTATION ROUND 1
O = zeros(rows,cols); 
Z = zeros(1,cols);
[~,idx1] = sort(Cx,1); 
for i = 1:rows
    for j = 1:cols
        Z(1,j) = Cx(idx1(i,j),j);
    end
    [~,idx2] = sort(Z,2);
    for k = 1:cols
        O(idx1(i,k),k) = circ_col_shift(idx1(i,idx2(1,k)),idx2(1,k));
    end
end
permutation = O;

%PHASE 4B: 2D-LSCM DIFFUSION ROUND 1
T = permutation;
D = zeros(size(T),'like',T);
D_ = zeros(size(T),'like',T);
T = double(T);
D = double(D);
D_ = double(D_);

%COLUMN DIFFUSION ROUND 1
D(1,:) = mod(T(1,:) + T(rows,:) + T(rows-1,:) + Cx(1,:),256);
D(2,:) = mod(T(2,:) + D(1,:) + T(rows,:) + Cx(2,:),256);
for i = 3:rows
    D(i,:) = mod(T(i,:) + D(i-1,:) + D(i-2,:) + Cx(i,:),256);
end

%ROW DIFFUSION ROUND 1
D_(:,1) = mod(D(:,1) + D(:,cols) + D(:,cols-1) + Cx(:,1),256);
D_(:,2) = mod(D(:,2) + D_(:,1) + D(:,cols) + Cx(:,2),256);
for i = 3:cols
    D_(:,i) = mod(D(:,i) + D_(:,i-1) + D_(:,i-2) + Cx(:,i),256);
end

diffusion = uint8(D_);

%PHASE 4C: 2D-LSCM PERMUTATION ROUND 2
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

%PHASE 4D: 2D-LSCM DIFFUSION ROUND 2
T1 = permutation1;
D1 = zeros(size(T1),'like',T1);
D1_ = zeros(size(T1),'like',T1);
T1 = double(T1);
D1 = double(D1);
D1_ = double(D1_);

%COLUMN DIFFUSION ROUND 2
D1(1,:) = mod(T1(1,:) + T1(rows,:) + T1(rows-1,:) + Cx(1,:),256);
D1(2,:) = mod(T1(2,:) + D1(1,:) + T1(rows,:) + Cx(2,:),256);
for i = 3:rows
    D1(i,:) = mod(T1(i,:) + D1(i-1,:) + D1(i-2,:) + Cx(i,:),256);
end

%ROW DIFFUSION ROUND 2
D1_(:,1) = mod(D1(:,1) + D1(:,cols) + D1(:,cols-1) + Cx(:,1),256);
D1_(:,2) = mod(D1(:,2) + D1_(:,1) + D1(:,cols) + Cx(:,2),256);
for i = 3:cols
    D1_(:,i) = mod(D1(:,i) + D1_(:,i-1) + D1_(:,i-2) + Cx(:,i),256);
end

enc_img = uint8(D1_);

%%%%%%%%%%%%%%%%%%%%%%DECRYPTION STARTS HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PHASE 4D: 2D-LSCM DIFFUSION ROUND 2
%ROW DIFFUSION ROUND 2
A = enc_img;
P = zeros(size(A),'like',A);
P_ = zeros(size(A),'like',A);
A = double(A);
P = double(P);
P_ = double(P_);

for i = cols:-1:3
    P(:,i) = mod(A(:,i)-A(:,i-1)-A(:,i-2)-Cx(:,i),256);
end
P(:,2) = mod(A(:,2)-P(:,cols)-A(:,1)-Cx(:,2),256);
P(:,1) = mod(A(:,1)-P(:,cols)-P(:,cols-1)-Cx(:,1),256);

%COLUMN DIFFUSION ROUND 2
for i = rows:-1:3
    P_(i,:) = mod(P(i,:)-P(i-1,:)-P(i-2,:)-Cx(i,:),256);
end
P_(2,:) = mod(P(2,:)-P_(rows,:)-P(1,:)-Cx(2,:),256);
P_(1,:) = mod(P(1,:)-P_(rows,:)-P_(rows-1,:)-Cx(1,:),256);

permutation_ = uint8(P_);

%PHASE 4C: 2D-LSCM PERMUATION ROUND 2
O_ = zeros(rows,cols);
Z_ = zeros(1,cols);
[~,idx1_] = sort(Cx,1);
for i = 1:rows
    for j = 1:cols
        Z_(1,j) = Cx(idx1_(i,j),j);
    end
    [~,idx2_] = sort(Z_,2);
    for k = 1:cols
        O_(idx1_(i,idx2_(1,k)),idx2_(1,k)) = permutation_(idx1_(i,k),k);
    end
end
circ_col_shift_ = O_;
circ_col_shift_ = uint8(circ_col_shift_);

%PHASE 4B: 2D-LSCM ROUND 1
%ROW DIFFUSION ROUND 1
A1 = circ_col_shift_;
P1 = zeros(size(A1),'like',A1);
P1_ = zeros(size(A1),'like',A1);
A1 = double(A1);
P1 = double(P1);
P1_ = double(P1_);
for i = cols:-1:3
    P1(:,i) = mod(A1(:,i)-A1(:,i-1)-A1(:,i-2)-Cx(:,i),256);
end
P1(:,2) = mod(A1(:,2)-P1(:,cols)-A1(:,1)-Cx(:,2),256);
P1(:,1) = mod(A1(:,1)-P1(:,cols)-P1(:,cols-1)-Cx(:,1),256);

%COLUMN DIFFUSION ROUND 1
for i = rows:-1:3
    P1_(i,:) = mod(P1(i,:)-P1(i-1,:)-P1(i-2,:)-Cx(i,:),256);
end
P1_(2,:) = mod(P1(2,:)-P1_(rows,:)-P1(1,:)-Cx(2,:),256);
P1_(1,:) = mod(P1(1,:)-P1_(rows,:)-P1_(rows-1,:)-Cx(1,:),256);

permutation1_ = double(P1_);
permutation1_ = uint8(permutation1_);

%PHASE 4A: 2D-LSCM PERMUTATION ROUND 1
O1_ = zeros(rows,cols);
Z1_ = zeros(1,cols);
[~,idx11_] = sort(Cx,1);
for i = 1:rows
    for j = 1:cols
        Z1_(1,j) = Cx(idx11_(i,j),j);
    end
    [~,idx21_] = sort(Z1_,2);
    for k = 1:cols
        O1_(idx11_(i,idx21_(1,k)),idx21_(1,k)) = permutation1_(idx11_(i,k),k);
    end
end
circ_col_shift1_ = O1_;
circ_col_shift1_ = uint8(circ_col_shift1_);

%PHASE 3: IMAGE SCRAMBLING

%PHASE 3D: CIRCULAR COLUMN SHIFT ACCORDING TO KEY
circ_row_shift_ = zeros(size(circ_col_shift1_), 'like', circ_col_shift1_);
for i = size(circ_col_shift1_,2):-1:1
    circ_row_shift_(:,i) = circshift(circ_col_shift1_(:,i), [-key1is(i), 0]);
end

%PHASE 3C: CIRCULAR ROW SHIFT ACCORDING TO KEY
col_swap_ = zeros(size(circ_row_shift_), 'like', circ_row_shift_);
for i = size(circ_row_shift_,1):-1:1
    col_swap_(i,:) = circshift(circ_row_shift_(i,:), [0, -key1is(i)]);
end

%PHASE 3B: COLUMN SWAPPING ACCORDING TO KEY
row_swap_ = col_swap_;
for i = numel(key1is):-1:1
    row_swap_(:,[1,key1is(i)]) = row_swap_(:,[key1is(i),1]);
end

%PHASE 3A: ROW SWAPPING ACCORDING TO KEY
shifted_image_ = row_swap_;
for i = numel(key1is):-1:1
    shifted_image_([1,key1is(i)],:) = shifted_image_([key1is(i),1],:);
end

%PHASE 2: BAKER'S MAP

%CHANGING DECRYPTION KEY BY 1 BIT (ONLY FOR KEY SENSITIVITY ANALYSIS)
x = 0.36; y = 0.91; n = 256; cp = 0.142;
xkey = zeros(n,1); ykey = zeros(n,1);
for i = 1:n
    xkey(i) = x;
    ykey(i) = y;
    condition = x<=0.5;
    x = condition.*(2*x) + (~condition).*(2*x-1);
    y = condition.*(cp*y) + (~condition).*(cp*y+0.5);
    x = x + rand()*0.001;
    y = y + rand()*0.001;
end

%PHASE 2C: SHIFITNG PIXELS ROW WISE
bakers_ = [];
for i = 1:rows
    bakers_ = cat(1,bakers_,circshift(shifted_image_(i,:),i));
end

%PHASE 2B: BAKER'S TRANSFORMATION (STRETCHING & STACKING)
key_binded_ = bakers_;
for i = 1:4
    r1 = key_binded_(1:cols/2,:,:);
    r2 = key_binded_(cols/2+1:end,:,:);
    r1_compressed = reshape(r1,[rows,floor(cols/2)]);
    r2_compressed = reshape(r2,[rows,floor(cols/2)]);
    key_binded_ = cat(2,r1_compressed,r2_compressed);
end
key_binded_ = uint8(key_binded_);

%PHASE 2A: KEY BINDING
blocks_ = mat2cell(key_binded_, block_size(1)*ones(1,rows/16), block_size(2)*ones(1,cols/16));
[rows_,cols_] = size(blocks_);
for i = 1:numel(blocks)
    og_image{i} = bitxor(blocks_{i}, key1t);
end
dec_img = cell2mat(reshape(og_image,[rows_,cols_]));

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
title({'Encrypted Image','using Proposed Algorithm'});

%decrypted image
subplot(1,3,3);
imshow(uint8(dec_img));
title({'Decrypted Image','using Proposed Algorithm'});

t = toc;            %stop timer
mem2 = memory;      %memory after

mem_used = (mem2.MemUsedMATLAB - mem1.MemUsedMATLAB)/(1024^2);

disp(['End-to-End Execution Time: ', num2str(t), ' seconds']);
disp(['End-to-End Memory Utilized: ', num2str(mem_used), ' MB']);
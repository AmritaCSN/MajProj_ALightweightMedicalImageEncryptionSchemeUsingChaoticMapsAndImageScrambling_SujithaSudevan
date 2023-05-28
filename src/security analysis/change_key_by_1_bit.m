clc;
clear;

%original seed value for key
seed = 0.18;

%converting the seed to binary
seed_binary = dec2bin(typecast(single(seed),'uint32'),32);

%changing random 1 bit in binary representation
bit_to_change = randi(numel(seed_binary));
seed_binary(bit_to_change) = num2str(~str2double(seed_binary(bit_to_change)));

%converting binary back to floating point
new_seed = typecast(uint32(bin2dec(seed_binary)),'single');

disp(['Original Seed: ', num2str(seed)]);
disp(['Modified Seed: ', num2str(new_seed)]);
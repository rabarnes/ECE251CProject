%Test loop with various wavelet types and lengths across a few diferent
%images
%note: use periodization dwt mode
clear;
clc;
close all;
s = rng(1);

imdata1 = phantom("Modified Shepp-Logan", 256);
% brain = imread('MRA_cartoid.jpeg');
% brain = imread('brain3.jpeg');
brain = imread('brain.jpg');
imdata2 = imresize(brain, [256 256]);
imdata2 = im2double(imdata2(:,:,1)); %black and white image, all layers same, conv to double [0 1]
% imdata2 = phantom('Modified Shepp-Logan', 256);
iter_length = 50;
% wnames = ["haar", "db2", "db8", "db16", "sym2", "sym8", "sym16", "bior2.2", "bior3.3", "bior4.4", ...
%     "fk4", "fk8", "fk18"];
wnames = ["haar", "bior1.1", "db2"];
% wnames = ["haar"];
w_mse = zeros(1,length(wnames));
w_psnr = zeros(1,length(wnames));

threshold_weight = 0.016;

for i=1:length(wnames)
    [mse, peak_snr_dif] = WaveletType_Loop(i, wnames(i),iter_length,imdata2, threshold_weight);
    w_mse(i) = mse;
    w_psnr(i) = peak_snr_dif;
end

figure;
X = categorical({'DWT', 'DT CDWT', 'DD DT CDWT'});
X = reordercats(X,{'DWT', 'DT CDWT', 'DD DT CDWT'});
bar(X, w_mse)
title("MSE vs DWT, Dual-Tree DWT, and Double Density DT DWT - Brain")


figure;
bar(X, w_psnr)
title('Peak SNR Increase vs Varying Wavelet Type and Filter Length - Brain');


% figure;
% X = categorical({'haar', 'db2,8,16', 'sym2,8,16', 'bior2.2,3.3,4.4', 'fk4,8,18'});
% X = reordercats(X,{'haar', 'db2,8,16', 'sym2,8,16', 'bior2.2,3.3,4.4', 'fk4,8,18'});
% bar(X,[w_psnr(1) 0 0; w_psnr(2:4); w_psnr(5:7); w_psnr(8:10); w_psnr(11:13)]);
% title('Peak SNR Increase vs Varying Wavelet Type and Filter Length - Brain');
% 
% figure;
% bar(X,[w_mse(1) 0 0; w_mse(2:4); w_mse(5:7); w_mse(8:10); w_mse(11:13)]);
% title("Mean Square Error vs varying Wavelet Type and Filter Length - Brain")





close all;
clear;
clc;

%%
%%% MRI Image
% imdata = phantom('Modified Shepp-Logan', 256);
%load("brain3Mod.mat");
%imdata = brain3Mod;
MRAneck = double(imread('MRA_cartoid.jpeg'));
% %imdata = imresize(MRAneck(:,:,1), [512 512]); % for 512x512 image

imdata = MRAneck(:,:,3)./255;
imdata = imdata(1:end-2,1:end-2);
% 
% minV = min(min(abs(imdata)));
% maxV = max(max(abs(imdata)));
% figure; imshow(abs(imdata), [minV maxV]); title('Original Image');



%%% Create Probability Density Function - PDF
PDF = create_PDF(imdata);
% figure; imshow(PDF); title("PDF");

rows = size(imdata,1);

%%% Create Sampling Mask
% mask = make_mask(rows, 4);
% mask = make_gauss_mask(rows, 1);
[mask, percent] = make_spiral_mask(rows, 1);
percent
figure; imshow(mask); title("Mask Image");


%%

%%% For Gaussian Settings
% iter_length = 75;
% threshold_weight = 0.014; 

% %%% For Cartesian Settings
% iter_length = 100;
% threshold_weight = 0.016; %for cart
iter_length = 20;
threshold_weight = 0.50; %for cart

% [imdata, im_og, im_final, MSE, peaksnr] = dddt_CS_loop(imdata,PDF, mask, iter_length, threshold_weight);
[imdata, im_og, im_final, MSE, peaksnr] = basic_CS_loop(imdata,PDF, mask, iter_length, threshold_weight);

minV = min(min(abs(imdata)));
maxV = max(max(abs(imdata)));

min_og = min(min(abs(im_og)));
max_og = max(max(abs(im_og)));

min_final = min(min(abs(im_final)));
max_final = max(max(abs(im_final)));

figure;
subplot(1,3,1); imshow(abs(imdata), [minV maxV]);
title("Orignial Image");
subplot(1,3,2); imshow(abs(im_og), [min_og max_og]);
title("Sparse Image");
subplot(1,3,3); imshow(abs(im_final), [min_final max_final]);
title("Final image");


figure;
subplot(1,2,1); plot(MSE);
subplot(1,2,2); plot(peaksnr);


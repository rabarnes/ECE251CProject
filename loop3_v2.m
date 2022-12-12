close all;
clear;
clc;

%%
%%% MRI Image
% imdata = phantom('Modified Shepp-Logan', 256);

load("brain3Mod.mat");
imdata = brain3Mod;
imdata = imresize(imdata, [256 256]);

% MRAneck = double(imread('MRA_cartoid.jpeg'));
% imdata = MRAneck(:,:,3)./255;
% imdata = imresize(imdata, [256 256]);


%%% imdata = imdata(1:end-2,1:end-2);
%%imdata = imresize(MRAneck(:,:,1), [512 512]); % for 512x512 image

% minV = min(min(abs(imdata)));
% maxV = max(max(abs(imdata)));
% figure; imshow(abs(imdata), [minV maxV]); title('Original Image');



%%% Create Probability Density Function - PDF
PDF = create_PDF(imdata);
% figure; imshow(PDF); title("PDF");

rows = size(imdata,1);

%%

%%% Create Sampling Mask
% [mask, percent] = make_mask(rows, 4); %4
% [mask, percent] = make_gauss_mask(rows, 4000);
[mask, percent] = make_gauss_mask(rows, 0.2);
% [mask, percent] = make_spiral_mask(rows, 1);
fprintf("Percent: %f \n", percent);
% figure; imshow(mask); title("Mask Image");


%%


%%% Best Cartesian Settings (basic_CS_loop)
% iter_length = 150;
% threshold_weight = 0.0119; %for Shepp-Logan - 77%
% iter_length = 150;
% threshold_weight = 0.0117; %for brain3Mod - 77%
% iter_length = 150;
% threshold_weight = 0.0127; %for MRA - 74%

%%% Best Cartesian Settings (dddt_CS_loop)
% iter_length = 150;
% threshold_weight = 0.0012; %for Shepp-Logan - 77%
% iter_length = 150;
% threshold_weight = 0.0234; %for brain3Mod - 77%
% iter_length = 150;
% threshold_weight = 0.00001; %for MRA - 74%


%%% Best spiral Settings (basic_CS_loop)
% iter_length = 150;
% threshold_weight = 0.004812; %for Shepp-Logan - 74%
% iter_length = 42;
% threshold_weight = 0.0240; %for brain3Mod - 74%
% iter_length = 48;
% threshold_weight = 0.04794; %for MRA - 74%
% iter_length = 37;
% threshold_weight = 0.195; %for MRA - 26.4% (628x628)


%%% Best spiral Settings (dddt_CS_loop)
% iter_length = 150;
% threshold_weight = 0.00188; %for Shepp-Logan - 74%
% iter_length = 35;
% threshold_weight = 0.01298; %for brain3Mod - 74%
% iter_length = 150;
% threshold_weight = 0.01905; %for MRA - 74%

%%% Best Gaussian Settings (basic_CS_loop)
% iter_length = 150;
% threshold_weight = 0.0187; %for Shepp-Logan - 77%
% iter_length = 150;
% threshold_weight = 0.0002; %for brain3Mod - 77%
% iter_length = 150;
% threshold_weight = 0.020; %for MRA - 74%


%%% Best Gaussian Settings (dddt_CS_loop)
% iter_length = 150;
% threshold_weight = 0.0076; %for Shepp-Logan - 77%
% iter_length = 150;
% threshold_weight = 0.0002; %for brain3Mod - 77%
% iter_length = 150;
% threshold_weight = 0.0391; %for MRA - 74%



wnames = 'haar';
% wnames = "";


% [imdata, im_og, im_final, MSE, peaksnr] = dddt_CS_loop(imdata,PDF, mask, iter_length, threshold_weight);
[imdata, im_og, im_final, MSE, peaksnr] = basic_CS_loop(imdata,PDF, mask, wnames, iter_length, threshold_weight);

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

% 
% figure;
% subplot(1,2,1); plot(MSE);
% subplot(1,2,2); plot(peaksnr);

%fprintf("Mean Square Error =  %e \n\n", MSE(end));

fprintf("Mean Square Error =  %e \n", min(MSE));
[M,indx] = min(MSE);
fprintf("MSE Index =  %i \n", indx);
fprintf("PSNR =  %e \n\n", max(peaksnr));

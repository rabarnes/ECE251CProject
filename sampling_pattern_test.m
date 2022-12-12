close all;
clear;
clc;

%%
%%% MRI Image
im_phantom = phantom('Modified Shepp-Logan', 256);

load("brain3Mod.mat");
im_brain = brain3Mod;
im_brain = imresize(im_brain, [256 256]);

im_MRA = double(imread('MRA_cartoid.jpeg'));
im_MRA = im_MRA(:,:,3)./255;
im_MRA = imresize(im_MRA, [256 256]);


%%% Create Probability Density Function - PDF
PDF = create_PDF(im_phantom);


%%% Create Sampling Mask
mask_size = 256;
[cartesian_mask, per1] = make_mask(mask_size, 4);
fprintf("Cartesian Mask Percent: %f \n", per1);
[gauss_mask, per2] = make_gauss_mask(mask_size, 4000);
fprintf("Gaussian Mask Percent: %f \n", per2);
[spiral_mask, per3] = make_spiral_mask(mask_size, 1);
fprintf("Spiral Mask Percent: %f \n", per3);




%%
iter_length = 150;

Shepp_Logan_basic_PSNR = zeros(1,3);

threshold_weight = 0.0119; %for Shepp-Logan
[~, ~, ~, ~, peaksnr] = basic_CS_loop(im_phantom,PDF, cartesian_mask, 'haar', iter_length, threshold_weight);
Shepp_Logan_basic_PSNR(1) = max(peaksnr);

threshold_weight = 0.004812; %for Shepp-Logan
[~, ~, ~, ~, peaksnr] = basic_CS_loop(im_phantom,PDF, spiral_mask, 'haar', iter_length, threshold_weight);
Shepp_Logan_basic_PSNR(2) = max(peaksnr);

threshold_weight = 0.0187; %for Shepp-Logan
[~, ~, ~, ~, peaksnr] = basic_CS_loop(im_phantom,PDF, gauss_mask, 'haar', iter_length, threshold_weight);
Shepp_Logan_basic_PSNR(3) = max(peaksnr);



brain_basic_PSNR = zeros(1,3);
threshold_weight = 0.0117; %for brain3Mod
[~, ~, ~, ~, peaksnr] = basic_CS_loop(im_brain,PDF, cartesian_mask, 'haar', iter_length, threshold_weight);
brain_basic_PSNR(1) = max(peaksnr);

threshold_weight = 0.0240; %for brain3Mod
[~, ~, ~, ~, peaksnr] = basic_CS_loop(im_brain,PDF, spiral_mask, 'haar', iter_length, threshold_weight);
brain_basic_PSNR(2) = max(peaksnr);

threshold_weight = 0.0002; %for brain3Mod
[~, ~, ~, ~, peaksnr] = basic_CS_loop(im_brain,PDF, gauss_mask, 'haar', iter_length, threshold_weight);
brain_basic_PSNR(3) = max(peaksnr);


MRA_basic_PSNR = zeros(1,3);
threshold_weight = 0.0127; %for MRA
[~, ~, ~, ~, peaksnr] = basic_CS_loop(im_MRA,PDF, cartesian_mask, 'haar', iter_length, threshold_weight);
MRA_basic_PSNR(1) = max(peaksnr);

threshold_weight = 0.04794; %for MRA
[~, ~, ~, ~, peaksnr] = basic_CS_loop(im_MRA,PDF, spiral_mask, 'haar', iter_length, threshold_weight);
MRA_basic_PSNR(2) = max(peaksnr);

threshold_weight = 0.020; %for MRA
[~, ~, ~, ~, peaksnr] = basic_CS_loop(im_MRA,PDF, gauss_mask, 'haar', iter_length, threshold_weight);
MRA_basic_PSNR(3) = max(peaksnr);


X = categorical({'Shepp Logan', 'brain', 'MRA'});
X = reordercats(X,{'Shepp Logan', 'brain', 'MRA'});

figure;
bar(X,[Shepp_Logan_basic_PSNR; brain_basic_PSNR; MRA_basic_PSNR]);
title("Peak SNR with Haar for Different Sampling Patterns");
exit
legend("Cartesian", "Spiral", "Gaussian");


%%
iter_length = 150;

Shepp_Logan_dddt_PSNR = zeros(1,3);
threshold_weight = 0.0012; %for Shepp-Logan
[~, ~, ~, ~, peaksnr] = dddt_CS_loop(im_phantom,PDF, cartesian_mask, iter_length, threshold_weight);
Shepp_Logan_dddt_PSNR(1) = max(peaksnr);

threshold_weight = 0.00188; %for Shepp-Logan
[~, ~, ~, ~, peaksnr] = dddt_CS_loop(im_phantom,PDF, spiral_mask, iter_length, threshold_weight);
Shepp_Logan_dddt_PSNR(2) = max(peaksnr);

threshold_weight = 0.0076; %for Shepp-Logan
[~, ~, ~, ~, peaksnr] = dddt_CS_loop(im_phantom,PDF, gauss_mask, iter_length, threshold_weight);
Shepp_Logan_dddt_PSNR(3) = max(peaksnr);


brain_dddt_PSNR = zeros(1,3);
threshold_weight = 0.0234; %for brain3Mod
[~, ~, ~, ~, peaksnr] = dddt_CS_loop(im_brain,PDF, cartesian_mask, iter_length, threshold_weight);
brain_dddt_PSNR(1) = max(peaksnr);

threshold_weight = 0.01298; %for brain3Mod
[~, ~, ~, ~, peaksnr] = dddt_CS_loop(im_brain,PDF, spiral_mask, iter_length, threshold_weight);
brain_dddt_PSNR(2) = max(peaksnr);

threshold_weight = 0.0002; %for brain3Mod
[~, ~, ~, ~, peaksnr] = dddt_CS_loop(im_brain,PDF, gauss_mask, iter_length, threshold_weight);
brain_dddt_PSNR(3) = max(peaksnr);


MRA_dddt_PSNR = zeros(1,3);
threshold_weight = 0.00001; %for MRA
[~, ~, ~, ~, peaksnr] = dddt_CS_loop(im_MRA,PDF, cartesian_mask, iter_length, threshold_weight);
MRA_dddt_PSNR(1) = max(peaksnr);

threshold_weight = 0.01905; %for MRA
[~, ~, ~, ~, peaksnr] = dddt_CS_loop(im_MRA,PDF, spiral_mask, iter_length, threshold_weight);
MRA_dddt_PSNR(2) = max(peaksnr);

threshold_weight = 0.0391; %for MRA
[~, ~, ~, ~, peaksnr] = dddt_CS_loop(im_MRA,PDF, gauss_mask, iter_length, threshold_weight);
MRA_dddt_PSNR(3) = max(peaksnr);


X = categorical({'Shepp Logan', 'brain', 'MRA'});
X = reordercats(X,{'Shepp Logan', 'brain', 'MRA'});

figure;
bar(X,[Shepp_Logan_dddt_PSNR; brain_dddt_PSNR; MRA_dddt_PSNR]);
title("Peak SNR with DDDT for Different Sampling Patterns");

legend("Cartesian", "Spiral", "Gaussian");



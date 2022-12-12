close all;
clear;
clc;

%%


im_phantom = phantom('Modified Shepp-Logan', 256);

load("brain3Mod.mat");
im_brain = brain3Mod;
im_brain = imresize(im_brain, [256 256]);

im_MRA = double(imread('MRA_cartoid.jpeg'));
im_MRA = im_MRA(:,:,3)./255;
im_MRA = imresize(im_MRA, [256 256]);



%% Create Sampling Mask
[Cartesian_mask, per1] = make_mask(256, 2); %4
[gauss_mask, per2] = make_gauss_mask(256, 1000);
[spiral_mask, per3] = make_spiral_mask(256, 6);

per1
per2
per3
%%
figure; 
subplot(3,2,1); imshow(im_phantom); title("Shepp-Logan Phantom");
subplot(3,2,3); imshow(im_brain); title("Brain Image");
subplot(3,2,5); imshow(im_MRA); title("MRA Image");
subplot(3,2,2); imshow(Cartesian_mask); title("Cartesian Mask");
subplot(3,2,4); imshow(spiral_mask); title("Spiral Mask");
subplot(3,2,6); imshow(gauss_mask); title("Gaussian Mask");

%%
figure; 
% subplot(3,2,1); imshow(im_phantom); title("Shepp-Logan Phantom");
% subplot(3,2,3); imshow(im_brain); title("Brain Image");
% subplot(3,2,5); imshow(im_MRA); title("MRA Image");
subplot(1,3,1); imshow(Cartesian_mask); title("Cartesian Mask");
subplot(1,3,2); imshow(spiral_mask); title("Spiral Mask");
subplot(1,3,3); imshow(gauss_mask); title("Gaussian Mask");


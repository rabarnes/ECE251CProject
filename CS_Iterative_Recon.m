% Tritai Nguyen
% iterative reconstruction for undersampled images
% parameters to tune: 
%   percent sample:     % of total pixels to sample gaussian and apply mask
%   n_iters:            # of loop iterations in iterative reconstruction
%   debug:              0/1 - display/don't display images in loop
%   N_x:                size of image
%   [LoD,HiD,LoR, HiR]: type of wavelet filters
%   thresh_coef:        top % of wavelet coef to keep

clear
close all
clc

N_x = 200;           % size of image
x   = phantom(N_x); % initialize phantom image
figure(); imshow(x,[0 1]); title('Original Image');

%% step 1 generate measured k space
% creating sampling mask and applying to image
dim_im = size(x);
num_pixels = dim_im(1) * dim_im(2);
percent_sample = 0.3;
num_subsample = round(num_pixels*percent_sample);
% random gaussian
u = [0 0];
sigma = [1 0; 0 1];
% generate subsamples of multivariate gaussian
R = mvnrnd(u, sigma, round(num_subsample));
% scale all values of R to indices of image
a = 1;
b = min(dim_im(1), dim_im(2)); % TODO fix this for rectangular images
min_R = min(R(:));
max_R = max(R(:));
R = ((b-a)*(R - min_R)) / (max_R-min_R) + a;
R = round(R);
% create mask by setting mask indices to 1
mask = zeros([dim_im(1) dim_im(2)]);
for i=1:num_subsample
    mask(R(i,1), R(i,2)) = 1;
end
% figure 
% imshow(gaussian_mask);
% title('Gaussian Mask');

debug = 0;
n_iters = 10;
thresh_coef = 0.3;

k = 0;
beta = 0.9;
T = 0.0095;
e = 0.001;
rho = 1;
z_k_new = x;
t_k_new = 1;
x_k_new = x;

[LoD,HiD,LoR, HiR] = wfilters('haar');

% apply fft to orig image, rescale from 0-1 and apply mask
y = fftshift(fft2(x));
y = ((y - min(y(:)))) / (max(y(:))-min(y(:)));
y = y.*mask;

figure;
y_im = abs(ifft2(ifftshift(y)));
y_im = ((y_im - min(y_im(:)))) / (max(y_im(:))-min(y_im(:)));
imshow(abs(y_im));
title("Received image");

while (k < n_iters) 
    x_k_old = x_k_new;
    t_k_old = t_k_new;
    z_k_old = z_k_new;
    
    % apply F_u'*(F_u*z_k-y)
    % grad in direction of minima
    % F_u' is ifft2, F_u is fft
    z_temp = fftshift(fft2(z_k_old));
    z_temp = ((z_temp - min(z_temp(:)))) / (max(z_temp(:))-min(z_temp(:)));
    z_temp = mask.*z_temp;
    
    grad_f = ifft2(ifftshift( (z_temp - y) ));
    x_grad = z_k_old - rho*(grad_f);
    
    [cA,cH,cV,cD] = dwt2(x_grad,LoD,HiD,'mode','symh');
    if (debug == 1)
        figure; imshow(abs(cA))
        title("wavelet transform of image before thresh")
    end

    %% thresholding
    % keep top thresh_coef values in wavelet transform, else set to 0
    m = sort(abs([cA(:) ; cH(:) ; cV(:) ; cD(:)]), 'descend');
    ndx = floor(length(m)*thresh_coef);
    thresh = m(ndx);

    cA_denoise = cA .* (abs(cA) > thresh);
    cD_denoise = cD .* (abs(cD) > thresh);
    cH_denoise = cH .* (abs(cH) > thresh);
    cV_denoise = cV .* (abs(cV) > thresh);

    if (debug == 1)
        figure; imshow(abs(cA_denoise))
        title("wavelet transform of image after thresh")
    end

    %% inv W trans
    % get denoised image
    x_k_new = idwt2(cA_denoise, cH_denoise, cV_denoise, cD_denoise, LoR, HiR);
    
    if (T > e) 
        T = T*beta;
        t_k_new = (1+sqrt(1+4*t_k_old^2)) / 2;
        z_k_new = x_k_new + ((t_k_old-1)/t_k_new)*(x_k_new - x_k_old);
    end
    k = k+1;

end

figure; imshow(abs(x_k_new));
title("denoised image");


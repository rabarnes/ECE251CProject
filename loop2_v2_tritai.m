close all;
clc;
clear;


%%% MRI Image
imdata = phantom('Modified Shepp-Logan', 256);
% figure; imshow(abs(imdata)); title('Shepp-Logan Image');


%%% Create Probability Density Function - PDF
PDF = create_PDF(imdata);
% figure; imshow(PDF); title("PDF");


rows = size(imdata,1);
cols = size(imdata,2);

% ft_weight = 1/(norm(imdata,"fro"));
ft_weight = 1/sqrt(size(imdata,1)*size(imdata,2));
% ft_weight = 1/(max(imdata,[],'all'));
F_imdata = fftshift(fft2(imdata).*ft_weight);


%%% Create Sampling Mask
mask = make_mask(rows, 4);
% mask = make_gauss_mask(rows, 1);
% figure; imshow(mask); title("Mask Image");


%%

F_imdata_sp = F_imdata.*mask;
im_sp = ifft2(ifftshift(F_imdata_sp))./ft_weight;
im_sp = im_sp./(abs(max(im_sp,[],'all')));
im_og = im_sp;
im_sp = im_sp./PDF;


%%% For Gaussian Settings
% iter_length = 75;
% threshold_weight = 0.014; 

%%% For Cartesian Settings
iter_length = 30;
threshold_weight = 0.016; %for cart


im_final = zeros(rows,cols);
mean_squared_error = ones(1, iter_length);
peaksnr = ones(1, iter_length);

%%% For iterative DWT
cA_hist = struct;
cH_hist = struct;
cV_hist = struct;
cD_hist = struct;
n_dwt = 7;
dwt_fname = 'bior1.1';


for n = 1:iter_length
    %wavelet transform
    cA = im_sp;
    for i = 1:n_dwt
        %[cA,cH,cV,cD] = dwt2(cA, LoD, HiD,'mode','symh');
        [cA,cH,cV,cD] = dwt2(cA, dwt_fname);
        cA_hist(i).data = cA;
        cH_hist(i).data = cH;
        cV_hist(i).data = cV;
        cD_hist(i).data = cD;
    end
    
    % only use threshold from first DWT
    thresh = abs(threshold_weight*max([cA_hist(1).data(:) ; cH_hist(1).data(:) ; cV_hist(1).data(:) ; cD_hist(1).data(:)], [], 'all'));

    for i = n_dwt:-1:1
        % thresholding
        % keep top thresh_coef values in wavelet transform, else set to 0
        cA_denoise = threshold2(cA_hist(i).data, thresh);
        cD_denoise = threshold2(cD_hist(i).data, thresh);
        cH_denoise = threshold2(cH_hist(i).data, thresh);
        cV_denoise = threshold2(cV_hist(i).data, thresh);
         
        % inv W trans
        % get denoised image
%         im_sp_th = idwt2(cA_denoise, cH_denoise, cV_denoise, cD_denoise, LoR, HiR);
        im_sp_th = idwt2(cA_denoise, cH_denoise, cV_denoise, cD_denoise,  dwt_fname);

        if (i > 1) 
            cA_hist(i-1).data = im_sp_th;
        end
    end

    %go to k space and downsample
    F_sp_th = fftshift(fft2(im_sp_th).*ft_weight);
    F_sp_th_masked = F_sp_th.*(1-mask);


    %find err in k space and ifft to get difference image (compare to original)
    F_err = F_sp_th_masked + F_imdata_sp; %Ax-y error between sampled original and wavelet recon in k space
    delta_im = ifft2(ifftshift(F_err))./ft_weight;

    %add diff im to original in image space
    new_im = im_sp + delta_im;
    im_sp = new_im./abs(max(new_im, [],'all'));
    im_final = im_final + delta_im;
    im_final = im_final./abs(max(im_final,[],'all'));
    im_final(im_final < 10^(-10)) = 0;
    mean_squared_error(n) = immse(imdata, abs(im_final));
    peaksnr(n) = psnr(imdata, abs(im_final));
end

minV = min(min(abs(imdata)));
maxV = max(max(abs(imdata)));

figure;
subplot(1,3,1); imshow(abs(imdata), [minV maxV]);
title("Orignial Image");
subplot(1,3,2); imshow(abs(im_og), [minV maxV]);
title("Sparse Image");
subplot(1,3,3); imshow(abs(im_final), [minV maxV]);
title("Final image");


figure;
subplot(1,2,1); plot(mean_squared_error);
subplot(1,2,2); plot(peaksnr);


function pdf = create_PDF(input_image)
    rows = size(input_image, 1);
    cols = size(input_image, 2);
    [x, y] = meshgrid(linspace(-1, 1, cols), linspace(-1, 1, rows));
    partial_sampling_factor = 0.5;
    
    r = sqrt(x.^2 + y.^2);
    r = r/max(abs(r(:)));
    pdf = (1-r).^2;
    
    minval = 0;
    maxval = 1;
    val = 0.5;
    PCTG = floor(partial_sampling_factor*rows*cols);
    
    while true
	    val = minval/2 + maxval/2;
	    pdf = (1-r).^2 + val; 
        pdf(pdf>1) = 1; 
    
	    N = floor(sum(pdf(:)));
	    if N > PCTG % infeasible
		    maxval = val;
	    end
	    if N < PCTG % not optimal
		    minval = val;
	    end
	    if N == PCTG % optimal
		    break;
	    end
    end

end

function mask = make_mask(size, m)
    
    mask = ones(size,size);
    mask(1:m:round(4*size/9),:) = zeros(length(1:m:round(4*size/9)),size);
    mask(round(5*size/9):m:end,:) = zeros(length(round(5*size/9):m:size),size);


end

function th = threshold2(A, thresh)
    th = (abs(A) > thresh).*(A.*(abs(A)-thresh)./(abs(A)+eps));
end

function [mse, peak_snr] = WaveletType_Loop(wname,iter_length,imdata)
close all

% brain = imread('brain.jpg');
% imdata = im2double(brain(:,:,1)); %black and white image, all layers same, conv to double [0 1]
% imdata = phantom('Modified Shepp-Logan', 256);

% size_dif = size(imdata,1)-size(imdata,2);
% imdata = [zeros(size(imdata,1), size_dif/2 +1) imdata zeros(size(imdata,1), size_dif/2 +1)]; %evenly pad zeros to square image, make divisble by 4
% imdata = [imdata;zeros(2,size(imdata,2))];
rows = size(imdata,1);
cols = size(imdata,2);


ft_weight = 1/sqrt(size(imdata,1)*size(imdata,2));
F_imdata = fftshift(fft2(imdata).*ft_weight);

%Use gaussian mask for testing
percent_sp = 0.2;

mse_iters = zeros(1,10);
psnr_iters = zeros(1,10);
for big_iter = 1:10 %run 10 sets of random masks and collect average data
    mask = make_gauss_mask(rows,percent_sp);
    %figure; imshow(abs(mask));title("Mask Image");

    F_imdata_sp = F_imdata.*mask;
    im_sp = ifft2(ifftshift(F_imdata_sp))./ft_weight;
    im_sp = im_sp./(abs(max(im_sp,[],'all')));
    im_og = im_sp;


    im_final = zeros(rows,cols);

    % iter_length = 55;
    mean_squared_error = ones(1,iter_length);
    snr = zeros(1,iter_length);
    % dwtmode('per');
    for n = 1:iter_length
        %wavelet transform
        %     wname = 'db16';

        [A1,H1,V1,D1] = dwt2(im_sp,wname); %take wavelet transform w/ haar filter, h,v,d are detail matrices, a is coarsest
        [A2, H2, V2, D2] = dwt2(A1, wname);

        thresh1 = abs(0.020*max([A1 H1 V1 D1],[],'all'));

        %     A1_th = C1(1:wav_dim1,1:wav_dim1);
        H1_th = threshold2(H1, thresh1);
        V1_th = threshold2(V1, thresh1);
        D1_th = threshold2(D1, thresh1);
        A2_th = threshold2(A2, thresh1);
        H2_th = threshold2(H2, thresh1);
        V2_th = threshold2(V2, thresh1);
        D2_th = threshold2(D2, thresh1);

        A1_th = idwt2(A2_th,H2_th,V2_th,D2_th, wname);
        im_sp_th = idwt2(A1_th,H1_th,V1_th,D1_th, wname); %reconstructed image after thresholding
        %go to k space and downsample
        F_sp_th = fftshift(fft2(im_sp_th).*ft_weight);
        F_sp_th_masked = F_sp_th.*(1-mask); %find k-space info to fill in gaps


        %fill in gaps, find iteration's "difference image"
        F_err = F_sp_th_masked + F_imdata_sp;
        delta_im = ifft2(ifftshift(F_err))./ft_weight;

        %add diff im to original in image space
        new_im = im_sp + delta_im;
        im_sp = new_im./abs(max(new_im, [],'all'));
        im_final = im_final + delta_im;
        im_final = im_final./abs(max(im_final,[],'all'));
        im_final(im_final < 10^(-10)) = 0;
        mean_squared_error(n) = immse(imdata, abs(im_final));
        snr(n) = psnr(abs(im_final),imdata);
    end

    % minV = min(min(abs(imdata)));
    % maxV = max(max(abs(imdata)));
    % figure;
    % subplot(1,3,1); imshow(abs(imdata), [minV maxV]);
    % title("Orignial Image");
    % subplot(1,3,2); imshow(abs(im_og), [minV maxV]);
    % title("Sparse Image");
    % subplot(1,3,3); imshow(abs(im_final), [minV maxV]);
    % title("Final image");
    %
    %
    % figure; plot(mean_squared_error)
    % title("Mean Square Error vs Iterations")
    %
    % figure; plot(snr)
    % title("Peak SNR vs iterations")
    %
    % disp("Final MSE: "+mean_squared_error(end));
    % disp("Final PSNR: "+snr(end));
    mse_iters(big_iter) = mean_squared_error(end);
    psnr_iters(big_iter) = snr(end)-snr(1);

end
mse = mean(mse_iters);
peak_snr = mean(psnr_iters);


end
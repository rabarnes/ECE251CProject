function [imdata, im_og, im_final, mean_squared_error, peaksnr] = dddt_CS_loop(imdata,PDF, mask, iter_length, threshold_weight)
rows = size(imdata,1);
cols = size(imdata,2);

ft_weight = 1/sqrt(size(imdata,1)*size(imdata,2));
F_imdata = fftshift(fft2(imdata).*ft_weight);


F_imdata_sp = F_imdata.*mask;
im_sp = ifft2(ifftshift(F_imdata_sp))./ft_weight;
im_sp = im_sp./(abs(max(im_sp,[],'all')));
im_og = im_sp;
im_sp = im_sp./PDF;

ft_weight = 1/sqrt(size(im_sp,1)*size(im_sp,2));
rows = size(im_sp, 1);
cols = size(im_sp, 2);

im_final = zeros(rows,cols);
mean_squared_error = ones(1, iter_length);
peaksnr = ones(1, iter_length);

%%% For iterative DWT
cA_hist = struct;
cH_hist = struct;
cV_hist = struct;
cD_hist = struct;
n_dwt = 2;
dwt_fname = 'bior1.1';

for n = 1:iter_length
    % choose type of wavelet transform
%     im_sp_th = dwtthresh(im_sp, 1, 'haar', threshold_weight);
%     im_sp_th = dualtreethresh(im_sp, 2, threshold_weight);
    im_sp_th = dddtreethresh(im_sp, 2, 'dddtf1', threshold_weight);

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


end
function [imdata, im_og, im_final, mean_squared_error, peaksnr] = basic_CS_loop(imdata,PDF, mask, iter_length, threshold_weight)
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

for n = 1:iter_length
    %wavelet transform
    [A1,H1,V1,D1] = dwt2(im_sp,'haar'); %take wavelet transform w/ haar filter, h,v,d are detail matrices, a is coarsest
    [A2, H2, V2, D2] = dwt2(A1, 'haar');

    % thresh_pct = 0.8; %set bottom x% to zero, use 5% value as threshold
    % C_sorted = sort(C,'ascend');
    % thresh = abs(C_sorted(round(thresh_pct*size(C,2))));
    C1 = [A1 H1 V1 D1];
    C2 = [A2 H2 V2 D2];
    wav_dim1 = size(A1,1);
    wav_dim2 = size(A2,1);
    thresh1 = abs(threshold_weight*max(C1,[],'all'));
    for i = 2:4
        wav_dim = wav_dim1;
        curr_wav = C1(1:wav_dim, wav_dim*(i-1)+1:i*wav_dim);
        C1(1:wav_dim, wav_dim*(i-1)+1:i*wav_dim) = threshold2(curr_wav,thresh1);

    end

%     A1_th = C1(1:wav_dim1,1:wav_dim1);
    H1_th = C1(1:wav_dim1,wav_dim1+1:2*wav_dim1);
    V1_th = C1(1:wav_dim1,2*wav_dim1+1:3*wav_dim1);
    D1_th = C1(1:wav_dim1,3*wav_dim1+1:4*wav_dim1);

    for i = 1:4
        wav_dim = wav_dim2;
        curr_wav = C2(1:wav_dim, wav_dim*(i-1)+1:i*wav_dim);
        C2(1:wav_dim, wav_dim*(i-1)+1:i*wav_dim) = threshold2(curr_wav,thresh1);

    end

    A2_th = C2(1:wav_dim2,1:wav_dim2);
    H2_th = C2(1:wav_dim2,wav_dim2+1:2*wav_dim2);
    V2_th = C2(1:wav_dim2,2*wav_dim2+1:3*wav_dim2);
    D2_th = C2(1:wav_dim2,3*wav_dim2+1:4*wav_dim2);

    A1_th = idwt2(A2_th,H2_th,V2_th,D2_th,'haar');
    im_sp_th = idwt2(A1_th,H1_th,V1_th,D1_th,'haar'); %reconstructed image after thresholding
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
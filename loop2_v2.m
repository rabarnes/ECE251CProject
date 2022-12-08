close all;



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
% mask = make_mask(rows, 4);
% mask = make_gauss_mask(rows, 1);
[mask, percent] = make_spiral_mask(rows, 1);
percent
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
iter_length = 100;
threshold_weight = 0.016; %for cart


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


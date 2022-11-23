close all
brain = imread('brain.jpg');
imdata = im2double(brain(:,:,1)); %black and white image, all layers same, conv to double [0 1]
size_dif = size(imdata,1)-size(imdata,2);
imdata = [zeros(size(imdata,1), size_dif/2 +1) imdata zeros(size(imdata,1), size_dif/2 +1)]; %evenly pad zeros to square image, make divisble by 4
imdata = [imdata;zeros(2,size(imdata,2))];
rows = size(imdata,1);
cols = size(imdata,2);

ft_weight = 1/(norm(imdata,"fro"));

F_imdata = fftshift(fft2(imdata)*ft_weight);

m= 4;
mask = make_mask(rows,m);
F_imdata_sp = F_imdata.*mask;
im_sp = ifft2(ifftshift(F_imdata_sp))/ft_weight;
im_sp = im_sp./(abs(max(im_sp,[],'all')));
im_og = im_sp;

im_final = zeros(rows,cols);

%Loop here
for n = 1:50
    %wavelet transform
    [A1,H1,V1,D1] = dwt2(im_sp,'haar'); %take wavelet transform w/ haar filter, h,v,d are detail matrices, a is coarsest
    [A2, H2, V2, D2] = dwt2(A1, 'haar');

    C1 = [A1 H1 V1 D1];
    C2 = [A2 H2 V2 D2];
    wav_dim1 = size(A1,1);
    wav_dim2 = size(A2,1);
    thresh1 = abs(0.01*max(C1,[],'all'));

    for i = 2:4
        wav_dim = wav_dim1;
        curr_wav = C1(1:wav_dim, wav_dim*(i-1)+1:i*wav_dim);
        if i==2
            C1(1:wav_dim, wav_dim*(i-1)+1:i*wav_dim) = threshold(curr_wav,10*thresh1);
        else
            C1(1:wav_dim, wav_dim*(i-1)+1:i*wav_dim) = threshold(curr_wav,thresh1);
        end

    end
    % A_th = C1(1:wav_dim1,1:wav_dim1);
    H1_th = C1(1:wav_dim1,wav_dim1+1:2*wav_dim1);
    V1_th = C1(1:wav_dim1,2*wav_dim1+1:3*wav_dim1);
    D1_th = C1(1:wav_dim1,3*wav_dim1+1:4*wav_dim1);

    for i = 1:4
        wav_dim = wav_dim2;
        curr_wav = C2(1:wav_dim, wav_dim*(i-1)+1:i*wav_dim);
        if i==2
            C2(1:wav_dim, wav_dim*(i-1)+1:i*wav_dim) = threshold(curr_wav,10*thresh1);
        else
            C2(1:wav_dim, wav_dim*(i-1)+1:i*wav_dim) = threshold(curr_wav,thresh1);
        end
    end

    A2_th = C2(1:wav_dim2,1:wav_dim2);
    H2_th = C2(1:wav_dim2,wav_dim2+1:2*wav_dim2);
    V2_th = C2(1:wav_dim2,2*wav_dim2+1:3*wav_dim2);
    D2_th = C2(1:wav_dim2,3*wav_dim2+1:4*wav_dim2);

    A1_th = idwt2(A2_th,H2_th,V2_th,D2_th,'haar');
    im_sp_th = idwt2(A1_th,H1_th,V1_th,D1_th,'haar'); %reconstructed image after thresholding

    %go to k space and downsample
    F_sp_th = fftshift(fft2(im_sp_th)*ft_weight);
    F_sp_th_masked = F_sp_th.*mask;


    %find err in k space and ifft to get difference image (compare to
    %original)
    F_err = F_sp_th_masked - F_imdata_sp; %Ax-y error between sampled original and wavelet recon in k space
    delta_im = ifft2(ifftshift(F_err))/ft_weight;

    %add diff im to original in image space
    new_im = im_sp + delta_im;
    im_sp = new_im./abs(max(new_im, [],'all'));
    im_final = im_final + delta_im;
    im_final = im_final./abs(max(im_final,[],'all'));

end

figure();
imshow(abs(im_sp));
title("im_sp");

figure();
imshow(abs(im_final));
title("im_final");

figure();
imshow(abs(im_og));
title("im og");
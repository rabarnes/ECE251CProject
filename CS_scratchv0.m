clear
close all
clc

%% MRI Image
im = phantom('Modified Shepp-Logan',200);
figure; imshow(abs(im));
title('Shepp-Logan Image');

%% step 1 generate measured k space
% creating sampling mask and applying to image
dim_im = size(im);
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
gaussian_mask = zeros([dim_im(1) dim_im(2)]);
for i=1:num_subsample
    gaussian_mask(R(i,1), R(i,2)) = 1;
end
% figure imshow(gaussian_mask);
% title('Gaussian Mask');

% step 1 still
F_im = fftshift(fft2(im));
F_im = ((1-0)*(F_im- min(F_im(:)))) / (max(F_im(:))-min(F_im(:)));
% figure; imshow(abs(F_im));
% title('Shepp-Logan K-space');


% create sparse image by multiplying mask & k-space
F_im = F_im .* gaussian_mask;
sparse_image = ifft2(ifftshift(F_im));
sparse_image = ((1-0)*(sparse_image - min(sparse_image(:)))) / (max(sparse_image(:))-min(sparse_image(:)));
figure; imshow(abs(sparse_image));
title("Received image");

[LoD,HiD,LoR, HiR] = wfilters('haar');

debug = 0;
n_max = 100;
n = 0;
epsilon = 10;
del = 10000;
difference_image = zeros(dim_im(1), dim_im(2));
error_graph = zeros([1 n_max]);

while (n < n_max && del > epsilon)
    
    %% step 2 && 10
    sparse_image = sparse_image + difference_image;
    sparse_image = ((1-0)*(sparse_image - min(sparse_image(:)))) / (max(sparse_image(:))-min(sparse_image(:)));
    if (debug == 1)
        figure; imshow(abs(sparse_image));
        title('reconstructed image after adding diff');
    end

    %% step 3
    [cA,cH,cV,cD] = dwt2(sparse_image,LoD,HiD,'mode','symh');
    if (debug == 1) 
        figure; imshow(abs(cA))
        title("low pass of wavelet transform of image")
    end
    
    %% step 4 thresholding
    m = sort(abs([cA(:) ; cH(:) ; cV(:) ; cD(:)]), 'descend');
    ndx = floor(length(m)*0.1);
    thresh = m(ndx);

    cA_denoise = cA .* (abs(cA) > thresh);
    cD_denoise = cD .* (abs(cD) > thresh);
    cH_denoise = cH .* (abs(cH) > thresh);
    cV_denoise = cV .* (abs(cV) > thresh);

    if (debug == 1)
        figure; imshow(abs(cA_denoise))
        title("wavelet transform of image after thresh")
    end

    %% step 5 inv W trans
    inv_c = idwt2(cA_denoise, cH_denoise, cV_denoise, cD_denoise, LoR, HiR);
%     inv_c = ((1-0)*(inv_c - min(inv_c(:)))) / (max(inv_c(:))-min(inv_c(:)));

    if (debug == 1)
        figure; imshow(abs(inv_c));
        title('Reconstructed Wavelet image after thresholding')
    end
    %% step 6 k space of denoised img
    F_im_2 = fftshift(fft2(inv_c));
    F_im_2 = ((1-0)*(F_im_2 - min(F_im_2(:)))) / (max(F_im_2(:))-min(F_im_2(:)));

    if (debug == 1)
        figure; imshow(abs(F_im_2));
        title('Shepp-Logan K-space after wavelet filtering');
    end
    
    %% step 7 apply mask from step 1
    F_im_2 = F_im_2 .* gaussian_mask;
    if (debug == 1)
        figure; imshow(abs(F_im_2));
        title('Shepp-Logan K-space after applying mask');
    end

    %% step 8 find diff kspace
    error_im = F_im_2 - F_im;

    if (debug == 1)
        figure; imshow(abs(error_im));
        title('error image');
    end

    %% step 9 find diff img
    difference_image = ifft2(ifftshift(error_im));
    if (debug == 1)
        figure; imshow(abs(difference_image));
        title('difference image');
    end
    % calc delta and determine if need to continue optimizing image
%     del = sum(abs(difference_image(:)));
    del = immse(difference_image, sparse_image);
    error_graph(n+1) = del; % save delta to plot at end
    n = n+1;
end

figure; plot(error_graph);
figure; imshow(abs(sparse_image), []); title("final image");


%% Wavelet Transform
% wt = dddtree2(typetree,x,level,fdf,df)
% wt = dddtree2(typetree,x,level,fname)
% wt = dddtree2(typetree,x,level,fname1,fname2)
% wt = dddtree2('ddt', abs(im), 1, 'filters1');
% 
% wt = dddtree2('cplxdt', abs(im), 1, 'FSfarras', 'qshift10');
% plotdt(wt);
% 
% 
% %% Check that DWT of works properly when not spare
% 
% xrec = idddtree2(wt);
% % Output the max error of any one pixel between DWT and original iamge
% max(max(abs(abs(im)-xrec))) % Should output a low number
% 
% %%
% 
% wt1 = wt;
% 
% %%
% [m,n,x,y,z] = size(wt.cfs{1});
% tempMatrix = zeros(m,n,x,y,z);
% 
% for i = 1:x
%     for j = 1:y
%         for k = 1:z
%             tempMatrix(:,:,i,j,k) = wt.cfs{1}(:,:,i,j,k);
%         end
%     end
% end
% %%
% 
% imp = sort(abs(tempMatrix(:)),'descend');
% idx = floor(length(imp)*20/100);
% threshold = imp(idx);
% 
% 
% %% Threshold
% 
% tempMatrix(abs(tempMatrix) < threshold) = 0;
% 
% for i = 1:x
%     for j = 1:y
%         for k = 1:z
%         
%             wt.cfs{1}(:,:,i,j,k) = tempMatrix(:,:,i,j,k);
% 
%         
%         end
%     end
% end
% 
% %% Show Thresholded Image and calculate max pixel error
% 
% xrec = idddtree2(wt);
% max(max(abs(abs(im)-xrec)))
% 
% 
% figure; imshow(xrec,[]);
% title('Recon Shepp-Logan Image');
% 
% 
% 
% 
% 
% 
% 
% %% Trash
% 
% 
% for i = 1:size(wt.cfs{1}, 3)
%     nexttile;
%     imagesc(wt.cfs{1}(:,:,i)); axis square;
% end
% 
% 
% % H1H1 = wt.cfs{1}(:,:,4);
% % H1H2 = wt.cfs{1}(:,:,5);
% % H2H1 = wt.cfs{1}(:,:,7);
% % H2H2 = wt.cfs{1}(:,:,8);
% 
% tiledlayout(2,4);
% 
% for i = 1:size(wt.cfs{1}, 3)
%     nexttile;
%     imagesc(wt.cfs{1}(:,:,i)); axis square;
% end
% colormap('gray');

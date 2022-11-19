clear
close all
clc

%% MRI Image
im = phantom('Modified Shepp-Logan',200);

figure; imshow(abs(im));
title('Shepp-Logan Image');


%% To k-space

ft_weight = 1/sqrt(size(im,1)*size(im,2));
F_im = fftshift(fft2(ifftshift(im)))*ft_weight;

figure; imshow(abs(F_im));
title('Shepp-Logan K-space');

%% create sampling masks
dim_im = size(im);
num_pixels = dim_im(1) * dim_im(2);
num_subsample = round(num_pixels*0.33);
%% random gaussian
u = [0 0];
sigma = [1 0; 0 1];
%% generate subsamples of multivariate gaussian
R = mvnrnd(u, sigma, round(num_subsample));
%% scale all values of R to indices of image
a = 1;
b = min(dim_im(1), dim_im(2)); % TODO fix this for rectangular images
min_R = min(R(:));
max_R = max(R(:));
R = ((b-a)*(R - min_R)) / (max_R-min_R) + a;
R = round(R);
%% create mask by setting mask indices to 1
gaussian_mask = zeros([dim_im(1) dim_im(2)]);
for i=1:num_subsample
    gaussian_mask(R(i,1), R(i,2)) = 1;
end

%% Back to image space
%% create sparse image by multiplying mask & k-space
F_im = F_im .* gaussian_mask;
sparse_k_space = F_im;
ift_weight = sqrt(size(im,1)*size(im,2));
sparse_image = fftshift(ifft2(ifftshift(F_im)))*ift_weight;

figure; imshow(abs(sparse_image));
title('Sparse Shepp-Logan');

[LoD,HiD] = wfilters('haar','d');
[cA,cH,cV,cD] = dwt2(sparse_image,LoD,HiD,'mode','symh');
subplot(2,2,1)
imagesc(abs(cA))
colormap gray
title('Approximation')
subplot(2,2,2)
imagesc(abs(cH))
colormap gray
title('Horizontal')
subplot(2,2,3)
imagesc(abs(cV))
colormap gray
title('Vertical')
subplot(2,2,4)
imagesc(abs(cD))
colormap gray
title('Diagonal')

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

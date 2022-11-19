%% MRI Image
im = phantom('Modified Shepp-Logan', 256);

im_even = padarray(im,[1 1], 0 ,'post'); 

figure; imshow(abs(im));
title('Shepp-Logan Image');


%% To k-space

ft_weight = 1/sqrt(size(im,1)*size(im,2));
F_im = fftshift(fft2(ifftshift(im)))*ft_weight;

figure; imshow(abs(F_im));
title('Shepp-Logan K-space');

%% https://people.eecs.berkeley.edu/~mlustig/Software.html <- Sparse MRI

% 
% TVWeight = 0.01; 	% Weight for TV penalty
% xfmWeight = 0.00;	% Weight for Transform L1 penalty
% Itnlim = 8;		% Number of iterations
% 
%  [pdf,val] = genPDF([256,256], 5, 0.33,1, 0, 1);

%%
L = 127;
rm = 10;


rm = rm+1;
div = 0;
for i = 1:rm
    div = div + (1/i);
end
a = ceil(L/div);

rm_lines = zeros(1, rm);
rm_cnt = linspace(rm,1,rm);


for i = 1:rm
%     rm_lines(i) = round(a*(1/rm_cnt(i)));
    if i == 1
        rm_lines(i) = round(a*(1/rm_cnt(i)));
    else
        rm_lines(i) = rm_lines(i-1) + round(a*(1/rm_cnt(i)));
    end
end
rm_lines = rm_lines(1:end-1);



kspace_choice = cat(2, 1:127, [ 129], flip(1:127));
kspace_choice = kspace_choice ./ 128;
kspace_choice(rm_lines) = 0;
kspace_choice(255-rm_lines) = 0;


figure;
plot(kspace_choice); axis tight;
title("Data Line Importance")


%% Create Sparse Image

sparse_cart_mask = ones(size(im,1),size(im,2));
sparse_cart_mask(rm_lines,:) = 0;
sparse_cart_mask(255-rm_lines,:) = 0;


figure; imshow(abs(sparse_cart_mask));
title('Sparce Image Mask');

%%
% Center row is (128,:)
sparse_k_space = F_im .* sparse_cart_mask;

figure; imshow(abs(sparse_k_space));
title('Shepp-Logan K-space');


%% Back to image space
ift_weight = sqrt(size(im,1)*size(im,2));
sparse_image = fftshift(ifft2(ifftshift(sparse_k_space)))*ift_weight;

sparse_image_even = padarray(sparse_image,[1 1], 0 ,'post'); 

figure; imshow(abs(sparse_image));
title('Sparse Shepp-Logan');




%%

y = sparse_k_space;

ift_weight = sqrt(size(y,1)*size(y,2));
x = fftshift(ifft2(ifftshift(y)))*ift_weight;

figure; imshow(abs(x));
title('Image with noise like artifacts');

for i = 1:1

        %[Wx_c, Wx_s] = wavedec2(x, 1, "db1");
        [cA,cH,cV,cD] = dwt2(x, "sym2");
        
        
        %figure; plot(abs(Wx_c));
        
%         figure;
%         subplot(2,2,1); imshow(abs(cA));
%         subplot(2,2,2); imshow(abs(cH));
%         subplot(2,2,3); imshow(abs(cV));
%         subplot(2,2,4); imshow(abs(cD));
        
        %Wx_c_thresh = Wx_c;
        %Wx_c_thresh(Wx_c_thresh <= 0.20) = 0;
        
        %dn_x = waverec2(Wx_c_thresh, Wx_s, "db1");

        cA_thresh = cA;
        %cA_thresh(cA_thresh < 0.02) = 0;
        cH_thresh = cH;
        cH_thresh(cH_thresh < 0.02) = 0;
        cV_thresh = cV;
        cV_thresh(cV_thresh < 0.02) = 0;
        cD_thresh = cD;
        cD_thresh(cD_thresh < 0.02) = 0;
%         figure; hold on;
%         plot(abs(cA)); plot(abs(cA_thresh));
%         hold off;

        dn_x = idwt2(cA_thresh,cH_thresh,cV_thresh,cD_thresh,"sym2");
        
        %figure; imshow(abs(dn_x));

        ft_weight = 1/sqrt(size(dn_x,1)*size(dn_x,2));
        dn_x_kspace = fftshift(fft2(ifftshift(dn_x)))*ft_weight;
        
        %figure; imshow(abs(dn_x_kspace));
        %title('K-Space of denoised image');
        
        Ax_kspace = dn_x_kspace .* sparse_cart_mask;
        %figure; imshow(abs(Ax_kspace));
        %title('K-Space of denoised image');

        diff_kspace = Ax_kspace - y;
        %figure; imshow(abs(diff_kspace));
        %title('Difference K-Space Ax-y');

        ift_weight = sqrt(size(diff_kspace,1)*size(diff_kspace,2));
        diff_image = fftshift(ifft2(ifftshift(diff_kspace)))*ift_weight;
        
        %figure; imshow(abs(diff_image));
        %title('Difference Image');

        
        x = x + diff_image;
        figure; imshow(abs(x));
        title('New X');

end

%% Wavelet Transform
% wt = dddtree2(typetree,x,level,fdf,df)
% wt = dddtree2(typetree,x,level,fname)
% wt = dddtree2(typetree,x,level,fname1,fname2)
% wt = dddtree2('ddt', abs(im), 1, 'filters1');

wt = dddtree2('cplxdt', abs(sparse_image_even), 1, 'FSfarras', 'qshift10');
plotdt(wt);


%% Check that DWT of works properly when not spare

xrec = idddtree2(wt);
% Output the max error of any one pixel between DWT and original iamge
max(max(abs(abs(im_even)-xrec))) % Should output a low number

%%

wt1 = wt;

%%
[m,n,x,y,z] = size(wt.cfs{1});
tempMatrix = zeros(m,n,x,y,z);

for i = 1:x
    for j = 1:y
        for k = 1:z
            tempMatrix(:,:,i,j,k) = wt.cfs{1}(:,:,i,j,k);
        end
    end
end


%%

imp = sort(abs(tempMatrix(:)),'descend');
idx = floor(length(imp)*20/100);
threshold = imp(idx);


%% Threshold

tempMatrix(abs(tempMatrix) < threshold) = 0;

for i = 1:x
    for j = 1:y
        for k = 1:z
        
            wt.cfs{1}(:,:,i,j,k) = tempMatrix(:,:,i,j,k);

        
        end
    end
end

%% Show Thresholded Image and calculate max pixel error

xrec = idddtree2(wt);
max(max(abs(abs(im_even)-xrec)))


figure; imshow(xrec,[]);
title('Recon Shepp-Logan Image');







%% Trash


for i = 1:size(wt.cfs{1}, 3)
    nexttile;
    imagesc(wt.cfs{1}(:,:,i)); axis square;
end


% H1H1 = wt.cfs{1}(:,:,4);
% H1H2 = wt.cfs{1}(:,:,5);
% H2H1 = wt.cfs{1}(:,:,7);
% H2H2 = wt.cfs{1}(:,:,8);

tiledlayout(2,4);

for i = 1:size(wt.cfs{1}, 3)
    nexttile;
    imagesc(wt.cfs{1}(:,:,i)); axis square;
end
colormap('gray');

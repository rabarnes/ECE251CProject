%% MRI Image
im = phantom('Modified Shepp-Logan', 255);

im_even = padarray(im,[1 1], 0 ,'post'); 

figure; imshow(abs(im));
title('Shepp-Logan Image');


%% To k-space

ft_weight = 1/sqrt(size(im,1)*size(im,2));
F_im = fftshift(fft2(ifftshift(im)))*ft_weight;

figure; imshow(abs(F_im));
title('Shepp-Logan K-space');


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



kspace_choice = cat(2, 1:127, [128], flip(1:127));
kspace_choice = kspace_choice ./ 128;
kspace_choice(rm_lines) = 0;
kspace_choice(255-rm_lines) = 0;


figure;
plot(kspace_choice); axis tight;
title("Data Line Importance")


%% Create Sparce Image
% Center row is (128,:)
sparse_k_space = F_im;
sparse_k_space(rm_lines,:) = 0;
sparse_k_space(255-rm_lines,:) = 0;

figure; imshow(abs(sparse_k_space));
title('Shepp-Logan K-space');


%% Back to image space
ift_weight = sqrt(size(im,1)*size(im,2));
sparse_image = fftshift(ifft2(ifftshift(sparse_k_space)))*ift_weight;

sparse_image_even = padarray(sparse_image,[1 1], 0 ,'post'); 

figure; imshow(abs(sparse_image));
title('Sparse Shepp-Logan');

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

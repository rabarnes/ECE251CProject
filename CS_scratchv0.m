%% MRI Image
im = phantom('Modified Shepp-Logan',200);

figure; imshow(abs(im));
title('Shepp-Logan Image');


%% To k-space

ft_weight = 1/sqrt(size(im,1)*size(im,2));
F_im = fftshift(fft2(ifftshift(im)))*ft_weight;

figure; imshow(abs(F_im));
title('Shepp-Logan K-space');

%% Back to image space

sparse_k_space = F_im;
ift_weight = sqrt(size(im,1)*size(im,2));
sparse_image = fftshift(ifft2(ifftshift(F_im)))*ift_weight;

figure; imshow(abs(sparse_image));
title('Sparse Shepp-Logan');

%% Wavelet Transform
% wt = dddtree2(typetree,x,level,fdf,df)
% wt = dddtree2(typetree,x,level,fname)
% wt = dddtree2(typetree,x,level,fname1,fname2)
% wt = dddtree2('ddt', abs(im), 1, 'filters1');

wt = dddtree2('cplxdt', abs(im), 1, 'FSfarras', 'qshift10');
plotdt(wt);


%%

xrec = idddtree2(wt);
max(max(abs(abs(im)-xrec)))

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


%%

tempMatrix(abs(tempMatrix) < threshold) = 0;

for i = 1:x
    for j = 1:y
        for k = 1:z
        
            wt.cfs{1}(:,:,i,j,k) = tempMatrix(:,:,i,j,k);

        
        end
    end
end

%%

xrec = idddtree2(wt);
max(max(abs(abs(im)-xrec)))
%%
figure; imshow(xrec,[]);
title('Recon Shepp-Logan Image');

%%
for i = 1:size(wt.cfs{1}, 3)
    nexttile;
    imagesc(wt.cfs{1}(:,:,i)); axis square;
end


%%

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

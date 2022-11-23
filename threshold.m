function th = threshold(A, thresh)
%Threshold matrix A by zeroing values less than thresh, and reducing all
%other values by thresh
th = zeros(size(A,1), size(A,2));
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            if(abs(A(i,j))<thresh)
                th(i,j) = 0;
            else
                th(i,j) = A(i,j) - thresh*exp(1j*angle(A(i,j))); %soft thresholding, subtract threshold val from all other pixels
            end
        end
    end


end
function th = threshold2(A, thresh)
    th = (abs(A) > thresh).*(A.*(abs(A)-thresh)./(abs(A)+eps));
end
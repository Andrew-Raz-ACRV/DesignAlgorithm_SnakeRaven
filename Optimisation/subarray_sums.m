function M = subarray_sums(A,sub_dim)
    %Splits matrix A into subarrays of size sub_dim and 
    %computes the sum of those subarrays in a new matrix
    %
    %e.g. A = ones(6,6,6), sub_dim = [2 2 1]
    % Find the sums amoung 2x2 sub arrays
    % 
    % M = subarray_sums(A,sub_dim)
    %
    % M becomes a matrix M = 4.*ones(6/2,6/2,6/1)
    % as each element is 4 the sum of every 2x2 ones array
    %
    %Written by Andrew Razjigaev 2019
    
    %subarray 3D matrix size
    a = sub_dim(1);
    b = sub_dim(2);
    c = sub_dim(3);
    
    %Matrix size
    x = size(A,1);
    y = size(A,2);
    z = size(A,3);
    
    %output size:
    n = x/a;
    m = y/b;
    o = z/c;
    
    M = zeros(n,m,o);

    %compute subarray sums
    for ii = 1:n
        for jj = 1:m
            for kk = 1:o
                %Extract subarray               
                rows = (ii-1)*a+1 : ii*a;
                cols = (jj-1)*b+1 : jj*b;
                dept = (kk-1)*c+1 : kk*c;
                subarray = A(rows,cols,dept);
                %sum of subarray
                M(ii,jj,kk) = sum(subarray(:));
            end
        end
    end
end
                
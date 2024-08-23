function cdf=measure2d_to_cdf(pdf)
 
n = size(pdf,1) * size(pdf,2);

cdf = reshape( pdf, [n,1] );

for i=2:n
    cdf(i) = cdf(i-1) + cdf(i);
end

cdf = cdf / cdf(n);
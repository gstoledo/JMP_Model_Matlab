function cdf=measure_to_cdf(pdf) 
 
n = length(pdf);

cdf=zeros(size(pdf));
cdf(1) = pdf(1);

for i=2:n
    cdf(i) = cdf(i-1) + pdf(i);
end

cdf =  cdf / cdf(n);

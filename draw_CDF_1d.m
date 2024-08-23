function j=draw_CDF_1d( cdf, u  ) %This will work for Unemp CDFs only in our case

n=length(cdf);  % length of the array CDF
j=sum(u>cdf);
j = min(j+1, n);
 


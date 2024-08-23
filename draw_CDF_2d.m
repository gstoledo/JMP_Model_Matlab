%Original
%function [j,k]=draw_CDF_2d( cdf, nh, u)
 
      % n = size(cdf,1);  
      % i=sum(u>cdf);  
      % i = min(i+1, n);    
      
      % tmp = (i-1) / nh;  

      % j = floor( tmp ) + 1;
      % k = i - nh * ( j-1 );

      % For our case
function [a,z]=draw_CDF_2d( cdf, ats,tpts, u)
 
      n = size(cdf,1);  
      i=sum(u>cdf);  
      i = min(i+1, n);
      [a,z]=ind2sub([ats,tpts],i); %Cant belive this function exists
          
      
      
 
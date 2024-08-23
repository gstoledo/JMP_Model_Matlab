%For our case
function [a,z,q]= draw_CDF_3d( cdf, ats,tpts, u)
  
      n = size(cdf,1);  
      i=sum(u>cdf);  
      i = min(i+1, n);
      [a,z,q]=ind2sub([ats,tpts,tpts],i); %Cant belive this function exists
          
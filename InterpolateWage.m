function wage=InterpolateWage( W_w, W, wage_w )
%Interpolate wage-- given new match values to worker W_w, new value to
%deliver W, and wage grid wage_w, returns the wage
 
n = length(W_w);

 
j=sum(W>W_w);
if j > 0  &&  j < n
    wage = wage_w(j)  + ( (W - W_w(j)) / (W_w(j+1) - W_w(j)) ) * ( wage_w(j+1) - wage_w(j) );
elseif  j == 0
    wage = wage_w(1);
else 
    wage = wage_w(n);
end

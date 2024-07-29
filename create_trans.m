function output=create_trans(down,stay,up,N)
%create_trans: output=create_trans(down,stay,up,N)
%Transition matrix with elements (i,j)=odds of going from state i-> state j
%down= the probability of transition down one state
%stay= prob of staying
%up = prob moving up one state
%N=number of states

trans=zeros(N,N); %Transition matrix
trans(1,:)=[down+stay,up,zeros(1,N-2)];
for ii=2:N-1
    trans(ii,:)=[zeros(1,ii-2),down,stay,up,zeros(1,N-(ii+1))];
end
trans(end,:)=[zeros(1,N-2),down,stay+up];


output=trans;

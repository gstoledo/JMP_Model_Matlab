% First attempt at computing the value functions for my model 


%Lets for now import some valeus from HLMP code
clear all
close all
%load('xmin_results_combined')
load('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit/xmin_results_combined.mat')
x=xmin;



lamu                    =x(1);
lam                    =x(2);
fcomp                   =x(3); %Complementarity in F
mubar                   =x(4);
del                     =x(5);
ulose                   =x(6);
sologain                =x(7);
team_learn_param_up  	=x(8);
mnew	                =x(9);
team_learn_param_up_sym	=x(10);
bpf                     =x(11);
homeprod                =x(12);
sololose                =0;
A                       =x(14);
team_learn_param_down   =0;
mnew_high               =floor(x(15)); 



%Fundamentals
bt   =1/(1.1^(1/12))   ; %beta, discount factor
death=1/(35*12)        ; %Probability agent dies  
bpw  =1-bpf            ; %bpw is bargaining power of worker
nfirm=1                ; %Measure of firms 
tpts =3                ; %type point space
ats=2                  ; %Productivity space 
spts =tpts+2           ; %state S for updating -- {u,0,{j}} -- dimension of that space is type+2

typemin=1; %Lowest type
typemax=mubar ; %Highest type

amin=1; %Lowest productivity
amax=2; %Highest productivity

type=[typemin:(typemax-typemin)/(tpts-1):typemax]; %Types

% typebirth=zeros(1,tpts);
% for i=1:mnew_high
%     typebirth(i)= mnew*exp(-type(i)*mnew)./sum(mnew*exp(-type(1:mnew_high)*mnew)); %Birth rate
% end


%Value Functions 
zero_tpts=zeros(1,tpts);
zero_tpts_tpts=zeros(tpts,tpts);
zero_ats=zeros(1,ats);
zero_tpts_ats=zeros(1,tpts,ats);
zero_tpts_tpts_ats=zeros(tpts,tpts,ats);

% Uh=zeros(1,tpts); %Unemployed value function
% Veh=zeros(1,ats);  %Empty firm value function
% Vmh=zeros(1,tpts,ats); %Firm with manager value function
% Vnh=zeros(1,tpts,ats); %Firm with no manager value function
% Vth=zeros(tpts,tpts,ats); %Firm with team value function



%Some number to test it out 
Veh=[0 1] ;
Vmh(:,:,1)=[1 2 3];
Vmh(:,:,2)=[4 5 6];
% Vmh(:,:,3)=[7 8 9];

Vnh(:,:,1)=[0 1 1];
Vnh(:,:,2)=[1 2 2];
% Vnh(:,:,3)=[2 3 3];

Vth(:,:,1)=[0 1 1; 1 8 2];
Vth(:,:,2)=[1 2 2; 2 3 3];
% Vth(:,:,3)=[2 3 3; 3 4 4; 4 5 5];
Uh=[0 1 1];

%Manager penalty toggle 
nm_penal=1;
%%%%Empty firm continuation value 

%Gains from trade from meeting unemployed
gt_meet_u= max(max(Vmh,Vnh) - repmat(reshape(Veh,1,1,ats),1,tpts) -repmat(Uh,1,1,ats), zero_tpts_ats)

%Gains from trade from meeting firm with manager 
%Pages are my own a
%a_tilda are rows
% Z tilda are columns
gt_meet_m= max(repmat(max(Vmh,Vnh),ats,1) - repmat(reshape(Veh,1,1,ats),ats,tpts) - out_opt_lose1(Vmh,Veh,tpts,ats), repmat(zero_tpts_ats,ats,1))

%Gains from trade from meeting firm with non-manager
%Pages are my own a
%a_tilda are rows
% Z tilda are columns
gt_meet_nm= max(repmat(max(Non_man_penalty(Vmh,tpts,ats,nm_penal),Non_man_penalty(Vnh,tpts,ats,nm_penal)),ats,1) - repmat(reshape(Veh,1,1,ats),ats,tpts)...
    - out_opt_lose1(Vnh,Veh,tpts,ats), repmat(zero_tpts_ats,ats,1))

%Gains from trade from meeting firm with team
%Pages are my own a
%Subpages are a_tilda
% Z tilda are rows
% Q tilda are columns
%Meet a team and get the manager 
gt_meet_t_m= max(permute(repmat(permute(max(Vmh,Vnh),[2 1 3]),1,tpts,1,ats),[1 2 4 3])...
    - permute(repmat(reshape(Veh,1,1,ats),tpts,tpts,1,ats), [1 2 4 3]) - repmat(Vth,1,1,1,ats)+repmat(Vnh,tpts,1,1,ats), repmat(zero_tpts_tpts_ats, 1,1,1,ats));      
%Meet a team and get the non-manager
alloc_nm=max(Non_man_penalty(Vmh,tpts,ats,nm_penal),Non_man_penalty(Vnh,tpts,ats,nm_penal));
gt_meet_t_nm= max(permute(repmat(alloc_nm,tpts,1,1,ats), [1 2 4 3])...
    - permute(repmat(reshape(Veh,1,1,ats),tpts,tpts,1,ats), [1 2 4 3]) - repmat(Vth,1,1,1,ats) + repmat(permute(Vmh, [2 1 3]),1,tpts,1,ats) , repmat(zero_tpts_tpts_ats, 1,1,1,ats));

%Final Gains from trade from meeting team
gt_meet_t= max(gt_meet_t_m,gt_meet_t_nm)


 

















clear all
close all
location="local";
% location="hpc";
if location == "local"
    cd('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
    addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
    addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit')
end
load('xmin_results_combined.mat')

%Baseline parameters and objects
   
%load('xmin_results_combined')
load('xmin_results_combined.mat')
x=xmin;


lamu                    =x(1);
% lamu=0.1;
lam                    =x(2);
% lam=0.1;
fcomp                   =x(3); %Complementarity in F
mubar                   =x(4);
del                     =x(5); %Delta 
ulose                   =x(6);
sologain                =x(7);
team_learn_param_up  	=x(8);
mnew	                =x(9);
team_learn_param_up_sym	=x(10);
bpf                     =x(11); %Firm bargaining power (1-gamma)
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
if location=="hpc"
    tpts=7;
    ats=5;
end
spts =tpts+2           ; %state S for updating -- {u,0,{j}} -- dimension of that space is type+2
cost_p=1                   ; %cost of promoting a non manager to manager
cost_d=1                   ; %cost of demoting a manager to non manager
alpha_m=1              ; %Manager returns
alpha_n=0.1              ; %Non manager returns



typemin=1; %Lowest type
typemax=mubar ; %Highest type

amin=1; %Lowest productivity
amax=10; %Highest productivity

type=[typemin:(typemax-typemin)/(tpts-1):typemax]; %Types
a_type=[amin:(amax-amin)/(ats-1):amax]; %Productivity


%Need to ajdust this later, pay attention to normalization to measure of workers 
mnew_high=tpts; %Number of types
typebirth=zeros(1,tpts);
for i=1:mnew_high
    typebirth(i)= mnew*exp(-type(i)*mnew)./sum(mnew*exp(-type(1:mnew_high)*mnew)); %Birth rate
end


soloprod=1    ;
pnew    =fcomp;
pold    =1    ;
 
fteam=zeros(ats,tpts,tpts);
fman=zeros(ats,tpts);
fnman=zeros(ats,tpts);
fe=zeros(1,ats);
for a=1:ats
    for z=1:tpts
        for q=1:tpts
            fteam(a,z,q)=A*(type(z).^alpha_m)*(type(q).^alpha_n);
        end
        fman(a,z)=A*a_type(a)*type(z).^alpha_m;
        fnman(a,z)=0;
    end
    fe(a)=0;
end

b=homeprod*fman(1,:) ; %Type home production vector
n=1; %mass of firms 


%Transtion matrices (They are not perfect yet just to write it all down)
%A transition
aup=0.01;
adown=0.3;
astay=1-aup-adown;
a_trans=create_trans(adown,astay,aup,ats);

%Q transition (Old one, does not depend on a)
% qup=0.02;
% qdown=0.01;
% qstay=1-qup-qdown;
% q_trans=create_trans(qdown,qstay,qup,tpts);

%Q transition (New one, depends on a), change the qup and down to depend, I am just introducing the notation
q_trans=zeros(tpts,tpts,ats);
qup=0.02*ones(ats,1);
qdown=0.01*ones(ats,1);
qstay=1-qup-qdown;

for a=1:ats
    q_trans(:,:,a)=create_trans(qdown(a),qstay(a),qup(a),tpts);
end



%Unemp transition
ugain=0.00           ; %Probability unemployed move up
ustay=1-ulose-ugain  ; %Probability unemployed stay
u_trans=create_trans(ulose,ustay,ugain,tpts);


%Toggles
dayy          ='1'   ; %Version of the code    
zero_tol      =1e-10 ; %Tolerance for zero
use_guess     ='y'   ;
tr=0               ; %Manager penalty toggle
speed=1              ; %Convergenge speed
% %Specification 0


% run run_sp0.m

% %Specification 1
% run run_sp1.m

% %Specification 2
% run run_sp2.m
% tic;
%Specification 3
run run_sp3.m
% toc;
% %Specification 4
% run run_sp4.m

% % Specification 5 
% run run_sp5.m 

% %Specification 6
% run run_sp6.m



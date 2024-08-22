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
tpts =7                ; %type point space
ats=5                  ; %Productivity space 
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
            fteam(a,z,q)=A*a_type(a)*(type(z).^alpha_m)*(type(q).^alpha_n);
        end
        fman(a,z)=A*a_type(a)*type(z).^alpha_m;
        fnman(a,z)=0;
    end
    fe(a)=0;
end

b=homeprod*fman(1,:) ; %Type home production vector
n=1; %mass of firms 


%Transtion matrices (They are not perfect yet just to write itw all down)
%A transition
aup=0.01;
adown=0.3;
astay=1-aup-adown;
a_trans=create_trans(adown,astay,aup,ats);

%Q transition
qup=0.02;
qdown=0.01;
qstay=1-qup-qdown;
q_trans=create_trans(qdown,qstay,qup,tpts);


%Unemp transition
ugain=0.00           ; %Probability unemployed move up
ustay=1-ulose-ugain  ; %Probability unemployed stay
u_trans=create_trans(ulose,ustay,ugain,tpts);




%Toggles
dayy          ='1'   ; %Version of the code    
zero_tol      =1e-10 ; %Tolerance for zero
use_guess     ='y'   ;
true=0               ; %Manager penalty toggle
speed=1              ; %Convergenge speed

%Join the loop
%Initial guesses for the LOM
%Alternative, start with all empty firms and mass 1 of unemployed
eplus_udist=(1/tpts)*ones(1,tpts);
eplus_edist=(n/ats)*ones(1,ats);
eplus_mdist=zeros(ats,tpts); %Distribution of firms with manager e_m(a,z)
eplus_ndist=zeros(ats,tpts); %Distribution of firms with no manager e_n(a,q)
eplus_tdist=zeros(ats,tpts,tpts); %Distribution of firms with team e_t(a,z,q)


e_udist=eplus_udist;
e_edist=eplus_edist;
e_mdist=eplus_mdist;
e_ndist=eplus_ndist;
e_tdist=eplus_tdist;




n0=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
pop0=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");

%VFs to initialize
%Inital guesses for value functions
Veini=zeros(1,ats);
Vmini=fman/(1-bt);
Vnini=fnman/(1-bt);
Vtini=fteam/(1-bt);
Uini=b/(1-bt);


%% VF Iteration
[Ve, Vm, Vn, Vt, U, Veh, Vmh, Vnh, Vth, Vetl, Vmtl, Vntl, Vttl, Utl]= vf_iterationV2(eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist,eplus_udist,Veini,Vmini,Vnini,Vtini,Uini,ats,tpts,cost_d,cost_p,true,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);
Uh=U;

%%  Wages
wmin  =  -6*min(fteam(:)); 
wmax = max(fteam(:));
wpts = 100;
 
wgrid = [wmin:(wmax - wmin)/(wpts - 1):wmax];

%% Wage iteration

%Inital guesses for wages
if exist('wages_a'+string(ats)+'_z'+string(tpts)+'.mat')==2 && use_guess(1)=='y'
    load('wages_a'+string(ats)+'_z'+string(tpts)+'.mat');
    Wmini=Wm;
    Wnini=Wn;
    Wtmini=Wtm;
    Wtnini=Wtn;
else
    Wmini=ones(wpts,ats,tpts)*wmin;
    Wnini=ones(wpts,ats,tpts)*wmin;
    Wtmini=ones(wpts,ats,tpts,tpts)*wmin;
    Wtnini=ones(wpts,ats,tpts,tpts)*wmin;
end


tic
[Wm,Wn,Wtm,Wtn]=wf_iteration(wpts,ats,tpts,Ve,Vm,Vn,Vt,U,Vmh,Vnh,Vth,Veh,Wmini,Wnini,Wtmini,Wtnini,...
speed,cost_d,cost_p,wgrid,bt,death,del,lam,lamu,bpw,n,a_trans,q_trans,e_udist,e_edist,e_mdist,e_ndist,e_tdist);
toc
save('wages_a'+string(ats)+'_z'+string(tpts)+'.mat','Wm','Wn','Wtm','Wtn')




            


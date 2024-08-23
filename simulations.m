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
            fteam(a,z,q)=A*a_type(a)*(type(z).^alpha_m)*(type(q).^alpha_n);
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

% %Q transition
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
true=0               ; %Manager penalty toggle
speed=1              ; %Convergenge speed


%Load some results that we will simullate over
load('model_solutions_sp3'+location)
load('wages_a'+string(ats)+'_z'+string(tpts))

e_udist=eplus_udist;
e_edist=eplus_edist;
e_mdist=eplus_mdist;
e_ndist=eplus_ndist;
e_tdist=eplus_tdist;

% Distributions 
%Firm type (unweighted)
firm_dist=[sum(e_edist),sum(e_mdist,"all"),sum(e_ndist,"all"),sum(e_tdist,"all")]/n;
firm_dist_pdf=firm_dist/sum(firm_dist);

%CDF os mass of each type of firm 
cdf_firm_dist=zeros(1,4);
for i=1:4
    cdf_firm_dist(i)=sum(firm_dist(1:i));
end
cdf_firm_dist(end)=1; %Just to make sure it is 1

%Employee type wighted for the arival rates. This will be useful to get the probability if a match
emp_dist=[lamu*sum(e_udist),lam*sum(e_mdist,"all"),lam*sum(e_ndist,"all"),lam*sum(e_tdist,"all")];
emp_dist_pdf=emp_dist/sum(emp_dist);

%Prob of matching 
prob_matching=sum(emp_dist);

% Now we do a cdf, considering the effective mass of each type and relative to the probability of matching
cdf_emp_dist=zeros(1,4);
for i=1:4
    cdf_emp_dist(i)=sum(emp_dist(1:i))/prob_matching;
end
cdf_emp_dist(end)=1; %Just to make sure it is 1

% Tranform measures into cdfs (multi-dimensional)
cdf_e_udist=measure_to_cdf(e_udist);
cdf_e_edist=measure_to_cdf(e_edist);
cdf_e_mdist=measure2d_to_cdf(e_mdist); %Will have to be careful here when making the draw. Do I need to map back to the original coordinates?
cdf_e_ndist=measure2d_to_cdf(e_ndist);
cdf_e_tdist=measure3d_to_cdf(e_tdist);


%% Adjust probabilities of shocks transitions



%% Small simulation for firms
n_months=10;
n_firms=15;
n_workers=n_firms;

%Random matrices
rng(1,"twister");

r.state=rand(n_months,n_firms); %Inside the meetings, which type of meeting is happening
r.event=rand(n_months,n_firms); % What kind of meeting is happening
r.ashock=rand(n_months,n_firms); %Productivity shock
r.qshock=rand(n_months,n_firms); %Learning shock


% Storage of firm status
firm_status=zeros(n_months,n_firms); %(1=empty, 2= Firm with manager, 3= Firm with non manager; 4= Firm with team)
m_ftp=zeros(n_months,n_firms); %Manager type
n_ftp=zeros(n_months,n_firms); %Non manager type
a_ftp=zeros(n_months,n_firms); %Productivity type
m_fwage=zeros(n_months,n_firms); %Manager wage
n_fwage=zeros(n_months,n_firms); %Non manager wage


%Inital conditions
t=1;

for i=1:n_firms  
    firm_status(t,i)=draw_CDF_1d(cdf_firm_dist,r.state(t,i));
    if (firm_status(t,i)==1) %Empty firm
        a_ftp(t,i)=draw_CDF_1d(cdf_e_edist,r.ashock(t,i));
        m_ftp(t,i)=0;
        n_ftp(t,i)=0;
        m_fwage(t,i)=0.0;
        n_fwage(t,i)=0.0;
    elseif (firm_status(t,i)==2) %Firm with manager
        [a_ftp(t,i),m_ftp(t,i)]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.event(t,i));
        n_ftp(t,i)=0;
        target_wage=(1-bpw)*U(m_ftp(t,i))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
        m_fwage(t,i)= InterpolateWage(Wm(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
        n_fwage(t,i)=0.0;
    elseif (firm_status(t,i)==3) %Firm with non manager
        [a_ftp(t,i),n_ftp(t,i)]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.event(t,i));
        m_ftp(t,i)=0;
        m_fwage(t,i)=0.0;
        target_wage=(1-bpw)*U(n_ftp(t,i))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
        n_fwage(t,i)= InterpolateWage(Wn(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
    else
        %Firm with team
        [a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.event(t,i));
        target_wage_m=(1-bpw)*U(m_ftp(t,i))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
        target_wage_n=(1-bpw)*U(n_ftp(t,i))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
        m_fwage(t,i)= InterpolateWage(Wtm(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage_m,wgrid);
        n_fwage(t,i)= InterpolateWage(Wtn(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage_n,wgrid);
    end
end    
  


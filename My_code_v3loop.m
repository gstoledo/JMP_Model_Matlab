
% First attempt at computing the value functions for my model 


%Lets for now import some valeus from HLMP code
clear all
close all
cd('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit')
%load('xmin_results_combined')
load('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit/xmin_results_combined.mat')
x=xmin;



lamu                    =x(1);
% lamu=0.1;
lam                    =x(2);
lam=lamu;
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

%% Simplifying to get the mistake

%Fundamentals
bt   =1/(1.1^(1/12))   ; %beta, discount factor
death=1/(35*12)        ; %Probability agent dies  
bpw  =1-bpf            ; %bpw is bargaining power of worker
nfirm=1                ; %Measure of firms 
tpts =5                ; %type point space
ats=2                ; %Productivity space 
spts =tpts+2           ; %state S for updating -- {u,0,{j}} -- dimension of that space is type+2
cost_p=1                   ; %cost of promoting a non manager to manager
cost_d=1                   ; %cost of demoting a manager to non manager
alpha_m=1              ; %Manager returns
alpha_n=0.1              ; %Non manager returns

typemin=1; %Lowest type
typemax=mubar ; %Highest type

amin=1; %Lowest productivity
amax=2; %Highest productivity

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


a_trans=1;
% q_trans=1;
u_trans=1;

%Transtion matrices (They are not perfect yet just to write it all down)
%A transition
aup=0.1;
adown=0.3;
astay=1-aup-adown;
a_trans=create_trans(adown,astay,aup,ats);

%Q transition
qup=0.2;
qdown=0.1;
qstay=1-qup-qdown;
q_trans=create_trans(qdown,qstay,qup,tpts);


%Unemp transition
ugain=0.00           ; %Probability unemployed move up
ustay=1-ulose-ugain  ; %Probability unemployed stay
u_trans=create_trans(ulose,ustay,ugain,tpts);

%Manager penalty toggle 
nm_penal=0;
%Convergenge speed
speed=1;

%Test with one dimension of everything
% a_trans=1;
% q_trans=1;
% u_trans=1;


% %Value Functions 
% zero_tpts=zeros(1,tpts);
% zero_tpts_tpts=zeros(tpts,tpts);
% zero_ats=zeros(1,ats);
% zero_tpts_ats=zeros(1,tpts,ats);
% zero_tpts_tpts_ats=zeros(tpts,tpts,ats);

% Uh=zeros(tpts); %Unemployed value function U(z)
% Veh=zeros(ats);  %Empty firm value function  Ve(a)
% Vmh=zeros(ats,tpts); %Firm with manager value function Vm(a,z)
% Vnh=zeros(ats,tpts); %Firm with no manager value function Vn(a,q)
% Vth=zeros(ats,tpts,tpts); %Firm with team value function Vt(a,z,q)

n=1; %mass of firms 


%Some Vectors for distributions
% e_udist=zeros(1,tpts); %Distribution of unemployed e_u(z)
% e_edist=zeros(1,ats); %Dist of empty firms e_e(a)
% e_mdist=zeros(ats,tpts); %Distribution of firms with manager e_m(a,z)
% e_ndist=zeros(ats,tpts); %Distribution of firms with no manager e_n(a,q)
% e_tdist=zeros(ats,tpts,tpts); %Distribution of firms with team e_t(a,z,q)

% %0.1 in each unemploted type
% e_udist=1*ones(1,tpts);
% u=sum(e_udist); %mass of unemployed

% % %We have (a)+2(a,z)x(a,z,q) type sof firm
% n_types=ats+2*(ats*tpts)+ats*tpts*tpts;
% %Split n into the n_types
% e_edist=(n/n_types)*ones(1,ats);
% e_mdist=(n/n_types)*ones(ats,tpts);
% e_ndist=(n/n_types)*ones(ats,tpts);
% e_tdist=(n/n_types)*ones(ats,tpts,tpts);

% %Need to guartee one of popoulation as well
% temp_u= 1 - sum(e_mdist,"all")+sum(e_ndist,"all")+2*sum(e_tdist,"all");
% e_udist=(temp_u/tpts)*ones(1,tpts);
% %Check
% n0=sum(e_edist)+sum(e_mdist,"all")+sum(e_ndist,"all")+sum(e_tdist,"all");
% %Split into
% %Population in the model 
% pop0=sum(e_udist)+sum(e_mdist,"all")+sum(e_ndist,"all")+2*sum(e_tdist,"all");




% %%%%Alternative, start with all empty firms and mass 1 of unemployed
e_udist=(1/tpts)*ones(1,tpts);
e_edist=(n/ats)*ones(1,ats);
e_mdist=zeros(ats,tpts); %Distribution of firms with manager e_m(a,z)
e_ndist=zeros(ats,tpts); %Distribution of firms with no manager e_n(a,q)
e_tdist=zeros(ats,tpts,tpts); %Distribution of firms with team e_t(a,z,q)

%Check
sum(e_edist)+sum(e_mdist,"all")+sum(e_ndist,"all")+sum(e_tdist,"all");

%Population in the model 
pop0=sum(e_udist)+sum(e_mdist,"all")+sum(e_ndist,"all")+2*sum(e_tdist,"all");


%For now lets tes with some random values for these distributions (values 0, .25 0.5 .75)
% e_udist= [0 0.25 0.5];
% %Mass of U
% u=sum(e_udist); %mass of unemployed
% e_edist= [0.25 0.5]; %Dist of empty firms   
% e_mdist(1,:)=[0 0.25 0.5];
% e_mdist(2,:)=[0.25 0.5 0.75];

% e_ndist(1,:)=[0 0.25 0.5];
% e_ndist(2,:)=[0.25 0.5 0.5];
% e_tdist(:,:,1)=[0 0.25 0.5; 0.25 0.5 0.75];
% e_tdist(:,:,2)=[0.25 0.5 0.75; 0.5 0.75 1];
% e_tdist(:,:,3)=[0.5 0.75 1; 0.75 1 1];



% %Some number to test it out 
% Veh=[0.5 1] ;
% Vmh(1,:)=[1 2 5];
% Vmh(2,:)=[4 5 2];
% % Vmh(:,:,3)=[7 8 9];

% Vnh(1,:)=[0 1 5];
% Vnh(2,:)=[1 2 2];
% % Vnh(:,:,3)=[2 3 3];

% Vth(:,:,1)=[0 1 1; 1 8 2];
% Vth(:,:,2)=[1 2 2; 2 5 3];
% Vth(:,:,3)=[2 5 3; 3 4 10];
% Uh=[0 2 3];





% % Initial VFs (test only)
% Ve=Veh;
% Vm=Vmh;
% Vn=Vnh;
% Vt=Vth;
% U=Uh;



%Inital guesses for value functions
Veini=zeros(1,ats);
Vmini=fman/(1-bt);
Vnini=fnman/(1-bt);
Vtini=fteam/(1-bt);
Uini=b/(1-bt);


%Value function iteration
[Ve, Vm, Vn, Vt, U, Veh, Vmh, Vnh, Vth, Vetl, Vmtl, Vntl, Vttl, Utl] = vf_iterationV2(e_edist,e_mdist,e_ndist,e_tdist,e_udist,Veini,Vmini,Vnini,Vtini,Uini,ats,tpts,cost_d,cost_p,nm_penal,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);

Uh=U;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Policy functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hiring policies, looking at origins not allocations yet 
[h_e_u, h_e_m, h_e_nm, h_e_t_m, h_e_t_nm, h_m_u, h_m_m, h_m_nm, h_m_t_m, h_m_t_nm, h_nm_u, h_nm_m, h_nm_nm, h_nm_t_m, h_nm_t_nm, h_t_u, h_t_m, h_t_nm, h_t_t_m, h_t_t_nm]...
=hire_policies(Vmh,Vnh,Veh,U,Vth,ats,tpts,nm_penal,cost_d,cost_p);

% Allocation Policies after hiring at SM
[p_e_m, p_e_n, p_m_m_u, p_m_m_d, p_m_n, p_n_n_u, p_n_n_p, p_n_m, p_t_m_u, p_t_m_d, p_t_n_u, p_t_n_p]...
=alloc_policies(Vmh,Vnh,U,Vth,ats,tpts,cost_d,cost_p);

%Reallocation Policies, before producing (using Vs not Vh)
% Notation here follows the paper, a bit wierd with previous notation in the code
[d_m, r_m, d_n, r_n, d_t_m, d_t_n, r_t_m, r_t_n, d_t_b]...
=reallocation_policies(Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p);



speed_dist=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distributions and LOMs 

% Distributions after the SM (implied by hiring and allocation policies)   
% e1 distributions
[e1_udist,e1_edist, e1_mdist, e1_ndist,e1_tdist] =e1_distV3(Veh,Vmh,Vnh,U,Vth,ats,tpts,cost_d,cost_p,e_udist,e_edist,e_mdist,e_ndist,e_tdist,lamu,lam,n,del);
%It is not really getting to one...
n1=sum(e1_edist)+sum(e1_mdist,"all")+sum(e1_ndist,"all")+sum(e1_tdist,"all"); % We are probabily missing a flow somewhere, that we will check later 
pop1=sum(e1_udist)+sum(e1_mdist,"all")+sum(e1_ndist,"all")+2*sum(e1_tdist,"all");

sum(e1_edist)-sum(e1_udist)-sum(e1_tdist,"all");

% %Fake update
% e_udist=(1-speed_dist)*e_udist+speed_dist*e1_udist;
% e_edist=(1-speed_dist)*e_edist+speed_dist*e1_edist;
% e_mdist=(1-speed_dist)*e_mdist+speed_dist*e1_mdist; 
% e_ndist=(1-speed_dist)*e_ndist+speed_dist*e1_ndist;
% e_tdist=(1-speed_dist)*e_tdist+speed_dist*e1_tdist;


%Lets assume these are right and focus on e1 
% Distributions before productions (reallocation policies)
%e2 distributions
[e2_udist,e2_edist,e2_mdist,e2_ndist,e2_tdist]=e2_dist(Ve,Vm,Vn,U,Vt,ats,tpts,cost_d,cost_p,e1_udist,e1_edist,e1_mdist,e1_ndist,e1_tdist);

%Check sum
n2=sum(e2_edist,"all")+sum(e2_mdist,"all")+sum(e2_ndist,"all")+sum(e2_tdist,"all"); %This is summing up to one if the fed in dist sums to 1 
pop2=sum(e2_udist)+sum(e2_mdist,"all")+sum(e2_ndist,"all")+2*sum(e2_tdist,"all");


%Disrtribtuions after shocks 
%e3 distributions
[e3_udist,e3_edist,e3_mdist,e3_ndist,e3_tdist]=e3_dist(ats,tpts,e2_udist,e2_edist,e2_mdist,e2_ndist,e2_tdist,u_trans,a_trans,q_trans);

%Check sum
n3=sum(e3_edist)+sum(e3_mdist,"all")+sum(e3_ndist,"all")+sum(e3_tdist,"all"); %This is summing up to one if the fed in dist sums to 1
pop3=sum(e3_udist)+sum(e3_mdist,"all")+sum(e3_ndist,"all")+2*sum(e3_tdist,"all");

%eplus Distributions
[eplus_udist,eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist]=eplus_dist(ats,tpts,e3_udist,e3_edist,e3_mdist,e3_ndist,e3_tdist,death,typebirth);

% Distributions After Entry and Exit 
%Check sum
nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all"); %This is summing up to one if the fed in dist sums to 1
popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");

% %Fake update
% e3_udist=(1-speed_dist)*e3_udist+speed_dist*eplus_udist;
% e3_edist=(1-speed_dist)*e3_edist+speed_dist*eplus_edist;
% e3_mdist=(1-speed_dist)*e3_mdist+speed_dist*eplus_mdist;
% e3_ndist=(1-speed_dist)*e3_ndist+speed_dist*eplus_ndist;
% e3_tdist=(1-speed_dist)*e3_tdist+speed_dist*eplus_tdist;

e_udist=(1-speed_dist)*e_udist+speed_dist*eplus_udist;
e_edist=(1-speed_dist)*e_edist+speed_dist*eplus_edist;
e_mdist=(1-speed_dist)*e_mdist+speed_dist*eplus_mdist;
e_ndist=(1-speed_dist)*e_ndist+speed_dist*eplus_ndist;
e_tdist=(1-speed_dist)*e_tdist+speed_dist*eplus_tdist;


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

% %Initial guesses for the LOM 
% eini_udist=0.1*ones(1,tpts);%These will be equal do the ini ones
% eini_edist=(n/n_types)*ones(1,ats);
% eini_mdist=(n/n_types)*ones(ats,tpts);
% eini_ndist=(n/n_types)*ones(ats,tpts);
% eini_tdist=(n/n_types)*ones(ats,tpts,tpts);
% pop0=sum(eini_udist)+sum(eini_mdist,"all")+sum(eini_ndist,"all")+2*sum(eini_tdist,"all");


%Alternative, start with all empty firms and mass 1 of unemployed
eini_udist=(1/tpts)*ones(1,tpts);
eini_edist=(n/ats)*ones(1,ats);
eini_mdist=zeros(ats,tpts); %Distribution of firms with manager e_m(a,z)
eini_ndist=zeros(ats,tpts); %Distribution of firms with no manager e_n(a,q)
eini_tdist=zeros(ats,tpts,tpts); %Distribution of firms with team e_t(a,z,q)
pop0=sum(eini_udist)+sum(eini_mdist,"all")+sum(eini_ndist,"all")+2*sum(eini_tdist,"all");


joint_it=0;
[e_udist, e_edist, e_mdist, e_ndist, e_tdist,store_n,store_p]=lom_iteration(eini_udist,eini_edist,eini_mdist,eini_ndist,eini_tdist,ats,tpts,Veh,Vmh,Vnh,U,Vth,Ve,Vm,Vn,Vt,cost_d,cost_p,nm_penal,lamu, lam, del, death, n,u_trans, a_trans,q_trans,typebirth,joint_it);
% %Need to fix the e1
% %Check sum
%Plot store_n and store_pop
figure 
plot(store_n)
title('Mass of Firms')

figure
plot(store_p)
title('Population')


n_post=sum(e_edist)+sum(e_mdist,"all")+sum(e_ndist,"all")+sum(e_tdist,"all");   
pop_post=sum(e_udist)+sum(e_mdist,"all")+sum(e_ndist,"all")+2*sum(e_tdist,"all");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Joint Loop
% Lets try the joint loop
%Initial guesses for the LOM
%Alternative, start with all empty firms and mass 1 of unemployed
eplus_udist=(1/tpts)*ones(1,tpts);
eplus_edist=(n/ats)*ones(1,ats);
eplus_mdist=zeros(ats,tpts); %Distribution of firms with manager e_m(a,z)
eplus_ndist=zeros(ats,tpts); %Distribution of firms with no manager e_n(a,q)
eplus_tdist=zeros(ats,tpts,tpts); %Distribution of firms with team e_t(a,z,q)

n0=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
pop0=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");


%VFs to initialize
%Inital guesses for value functions
Veini=zeros(1,ats);
Vmini=fman/(1-bt);
Vnini=fnman/(1-bt);
Vtini=fteam/(1-bt);
Uini=b/(1-bt);

[Ve, Vm, Vn, Vt, U, Veh, Vmh, Vnh, Vth, Vetl, Vmtl, Vntl, Vttl, Utl]= vf_iterationV2(eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist,eplus_udist,Veini,Vmini,Vnini,Vtini,Uini,ats,tpts,cost_d,cost_p,nm_penal,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);


update_speed_v=0.3;
update_speed=0.3;

diff_joint_lag1=0;
diff_joint_lag2=0;

%Iterate on value functions and type distribution
diff_joint    =100;
diff_joint_max=1e-4; %Value func/distribution max tolerance
  
it_joint      =0;
it_joint_min  =250;
it_joint_max  =2000;

store_n_outer=zeros(1,it_joint_min);
store_pop_outer=zeros(1,it_joint_min);
store_e_outer=zeros(1,it_joint_min);
store_u_outer=zeros(1,it_joint_min);
store_m_outer=zeros(1,it_joint_min);
store_nm_outer=zeros(1,it_joint_min);
store_t_outer=zeros(1,it_joint_min);


figure 
plot(store_n);
title('Mass of Firms');
hold on;
%Start timer
tic;
while (diff_joint>diff_joint_max | diff_joint_lag1>diff_joint_max  | diff_joint_lag2>diff_joint_max | it_joint<it_joint_min ) &&  it_joint<it_joint_max
    it_joint=it_joint+1;
    if it_joint==1
        e_udist=eplus_udist;
        e_edist=eplus_edist;
        e_mdist=eplus_mdist;
        e_ndist=eplus_ndist;
        e_tdist=eplus_tdist;
    end
    if it_joint>1
        e_udist=(1-update_speed)*e_udist+update_speed*eplus_udist;
        e_edist=(1-update_speed)*e_edist+update_speed*eplus_edist;
        e_mdist=(1-update_speed)*e_mdist+update_speed*eplus_mdist;
        e_ndist=(1-update_speed)*e_ndist+update_speed*eplus_ndist;
        e_tdist=(1-update_speed)*e_tdist+update_speed*eplus_tdist;
    end

    %Value function iteration
    [Veplus, Vmplus, Vnplus, Vtplus, Uplus, Vehplus, Vmhplus, Vnhplus, Vthplus, Vetlplus, Vmtlplus, Vntlplus, Vttlplus, Utlplus] = vf_iterationV2(e_edist,e_mdist,e_ndist,e_tdist,e_udist,Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p,nm_penal,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,update_speed_v);
    %Update the values
    Ve=(1-update_speed_v)*Ve+update_speed_v*Veplus;
    Vm=(1-update_speed_v)*Vm+update_speed_v*Vmplus;
    Vn=(1-update_speed_v)*Vn+update_speed_v*Vnplus;
    Vt=(1-update_speed_v)*Vt+update_speed_v*Vtplus;
    U=(1-update_speed_v)*U+update_speed_v*Uplus;
    Veh=(1-update_speed_v)*Veh+update_speed_v*Vehplus;
    Vmh=(1-update_speed_v)*Vmh+update_speed_v*Vmhplus;
    Vnh=(1-update_speed_v)*Vnh+update_speed_v*Vnhplus;
    Vth=(1-update_speed_v)*Vth+update_speed_v*Vthplus;



    %Keep the outer dist 
    eouter_udist=e_udist;
    eouter_edist=e_edist;
    eouter_mdist=e_mdist;
    eouter_ndist=e_ndist;
    eouter_tdist=e_tdist;

    %Iterate on masses of workers
    [eplus_udist,eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist,store_n,store_p]=lom_iteration(e_udist,e_edist,e_mdist,e_ndist,e_tdist,ats,tpts,Veh,Vmh,Vnh,U,Vth,Ve,Vm,Vn,Vt,cost_d,cost_p,nm_penal,lamu, lam, del, death, n,u_trans, a_trans,q_trans,typebirth,it_joint);
    nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
    popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");
    store_n_outer(it_joint)=nplus;
    store_pop_outer(it_joint)=popplus;
    store_e_outer(it_joint)=sum(eplus_edist);
    store_u_outer(it_joint)=sum(eplus_udist);
    store_m_outer(it_joint)=sum(eplus_mdist,"all");
    store_nm_outer(it_joint)=sum(eplus_ndist,"all");
    store_t_outer(it_joint)=sum(eplus_tdist,"all");


    %Compare Outer and Inner
    diff_joint_store=max([max(abs(eplus_udist-eouter_udist)),max(abs(eplus_edist-eouter_edist)),max(abs(eplus_mdist-eouter_mdist),[],[1 2]),max(abs(eplus_ndist-eouter_ndist),[],[1 2]),max(abs(eplus_tdist-eouter_tdist),[],[1 2 3])]);
    % diff_joint_lag2=diff_joint_lag1;
    % diff_joint_lag1=diff_joint;
    diff_joint=diff_joint_store;

    %Print every 100 iterations
    if mod(it_joint,10)==0
        fprintf('Joint Iteration %d, error %f \n', it_joint, diff_joint)
        plot(store_n, 'DisplayName', ['Iteration ' num2str(it_joint)]);
        drawnow;
    end
    if it_joint==it_joint_max
        fprintf(2,'Joint Failed to converge\n')
    end
    if diff_joint<diff_joint_max
        cprintf('green','Joint Converged in %d iterations\n',it_joint)
    end
    % if it_joint==46
    %     break
    % end
end
%End timer
toc




%Check sum of the eplus final distributions
nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");


%Plot store_n and store_pop
% figure
% plot(store_n_outer)
% title('Mass of Firms')

% figure
% plot(store_pop_outer)
% title('Population')

%Plot together
figure
plot(store_n_outer)
hold on
plot(store_pop_outer)
title('Mass of Firms and Population')
xlabel('Iterations')
legend('Mass of Firms','Population')

%Plot all together
figure
plot(store_e_outer)
hold on
% plot(store_pop_outer)
% plot(store_e_outer)
plot(store_u_outer)
plot(store_m_outer)
plot(store_nm_outer)
plot(store_t_outer)
title('Mass of Firms and Population')
xlabel('Iterations')
legend('Empty Firms','Unemployed','Managers','Non Managers','Teams')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Policy functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hiring policies, looking at origins not allocations yet 
[h_e_u, h_e_m, h_e_nm, h_e_t_m, h_e_t_nm, h_m_u, h_m_m, h_m_nm, h_m_t_m, h_m_t_nm, h_nm_u, h_nm_m, h_nm_nm, h_nm_t_m, h_nm_t_nm, h_t_u, h_t_m, h_t_nm, h_t_t_m, h_t_t_nm]...
=hire_policies(Vmh,Vnh,Veh,U,Vth,ats,tpts,nm_penal,cost_d,cost_p);

% Allocation Policies after hiring at SM
[p_e_m, p_e_n, p_m_m_u, p_m_m_d, p_m_n, p_n_n_u, p_n_n_p, p_n_m, p_t_m_u, p_t_m_d, p_t_n_u, p_t_n_p]...
=alloc_policies(Vmh,Vnh,U,Vth,ats,tpts,cost_d,cost_p);

%Reallocation Policies, before producing (using Vs not Vh)
% Notation here follows the paper, a bit wierd with previous notation in the code
[d_m, r_m, d_n, r_n, d_t_m, d_t_n, r_t_m, r_t_n, d_t_b]...
=reallocation_policies(Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p);

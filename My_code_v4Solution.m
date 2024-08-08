
% First attempt at computing the value functions for my model 


%Lets for now import some valeus from HLMP code
clear all
close all
cd('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit')
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
tpts =7                ; %type point space
ats=5                  ; %Productivity space 
spts =tpts+2           ; %state S for updating -- {u,0,{j}} -- dimension of that space is type+2
cost_p=0                   ; %cost of promoting a non manager to manager
cost_d=0                   ; %cost of demoting a manager to non manager
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
n=1; %mass of firms 


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
true=0;
%Convergenge speed
speed=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Joint loop
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

[Ve, Vm, Vn, Vt, U, Veh, Vmh, Vnh, Vth, Vetl, Vmtl, Vntl, Vttl, Utl]= vf_iterationV2(eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist,eplus_udist,Veini,Vmini,Vnini,Vtini,Uini,ats,tpts,cost_d,cost_p,true,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);


update_speed_v=1;
update_speed=1;

diff_joint_lag1=0;
diff_joint_lag2=0;

%Iterate on value functions and type distribution
diff_joint    =100;
diff_joint_max=1e-4; %Value func/distribution max tolerance
  
it_joint      =0;
it_joint_min  =250;
it_joint_max  =2000;

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
    [Veplus, Vmplus, Vnplus, Vtplus, Uplus, Vehplus, Vmhplus, Vnhplus, Vthplus, Vetlplus, Vmtlplus, Vntlplus, Vttlplus, Utlplus] = vf_iterationV2(e_edist,e_mdist,e_ndist,e_tdist,e_udist,Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p,true,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);
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
    [eplus_udist,eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist]=lom_iteration(e_udist,e_edist,e_mdist,e_ndist,e_tdist,ats,tpts,Veh,Vmh,Vnh,U,Vth,Ve,Vm,Vn,Vt,cost_d,cost_p,true,lamu, lam, del, death, n,u_trans, a_trans,q_trans,typebirth,it_joint);

    %Compare Outer and Inner
    diff_joint_store=max([max(abs(eplus_udist-eouter_udist)),max(abs(eplus_edist-eouter_edist)),max(abs(eplus_mdist-eouter_mdist),[],[1 2]),max(abs(eplus_ndist-eouter_ndist),[],[1 2]),max(abs(eplus_tdist-eouter_tdist),[],[1 2 3])]);
    diff_joint_lag2=diff_joint_lag1;
    diff_joint_lag1=diff_joint;
    diff_joint=diff_joint_store;

    %Print every 100 iterations
    if mod(it_joint,10)==0
        fprintf('Joint Iteration %d, error %f \n', it_joint, diff_joint)
    end
    if it_joint==it_joint_max
        fprintf(2,'Joint Failed to converge\n')
    end
    if diff_joint<diff_joint_max
        cprintf('green','Joint Converged in %d iterations\n',it_joint)
    end

end
%End timers
toc
    
%Check sum of the eplus final distributions
nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");

 save('model_solutions','Ve','Vm','Vn','Vt','U','Veh','Vmh','Vnh','Vth','eplus_udist','eplus_edist','eplus_mdist','eplus_ndist','eplus_tdist','nplus','popplus');


        
%% Alternative Specifications
%% Specification 1
%Team production function has q with return 1 and a with return 1
%But firm with no manager does not produce

fteam=zeros(ats,tpts,tpts);
fman=zeros(ats,tpts);
fnman=zeros(ats,tpts);
fe=zeros(1,ats);
for a=1:ats
    for z=1:tpts
        for q=1:tpts
            fteam(a,z,q)=A*a_type(a)*(type(z).^alpha_m)*(type(q));
        end
        fman(a,z)=A*a_type(a)*type(z).^alpha_m;
        fnman(a,z)=0;
    end
    fe(a)=0;
end

%Join the loop
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

[Ve, Vm, Vn, Vt, U, Veh, Vmh, Vnh, Vth, Vetl, Vmtl, Vntl, Vttl, Utl]= vf_iterationV2(eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist,eplus_udist,Veini,Vmini,Vnini,Vtini,Uini,ats,tpts,cost_d,cost_p,true,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);


update_speed_v=1;
update_speed=1;

diff_joint_lag1=0;
diff_joint_lag2=0;

%Iterate on value functions and type distribution
diff_joint    =100;
diff_joint_max=1e-4; %Value func/distribution max tolerance
  
it_joint      =0;
it_joint_min  =250;
it_joint_max  =2000;

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
    [Veplus, Vmplus, Vnplus, Vtplus, Uplus, Vehplus, Vmhplus, Vnhplus, Vthplus, Vetlplus, Vmtlplus, Vntlplus, Vttlplus, Utlplus] = vf_iterationV2(e_edist,e_mdist,e_ndist,e_tdist,e_udist,Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p,true,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);
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
    [eplus_udist,eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist]=lom_iteration(e_udist,e_edist,e_mdist,e_ndist,e_tdist,ats,tpts,Veh,Vmh,Vnh,U,Vth,Ve,Vm,Vn,Vt,cost_d,cost_p,true,lamu, lam, del, death, n,u_trans, a_trans,q_trans,typebirth,it_joint);

    %Compare Outer and Inner
    diff_joint_store=max([max(abs(eplus_udist-eouter_udist)),max(abs(eplus_edist-eouter_edist)),max(abs(eplus_mdist-eouter_mdist),[],[1 2]),max(abs(eplus_ndist-eouter_ndist),[],[1 2]),max(abs(eplus_tdist-eouter_tdist),[],[1 2 3])]);
    diff_joint_lag2=diff_joint_lag1;
    diff_joint_lag1=diff_joint;
    diff_joint=diff_joint_store;

    %Print every 100 iterations
    if mod(it_joint,10)==0
        fprintf('Joint Iteration %d, error %f \n', it_joint, diff_joint)
    end
    if it_joint==it_joint_max
        fprintf(2,'Joint Failed to converge\n')
    end
    if diff_joint<diff_joint_max
        cprintf('green','Sp1 Joint Converged in %d iterations\n',it_joint)
    end

end
%End timers
toc
    
%Check sum of the eplus final distributions
nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");

save('model_solutions_sp1','Ve','Vm','Vn','Vt','U','Veh','Vmh','Vnh','Vth','eplus_udist','eplus_edist','eplus_mdist','eplus_ndist','eplus_tdist','nplus','popplus');


%% Specification 2
%Only team produces, with q having a smaller return than z 


fteam=zeros(ats,tpts,tpts);
fman=zeros(ats,tpts);
fnman=zeros(ats,tpts);
fe=zeros(1,ats);
for a=1:ats
    for z=1:tpts
        for q=1:tpts
            fteam(a,z,q)=A*a_type(a)*(type(z).^alpha_m)*(type(q)^alpha_n);
        end
        fman(a,z)=0;
        fnman(a,z)=0;
    end
    fe(a)=0;
end



%Join the loop
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

[Ve, Vm, Vn, Vt, U, Veh, Vmh, Vnh, Vth, Vetl, Vmtl, Vntl, Vttl, Utl]= vf_iterationV2(eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist,eplus_udist,Veini,Vmini,Vnini,Vtini,Uini,ats,tpts,cost_d,cost_p,true,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);


update_speed_v=1;
update_speed=1;

diff_joint_lag1=0;
diff_joint_lag2=0;

%Iterate on value functions and type distribution
diff_joint    =100;
diff_joint_max=1e-4; %Value func/distribution max tolerance
  
it_joint      =0;
it_joint_min  =250;
it_joint_max  =2000;

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
    [Veplus, Vmplus, Vnplus, Vtplus, Uplus, Vehplus, Vmhplus, Vnhplus, Vthplus, Vetlplus, Vmtlplus, Vntlplus, Vttlplus, Utlplus] = vf_iterationV2(e_edist,e_mdist,e_ndist,e_tdist,e_udist,Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p,true,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);
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
    [eplus_udist,eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist]=lom_iteration(e_udist,e_edist,e_mdist,e_ndist,e_tdist,ats,tpts,Veh,Vmh,Vnh,U,Vth,Ve,Vm,Vn,Vt,cost_d,cost_p,true,lamu, lam, del, death, n,u_trans, a_trans,q_trans,typebirth,it_joint);

    %Compare Outer and Inner
    diff_joint_store=max([max(abs(eplus_udist-eouter_udist)),max(abs(eplus_edist-eouter_edist)),max(abs(eplus_mdist-eouter_mdist),[],[1 2]),max(abs(eplus_ndist-eouter_ndist),[],[1 2]),max(abs(eplus_tdist-eouter_tdist),[],[1 2 3])]);
    diff_joint_lag2=diff_joint_lag1;
    diff_joint_lag1=diff_joint;
    diff_joint=diff_joint_store;

    %Print every 100 iterations
    if mod(it_joint,10)==0
        fprintf('Joint Iteration %d, error %f \n', it_joint, diff_joint)
    end
    if it_joint==it_joint_max
        fprintf(2,'Joint Failed to converge\n')
    end
    if diff_joint<diff_joint_max
        cprintf('green','Sp2 Joint Converged in %d iterations\n',it_joint)
    end

end
%End timers
toc
    
%Check sum of the eplus final distributions
nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");

 save('model_solutions_sp2','Ve','Vm','Vn','Vt','U','Veh','Vmh','Vnh','Vth','eplus_udist','eplus_edist','eplus_mdist','eplus_ndist','eplus_tdist','nplus','popplus');


 %% Specification 3
 % Sp1 but with higher transition to upper shocks in a and lower learning

fteam=zeros(ats,tpts,tpts);
fman=zeros(ats,tpts);
fnman=zeros(ats,tpts);
fe=zeros(1,ats);
for a=1:ats
    for z=1:tpts
        for q=1:tpts
            fteam(a,z,q)=A*a_type(a)*(type(z).^alpha_m)*(type(q));
        end
        fman(a,z)=A*a_type(a)*type(z).^alpha_m;
        fnman(a,z)=0;
    end
    fe(a)=0;
end

%A transition
aup=0.3;
adown=0.1;
astay=1-aup-adown;
a_trans=create_trans(adown,astay,aup,ats);

%Q transition
qup=0.1;
qdown=0.1;
qstay=1-qup-qdown;
q_trans=create_trans(qdown,qstay,qup,tpts);


%Unemp transition
ugain=0.00           ; %Probability unemployed move up
ustay=1-ulose-ugain  ; %Probability unemployed stay
u_trans=create_trans(ulose,ustay,ugain,tpts);


%Join the loop
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

[Ve, Vm, Vn, Vt, U, Veh, Vmh, Vnh, Vth, Vetl, Vmtl, Vntl, Vttl, Utl]= vf_iterationV2(eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist,eplus_udist,Veini,Vmini,Vnini,Vtini,Uini,ats,tpts,cost_d,cost_p,true,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);


update_speed_v=1;
update_speed=1;

diff_joint_lag1=0;
diff_joint_lag2=0;

%Iterate on value functions and type distribution
diff_joint    =100;
diff_joint_max=1e-4; %Value func/distribution max tolerance
  
it_joint      =0;
it_joint_min  =250;
it_joint_max  =2000;

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
    [Veplus, Vmplus, Vnplus, Vtplus, Uplus, Vehplus, Vmhplus, Vnhplus, Vthplus, Vetlplus, Vmtlplus, Vntlplus, Vttlplus, Utlplus] = vf_iterationV2(e_edist,e_mdist,e_ndist,e_tdist,e_udist,Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p,true,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);
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
    [eplus_udist,eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist]=lom_iteration(e_udist,e_edist,e_mdist,e_ndist,e_tdist,ats,tpts,Veh,Vmh,Vnh,U,Vth,Ve,Vm,Vn,Vt,cost_d,cost_p,true,lamu, lam, del, death, n,u_trans, a_trans,q_trans,typebirth,it_joint);

    %Compare Outer and Inner
    diff_joint_store=max([max(abs(eplus_udist-eouter_udist)),max(abs(eplus_edist-eouter_edist)),max(abs(eplus_mdist-eouter_mdist),[],[1 2]),max(abs(eplus_ndist-eouter_ndist),[],[1 2]),max(abs(eplus_tdist-eouter_tdist),[],[1 2 3])]);
    diff_joint_lag2=diff_joint_lag1;
    diff_joint_lag1=diff_joint;
    diff_joint=diff_joint_store;

    %Print every 100 iterations
    if mod(it_joint,10)==0
        fprintf('Joint Iteration %d, error %f \n', it_joint, diff_joint)
    end
    if it_joint==it_joint_max
        fprintf(2,'Joint Failed to converge\n')
    end
    if diff_joint<diff_joint_max
        cprintf('green','Sp3 Joint Converged in %d iterations\n',it_joint)
    end

end
%End timers
toc
    
%Check sum of the eplus final distributions
nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");

save('model_solutions_sp3','Ve','Vm','Vn','Vt','U','Veh','Vmh','Vnh','Vth','eplus_udist','eplus_edist','eplus_mdist','eplus_ndist','eplus_tdist','nplus','popplus');


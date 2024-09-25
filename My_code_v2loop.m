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
% typebirth=zeros(1,tpts);
% for i=1:mnew_high
%     typebirth(i)= mnew*exp(-type(i)*mnew)./sum(mnew*exp(-type(1:mnew_high)*mnew)); %Birth rate
% end


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
nm_penal=1;



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
e_udist=zeros(1,tpts); %Distribution of unemployed e_u(z)
e_edist=zeros(1,ats); %Dist of empty firms e_e(a)
e_mdist=zeros(ats,tpts); %Distribution of firms with manager e_m(a,z)
e_ndist=zeros(ats,tpts); %Distribution of firms with no manager e_n(a,q)
e_tdist=zeros(ats,tpts,tpts); %Distribution of firms with team e_t(a,z,q)

% 0.1 in each unemploted type
e_udist=0.1*ones(1,tpts);

%We have (a)+2(a,z)x(a,z,q) type sof firm
n_types=ats+2*(ats*tpts)+ats*tpts*tpts;
%Split n into the n_types
e_edist=(n/n_types)*ones(1,ats);
e_mdist=(n/n_types)*ones(ats,tpts);
e_ndist=(n/n_types)*ones(ats,tpts);
e_tdist=(n/n_types)*ones(ats,tpts,tpts);
%Check
sum(e_edist)+sum(e_mdist,"all")+sum(e_ndist,"all")+sum(e_tdist,"all");



% %For now lets tes with some random values for these distributions (values 0, .25 0.5 .75)
% e_udist= [0 0.25 0.5];
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



%Use the function
[Ve, Vm, Vn, Vt, U, Veh, Vmh, Vnh, Vth, Vetl, Vmtl, Vntl, Vttl, Utl] = vf_iteration(e_edist,e_mdist,e_ndist,e_tdist,e_udist,n,Veini,Vmini,Vnini,Vtini,Uini,ats,tpts,cost_d,cost_p,nm_penal,lamu, lam, del, bt, death, bpf, bpw, nfirm, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Here is the start of what will be the function to compute the value functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Use initial guesses for value functions
% Veup = Veini;
% Vmup = Vmini;
% Vnup = Vnini;
% Vtup = Vtini;
% Uup = Uini;

% %Mass of U
% u=sum(e_udist); %mass of unemployed

% %Iteration parameters
% diff=100;
% diffmax=.5e-4;
 
% it=0; 
% itmax=100000;

% speed=1; %Speed of convergence
 
% while diff>diffmax & it<itmax
%  it=it+1;
%  if it==1
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Update stage
%     Ve=Veup;
%     Vm=Vmup;
%     Vn=Vnup;
%     Vt=Vtup;
%     U=Uup;
%  end

%  if it>1
%     Ve=(1-speed)*Ve + speed*Veup;
%     Vm= (1-speed)*Vm + speed*Vmup;
%     Vn= (1-speed)*Vn + speed*Vnup;
%     Vt= (1-speed)*Vt + speed*Vtup;
%     U= (1-speed)*U + speed*Uup;

    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Fire/Promotion VFs
%     %Empty firm
%     Veh=Ve;
    
%     %Firm with manager (a,z)
%     Vm_alloc(:,:,1)=Vm; %Keep as it is
%     Vm_alloc(:,:,2)=Vn-cost_d; %Demote z
%     Vm_alloc(:,:,3)=repmat(Ve',1, tpts)+repmat(U,ats,1);%Fire z
%     Vmh=max(Vm_alloc,[],3);
    
    
%     %Firm with non-manager (a,q) 
%     Vn_alloc(:,:,1)=Vn; %Keep as it is
%     Vn_alloc(:,:,2)=Vm-cost_p; %Promote q
%     Vn_alloc(:,:,3)=repmat(Ve',1, tpts)+repmat(U,ats,1); %Fire q
%     Vnh=max(Vn_alloc,[],3);
    
%     %Firm with team (a,z,q)
%     Vt_alloc(:,:,:,1)=Vt; %Keep as it is
%     Vt_alloc(:,:,:,2)=permute(repmat(Vm,1,1,tpts),[1 3 2])+repmat(U,ats,1,tpts)-cost_p; %Fire z and promote q
%     Vt_alloc(:,:,:,3)=repmat(Vn,1,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]) - cost_d; %Fire q and demote z
%     Vt_alloc(:,:,:,4)=repmat(Vm,1,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]); %Fire q 
%     Vt_alloc(:,:,:,5)=permute(repmat(Vn,1,1,tpts),[1 3 2])+repmat(U,ats,1,tpts); %Fire z
%     Vt_alloc(:,:,:,6)=permute(Vt, [1 3 2])-cost_d-cost_p; %Swap z and q with promotion and demotion costs
%     Vt_alloc(:,:,:,7)=permute(repmat(Ve',1,tpts,tpts),[1 3 2])+repmat(U,ats,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]); %Fire z and q
%     Vth=max(Vt_alloc,[],4);
    
%     %Unemployed
%     Uh=U;
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%Empty firm continuation value Ve(a)
%     %Gains from trade 
%     [gt_ef_meet_u, gt_ef_meet_m, gt_ef_meet_nm, gt_ef_meet_t_m, gt_ef_meet_t_nm, gt_ef_meet_t]=ef_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,nm_penal);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% Firm with manager continuation value Vm(a,z)
%     %Gains from trade
%     [gt_em_meet_u, gt_em_meet_m, gt_em_meet_nm, gt_em_meet_t_m, gt_em_meet_t_nm, gt_em_meet_t]=mf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,nm_penal,cost_d);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% Firm with no manager continuation value Vn(a,q)
%     %Gains from trade
%     [gt_en_meet_u, gt_en_meet_m, gt_en_meet_nm, gt_en_meet_t_m, gt_en_meet_t_nm, gt_en_meet_t]=nf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,nm_penal,cost_p);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% Firm with team continuation value Vt(a,z,q)
%     %Gains from trade
%     [gt_tf_meet_u, gt_tf_meet_m, gt_tf_meet_nm, gt_tf_meet_t_m, gt_tf_meet_t_nm, gt_tf_meet_t]=tf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,nm_penal,cost_p,cost_d);
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Expected Gains of poaching 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Expected gains from trade for empty firm
%     %Empty firm (a) poaching from unemplyment
%     Egt_e_u=e_udist*gt_ef_meet_u; %From unemployment to empty firm
    
%     %Empty firm (a) poaching from firm with manager
%     Egt_e_m=reshape(sum(e_mdist.*gt_ef_meet_m, [1 2]),[1 ats]); %From manager to empty firm
    
%     %Empty firm (a) poaching from firm with non manager
%     Egt_e_n=reshape(sum(e_ndist.*gt_ef_meet_nm, [1 2]), [1 ats]);
    
%     %Empty firm (a) poaching from firm with team
%     Egt_e_t=reshape(sum(e_tdist.*gt_ef_meet_t, [1 2 3]), [1 ats]);
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Expected Gains of poaching for firm with manager
%     %Firm (a,z) Poaching from unemployment
%     Egt_m_u=reshape(permute(sum(e_udist.*permute(gt_em_meet_u, [2 1 3]),[2]),[2 1 3]),ats,tpts); %From unemployment to firm with manager
    
%     %Firm (a,z) Poaching from firm with manager
%     Egt_m_m= reshape(sum(e_mdist.*gt_em_meet_m, [1 2]), ats, tpts); %From firm with manager manager to firm with manager
    
%     %Firm (a,z) Poaching from firm with no manager
%     Egt_m_n=reshape(sum(e_ndist.*gt_em_meet_nm, [1 2]), ats,tpts); %From firm with no manager to firm with manager
    
%     %Firm (a,z) Poaching from firm with team
%     Egt_m_t=reshape(sum(e_tdist.*gt_em_meet_t, [1 2 3]), ats, tpts); %From firm with team to firm with manager
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Expected Gains of poaching for firm with no manager
%     %Firm (a,q) Poaching from unemployment
%     Egt_n_u=reshape(permute(sum(e_udist.*permute(gt_en_meet_u, [2 1 3]),[2]),[2 1 3]),ats,tpts); %From unemployment to firm with no manager
    
%     %Firm (a,q) Poaching from firm with manager
%     Egt_n_m= reshape(sum(e_mdist.*gt_en_meet_m, [1 2]), ats, tpts); %From firm with manager to firm with no manager
    
%     %Firm (a,q) Poaching from firm with no manager
%     Egt_n_n=reshape(sum(e_ndist.*gt_en_meet_nm, [1 2]), ats,tpts); %From firm with no manager to firm with no manager
    
%     %Firm (a,q) Poaching from firm with team
%     Egt_n_t=reshape(sum(e_tdist.*gt_en_meet_t, [1 2 3]), ats, tpts); %From firm with team to firm with no manager
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Expected Gains of poaching for firm with team
%     %Firm (a,z,q) Poaching from unemployment
%     Egt_t_u=reshape(permute(sum(e_udist.*permute(gt_tf_meet_u, [2 1 3 4]),[2]),[2 1 3 4]),ats,tpts,tpts); %From unemployment to firm with team
    
%     %Firm (a,z,q) Poaching from firm with manager
%     Egt_t_m= reshape(sum(e_mdist.*gt_tf_meet_m, [1 2]), ats, tpts, tpts); %From firm with manager to firm with team
    
%     %Firm (a,z,q) Poaching from firm with no manager
%     Egt_t_n=reshape(sum(e_ndist.*gt_tf_meet_nm, [1 2]), ats,tpts, tpts); %From firm with no manager to firm with team
    
%     %Firm (a,z,q) Poaching from firm with team
%     Egt_t_t=reshape(sum(e_tdist.*gt_tf_meet_t, [1 2 3]), ats, tpts, tpts); %From firm with team to firm with team
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Expected Gains of being poached
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Expected gains of being poached for unemployed worker 
%     %Unemp z Poached from empty firm
%     Epgt_u_e=e_edist*permute(gt_ef_meet_u, [2 1]); %From unemployment to empty firm
    
%     %Unemp z Poached from firm with manager
%     Epgt_u_m=reshape(sum(e_mdist.*permute(gt_em_meet_u, [2 3 1]),[1 2]),1,tpts); %From unemployment to firm with manager
    
%     %Unemp z Poached from firm with no manager
%     Epgt_u_n=reshape(sum(e_ndist.*permute(gt_en_meet_u, [2 3 1]), [1 2]),1,tpts); %From unemployment to firm with no manager
    
%     %Unemp z Poached from firm with team
%     Epgt_u_t=reshape(sum(e_tdist.*permute(gt_tf_meet_u, [2 3 4 1]),[1 2 3]),1,tpts); %From unemployment to firm with team
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Expected Gains of being poached for firm with manager
%     %Firm (a,z) Poached from empty firm
%     Epgt_m_e=sum(e_edist'.*(gt_ef_meet_m),3); %From manager to empty firm 
    
%     %Firm (a,z) Poached from firm with manager
%     Epgt_m_m=reshape(sum(e_mdist.*permute(gt_em_meet_m, [3 4 1 2]),[1 2]),ats,tpts); %From manager to manager
    
%     %Firm (a,z) Poached from firm with no manager
%     Epgt_m_n=reshape(sum(e_ndist.*permute(gt_en_meet_m, [3 4 1 2]), [1 2]),ats,tpts); %From manager to no manager
    
%     %Firm (a,z) Poached from firm with team
%     Epgt_m_t=reshape(sum(e_tdist.*permute(gt_tf_meet_m, [3 4 5 1 2]),[1 2 3]),ats,tpts); %From manager to team
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Expected Gains of being poached for firm with no manager
%     %Firm (a,q) Poached from empty firm
%     Epgt_n_e=sum(e_edist'.*(gt_ef_meet_nm),3); %From no manager to empty firm
    
%     %Firm (a,q) Poached from firm with manager
%     Epgt_n_m=reshape(sum(e_mdist.*permute(gt_em_meet_nm, [3 4 1 2]),[1 2]),ats,tpts); %From no manager to manager
    
%     %Firm (a,q) Poached from firm with no manager
%     Epgt_n_n=reshape(sum(e_ndist.*permute(gt_en_meet_nm, [3 4 1 2]), [1 2]),ats,tpts); %From no manager to no manager
    
%     %Firm (a,q) Poached from firm with team
%     Epgt_n_t=reshape(sum(e_tdist.*permute(gt_tf_meet_nm, [3 4 5 1 2]),[1 2 3]),ats,tpts); %From no manager to team
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Expected Gains of being poached for firm with team
%     %Firm (a,z,q) Poached from empty firm
%     Epgt_t_e=sum(e_edist'.*(gt_ef_meet_t),3); %From team to empty firm
    
%     %Firm (a,z,q) Poached from firm with manager
%     Epgt_t_m=reshape(sum(e_mdist.*permute(gt_em_meet_t, [4 5 1 2 3]),[1 2]),ats,tpts, tpts); %From team to manager   
    
%     %Firm (a,z,q) Poached from firm with no manager
%     Epgt_t_n=reshape(sum(e_ndist.*permute(gt_en_meet_t, [4 5 1 2 3]), [1 2]),ats,tpts, tpts); %From team to no manager
    
%     %Firm (a,z,q) Poached from firm with team
%     Epgt_t_t=reshape(sum(e_tdist.*permute(gt_tf_meet_t, [4 5 6 1 2 3]),[1 2 3]),ats,tpts, tpts); %From team to team
    
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % VF after the S&M
    
%     % Empty firm
%     Vetl=Veh+ (lamu*bpf/u)*Egt_e_u + (lam*bpf/n)*Egt_e_m + (lam*bpf/n)*Egt_e_n + (lam*bpf/n)*Egt_e_t;
    
%     % Firm with manager
%     Vmtl=Vmh+ del*(repmat(Veh',1, tpts)+repmat(Uh,ats,1)- Vmh) + (lamu*bpf/u)*Egt_m_u + (lam*bpf/n)*Egt_m_m + (lam*bpf/n)*Egt_m_n + (lam*bpf/n)*Egt_m_t...
%         + (lam*bpw/n)*(Epgt_m_e + Epgt_m_m + Epgt_m_n + Epgt_m_t);
    
%     % Firm with no manager
%     Vntl=Vnh+ del*(repmat(Veh',1, tpts)+repmat(Uh,ats,1)- Vnh) + (lamu*bpf/u)*Egt_n_u + (lam*bpf/n)*Egt_n_m + (lam*bpf/n)*Egt_n_n + (lam*bpf/n)*Egt_n_t...
%         + (lam*bpw/n)*(Epgt_n_e + Epgt_n_m + Epgt_n_n + Epgt_n_t);
    
%     %Firm with team
%     Vttl=Vth+ del*(repmat(Vmh,1,1,tpts)+repmat(Uh,ats,1,tpts)-Vth) + del*(repmat(Vnh,1,1,tpts)+repmat(Uh,ats,1,tpts)-Vth) + (lamu*bpf/u)*Egt_t_u + (lam*bpf/n)*Egt_t_m + (lam*bpf/n)*Egt_t_n + (lam*bpf/n)*Egt_t_t...
%         + (lam*bpw/n)*(Epgt_t_e + Epgt_t_m + Epgt_t_n + Epgt_t_t);
    
%     %Unemployed
%     Utl= (lamu*bpw/u)*(Epgt_u_e + Epgt_u_m + Epgt_u_n + Epgt_u_t);  %Not dead sure this is right
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Vup - Productions and shocks VFs
    
%     %Empty firm
%     Veup= fe + bt*(a_trans*Vetl')';
    
%     %Firm with manager
%     Vmup=fman + bt*death*repmat(a_trans*Vetl',1,tpts)+ bt*(1-death)*(a_trans*Vmtl);
    
%     %Firm with no manager (This is q-trans do not depending on a. If it depends will have to include pages)
%     Vnup=fnman+bt*death*repmat(a_trans*Vetl',1,tpts)+ bt*(1-death)*(q_trans*(a_trans*Vntl)')';
    
%     %Firm with team
%     %Transitions for full firm
%     %Conditional on q, matrix mult to get the expected value for (a,z). Do that for every q 
%     Va_prime=zeros(ats,tpts,tpts);
%     Va_prime(:,:,1)=a_trans*Vttl(:,:,1);
%     Va_prime(:,:,2)=a_trans*Vttl(:,:,2);
%     Va_prime(:,:,3)= a_trans*Vttl(:,:,3);
%     %Permute to multiply it by q_trans
%     Va_prime=permute(Va_prime, [3 2 1]);
    
%     %Conditional on a, matrix mult to get the expected value for (z,q). Do that for every a
%     Vq_prime=zeros(tpts,tpts,ats);
%     Vq_prime(:,:,1)=q_trans*Va_prime(:,:,1);
%     Vq_prime(:,:,2)=q_trans*Va_prime(:,:,2);
%     Vqa_prime=permute(Vq_prime, [3 2 1]);
    
%     Vtup=fteam+ bt*death.^2*repmat(a_trans*Vetl',1,tpts,tpts)+ bt*death*(1-death)*repmat((a_trans*Vmtl),1,1,tpts)...
%         + bt*death*(1-death)*permute(repmat((q_trans*(a_trans*Vntl)')',1,1,tpts),[1 3 2])...
%         + bt*(1-death).^2*Vqa_prime;
    
%     %Unemployed worker
%     Uplus= Uh + death*(-Uh) + (1-death)*Utl;
%     Uup=  b + bt*(u_trans*Uplus')';


%     %diff=max([max(abs(V0 - V0up)),max(max(abs(V1 - V1up))),max(max(abs(V2 - V2up))),max(abs(U-Uup))]);
%     diff=max([max(abs(Ve-Veup)), max(abs(Vm-Vmup),[],[1 2]), max(abs(Vn-Vnup),[],[1 2]), max(abs(Vt-Vtup),[],[1 2 3]), max(abs(U-Uup))]);  

%     %Print every 50 iterations the difference and some text
%     if mod(it,100)==0
%         fprintf('Iteration %d, error %f \n', it, diff)
%     end
%     %In red failed to converge
%     if it==itmax
%         fprintf('Failed to converge\n')
%     end
%     if diff<diffmax
%         fprintf('Converged in %d iterations \n',it)
%     end
% end
% end 

% Ve=Veup;
% Vm=Vmup;
% Vn=Vnup;
% Vt=Vtup;
% U=Uup;

% %Export the hats as well    %Fire/Promotion VFs
% % %Empty firm
% Veh=Ve;    
% %Firm with manager (a,z)
% Vm_alloc(:,:,1)=Vm; %Keep as it is
% Vm_alloc(:,:,2)=Vn-cost_d; %Demote z
% Vm_alloc(:,:,3)=repmat(Ve',1, tpts)+repmat(U,ats,1);%Fire z
% Vmh=max(Vm_alloc,[],3);


% %Firm with non-manager (a,q) 
% Vn_alloc(:,:,1)=Vn; %Keep as it is
% Vn_alloc(:,:,2)=Vm-cost_p; %Promote q
% Vn_alloc(:,:,3)=repmat(Ve',1, tpts)+repmat(U,ats,1); %Fire q
% Vnh=max(Vn_alloc,[],3);

% %Firm with team (a,z,q)
% Vt_alloc(:,:,:,1)=Vt; %Keep as it is
% Vt_alloc(:,:,:,2)=permute(repmat(Vm,1,1,tpts),[1 3 2])+repmat(U,ats,1,tpts)-cost_p; %Fire z and promote q
% Vt_alloc(:,:,:,3)=repmat(Vn,1,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]) - cost_d; %Fire q and demote z
% Vt_alloc(:,:,:,4)=repmat(Vm,1,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]); %Fire q 
% Vt_alloc(:,:,:,5)=permute(repmat(Vn,1,1,tpts),[1 3 2])+repmat(U,ats,1,tpts); %Fire z
% Vt_alloc(:,:,:,6)=permute(Vt, [1 3 2])-cost_d-cost_p; %Swap z and q with promotion and demotion costs
% Vt_alloc(:,:,:,7)=permute(repmat(Ve',1,tpts,tpts),[1 3 2])+repmat(U,ats,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]); %Fire z and q
% Vth=max(Vt_alloc,[],4);

% %Unemployed
% Uh=U;


    


     
    
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

type=typemin:(typemax-typemin)/(tpts-1):typemax; %Types
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
nm_penal=0               ; %Manager penalty toggle
speed=1              ; %Convergenge speed






%Load some results that we will simullate over
load('model_solutions_sp3'+location)
load('wages_a'+string(ats)+'_z'+string(tpts))


%Lets create classes
%Parameters
p=struct('lamu',lamu,'lam',lam,'fcomp',fcomp,'del',del,'death',death,...
'ulose',ulose,'mnew',mnew,'bpf',bpf,'bpw',bpw, ...
'homeprod',homeprod,'A',A,'mnew_high',mnew_high,'typemin',typemin,'typemax',typemax,'amin',amin,'amax',amax,'type',type,'a_type',a_type,'typebirth',typebirth,'pnew',pnew,'pold',pold,...
'fteam',fteam,'fman',fman,'fnman',fnman,'fe',fe,'b',b,'a_trans',a_trans,'q_trans',q_trans,'u_trans',u_trans,'cost_p',cost_p,'cost_d',cost_d,'alpha_m',alpha_m,'alpha_n',alpha_n,'n',n,'tpts',tpts,'ats',ats);

%Value Functions
v=struct('Ve',Ve,'Vm',Vm,'Vn',Vn,'Vt',Vt,'U',U,'Veh',Veh,'Vmh',Vmh,'Vnh',Vnh,'Vth',Vth);

%Distributions
d=struct('eplus_udist',eplus_udist,'eplus_edist',eplus_edist,'eplus_mdist',eplus_mdist,'eplus_ndist',eplus_ndist,'eplus_tdist',eplus_tdist);

%Wage Frunctiions
w=struct('Wm',Wm,'Wn',Wn,'Wtm',Wtm,'Wtn',Wtn,'wgrid',wgrid,'Wmh',Wmh,'Wnh',Wnh,'Wtmh',Wtmh,'Wtnh',Wtnh);

%Simulation as a function
n_months=240;
n_firms=200000;
n_workers=n_firms;
seed=1;
%For firms
tic
[firm_status,a_ftp,m_ftp,n_ftp,m_fwage,n_fwage,hire_m,hire_n,fire_m,fire_n,prom_sm,prom_realloc]=SimulateFirm(n_months,n_firms,seed,ats,tpts,p,v,d,w);
toc

fs=struct('firm_status',firm_status,'a_ftp',a_ftp,'m_ftp',m_ftp,'n_ftp',n_ftp,'m_fwage',m_fwage,'n_fwage',n_fwage,'hire_m',hire_m,'hire_n',hire_n,'fire_m',fire_m,'fire_n',fire_n,'prom_sm',prom_sm,'prom_realloc',prom_realloc);

%For workers
tic
[worker_status,self_z,co_z,a_wtp,self_wage,co_wage,self_age,self_tenure,co_tenure,hire_m,hire_n,fire_m,fire_n,prom_sm,prom_realloc,fire_m_realloc,fire_n_realloc]=...
    SimulateWorker(n_months,n_workers,seed,ats,tpts,p,v,d,w,fs);
toc
%Whole simulation is taking about 6 minutes, which is not bad.

% if any(firm_status==0)
%     disp('There are mistakes')
% else
%     disp('Everything is fine')
% end




%% Simulations original codes

% %Policy functions
% %Get the policy functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Policy functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Hiring policies, looking at origins not allocations yet 
% [h_e_u, h_e_m, h_e_nm, h_e_t_m, h_e_t_nm, h_m_u, h_m_m, h_m_nm, h_m_t_m, h_m_t_nm, h_nm_u, h_nm_m, h_nm_nm, h_nm_t_m, h_nm_t_nm, h_t_u, h_t_m, h_t_nm, h_t_t_m, h_t_t_nm]...
%     =hire_policies(Vmh,Vnh,Veh,U,Vth,ats,tpts,0,cost_d,cost_p);
    
% % Allocation Policies after hiring at SM
% [p_e_m, p_e_n, p_m_m_u, p_m_m_d, p_m_n, p_n_n_u, p_n_n_p, p_n_m, p_t_m_u, p_t_m_d, p_t_n_u, p_t_n_p]...
%     =alloc_policies(Vmh,Vnh,U,Vth,ats,tpts,cost_d,cost_p);

% %Reallocation policies
% [d_m, r_m, d_n, r_n, d_t_m, d_t_n, r_t_m, r_t_n, d_t_b]...
%     =reallocation_policies(Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p);



% e_udist=eplus_udist;
% e_edist=eplus_edist;
% e_mdist=eplus_mdist;
% e_ndist=eplus_ndist;
% e_tdist=eplus_tdist;

% % Distributions 
% %Firm type (unweighted)
% firm_dist=[sum(e_edist),sum(e_mdist,"all"),sum(e_ndist,"all"),sum(e_tdist,"all")]/n;
% firm_dist_pdf=firm_dist/sum(firm_dist);

% %CDF os mass of each type of firm 
% cdf_firm_dist=zeros(1,4);
% for i=1:4
%     cdf_firm_dist(i)=sum(firm_dist(1:i));
% end
% cdf_firm_dist(end)=1; %Just to make sure it is 1

% %Employee type wighted for the arival rates. This will be useful to get the probability if a match
% emp_dist=[lamu*sum(e_udist),lam*sum(e_mdist,"all"),lam*sum(e_ndist,"all"),lam*sum(e_tdist,"all")];
% emp_dist_pdf=emp_dist/sum(emp_dist);

% %Prob of matching 
% prob_matching=sum(emp_dist);

% % Now we do a cdf, considering the effective mass of each type and relative to the probability of matching
% cdf_emp_dist=zeros(1,4);
% for i=1:4
%     cdf_emp_dist(i)=sum(emp_dist(1:i))/prob_matching;
% end
% cdf_emp_dist(end)=1; %Just to make sure it is 1

% % Tranform measures into cdfs (multi-dimensional)
% cdf_e_udist=measure_to_cdf(e_udist);
% cdf_e_edist=measure_to_cdf(e_edist);
% cdf_e_mdist=measure2d_to_cdf(e_mdist); %Will have to be careful here when making the draw. Do I need to map back to the original coordinates?
% cdf_e_ndist=measure2d_to_cdf(e_ndist);
% cdf_e_tdist=measure3d_to_cdf(e_tdist);


% %% Adjust probabilities of shocks transitions
% % For a
% %Given a, vector of probabilities to go down, stay or up
% down_a=zeros(1,ats);
% stay_a=zeros(1,ats);
% up_a=zeros(1,ats);
% for a=1:ats
%     if a==1
%         down_a(a)=0;
%     else
%         down_a(a)=a_trans(a,a-1);
%     end
%     stay_a(a)=a_trans(a,a);
%     if a==ats
%         up_a(a)=0;
%     else
%         up_a(a)=a_trans(a,a+1);
%     end
% end

% % For q
% down_q=zeros(ats,tpts);
% stay_q=zeros(ats,tpts);
% up_q=zeros(ats,tpts);
% for a=1:ats
%     for q=1:tpts
%         if q==1
%             down_q(a,q)=0;
%         else
%             down_q(a,q)=q_trans(q,q-1,a);
%         end
%         stay_q(a,q)=q_trans(q,q,a);
%         if q==tpts
%             up_q(a,q)=0;
%         else
%             up_q(a,q)=q_trans(q,q+1,a);
%         end
%     end
% end

% % For u
% down_u=zeros(1,tpts);
% stay_u=zeros(1,tpts);
% up_u=zeros(1,tpts);
% for u=1:tpts
%     if u==1
%         down_u(u)=0;
%     else
%         down_u(u)=u_trans(u,u-1);
%     end
%     stay_u(u)=u_trans(u,u);
%     if u==tpts
%         up_u(u)=0;
%     else
%         up_u(u)=u_trans(u,u+1);
%     end
% end


% tic
% %% Simulation for firms
% n_months=10;
% n_firms=2000;
% n_workers=n_firms;

% % Random matrices
% rng(1,"twister");

% r.state=rand(n_months,n_firms); %Inside the meetings, which type of meeting is happening
% r.event=rand(n_months,n_firms); % What kind of meeting is happening
% r.ashock=rand(n_months,n_firms); %Productivity shock
% r.qshock=rand(n_months,n_firms); %Learning shock
% r.type=rand(n_months,n_firms); %Type of match


% % Storage of firm status
% firm_status=zeros(n_months,n_firms); %(1=empty, 2= Firm with manager, 3= Firm with non manager; 4= Firm with team)
% m_ftp=zeros(n_months,n_firms); %Manager type
% n_ftp=zeros(n_months,n_firms); %Non manager type
% a_ftp=zeros(n_months,n_firms); %Productivity type
% m_fwage=zeros(n_months,n_firms); %Manager wage
% n_fwage=zeros(n_months,n_firms); %Non manager wage
% %Flagging the events
% hire_m=zeros(n_months,n_firms); %Hiring of manager ( 1=hire)
% hire_n=zeros(n_months,n_firms); %Hiring of non manager ( 1=hire)
% prom_sm=zeros(n_months,n_firms); %Promotion at SM (-1=demote, 0=stay, 1=promote)
% prom_realloc=zeros(n_months,n_firms); %Promotion at reallocation (-1=demote, 0=stay, 1=promote)
% fire_m=zeros(n_months,n_firms); %Firing of manager at reallocation ( 1=fire)
% fire_n=zeros(n_months,n_firms); %Firing of non manager at reallocation ( 1=fire)


% %Inital conditions
% t=1;
% for i=1:n_firms  
%     firm_status(t,i)=draw_CDF_1d(cdf_firm_dist,r.state(t,i));
%     if (firm_status(t,i)==1) %Empty firm
%         a_ftp(t,i)=draw_CDF_1d(cdf_e_edist,r.ashock(t,i));
%         m_ftp(t,i)=0;
%         n_ftp(t,i)=0;
%         m_fwage(t,i)=0.0;
%         n_fwage(t,i)=0.0;
%     elseif (firm_status(t,i)==2) %Firm with manager
%         [a_ftp(t,i),m_ftp(t,i)]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.event(t,i));
%         n_ftp(t,i)=0;
%         target_wage=(1-bpw)*U(m_ftp(t,i))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
%         m_fwage(t,i)= InterpolateWage(Wm(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
%         n_fwage(t,i)=0.0;
%     elseif (firm_status(t,i)==3) %Firm with non manager
%         [a_ftp(t,i),n_ftp(t,i)]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.event(t,i));
%         m_ftp(t,i)=0;
%         m_fwage(t,i)=0.0;
%         target_wage=(1-bpw)*U(n_ftp(t,i))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
%         n_fwage(t,i)= InterpolateWage(Wn(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%     else
%         %Firm with team
%         [a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.event(t,i));
%         target_wage_m=(1-bpw)*U(m_ftp(t,i))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
%         target_wage_n=(1-bpw)*U(n_ftp(t,i))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
%         m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage_m,wgrid);
%         n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage_n,wgrid);
%     end
% end    

% %Other months
% for i=1:n_firms
%     for t=2:n_months
%         % fprintf('Firm %d, Month %d\n',i,t)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %% Shocks and SM
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % fprintf('SM for status 1')
%         if (firm_status(t-1,i)==1) %Empty firm
%             %Check productivity shock
%             if (r.ashock(t,i)<down_a(a_ftp(t-1,i))) %Firm go down in a
%                 a_ftp(t,i)=a_ftp(t-1,i)-1;
%             elseif (r.ashock(t,i)< down_a(a_ftp(t-1,i))+stay_a(a_ftp(t-1,i)) && r.ashock(t,i)>=down_a(a_ftp(t-1,i))) %Firm stay in a
%                 a_ftp(t,i)=a_ftp(t-1,i);
%             else %Firm go up in a
%                 a_ftp(t,i)=a_ftp(t-1,i)+1;
%             end
%            % Does it match someone?
           
%            if (r.event(t,i)<prob_matching)
%             %Match unemp
%             if (r.state(t,i)<=cdf_emp_dist(1))
%                 zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
%                 % Hire unemp?
%                 if (h_e_u(zu,a_ftp(t,i))>0)
%                     %Where does it go?
%                     if (p_e_m(a_ftp(t,i),zu)>0) %Put as manager
%                         firm_status(t,i)=2;
%                         hire_m(t,i)=1;
%                         m_ftp(t,i)=zu;
%                         n_ftp(t,i)=0;
%                         target_wage=(1-bpw)*U(m_ftp(t,i))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
%                         m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
%                         n_fwage(t,i)=0.0;
%                     else %Put as non manager
%                         firm_status(t,i)=3;
%                         hire_n(t,i)=1;
%                         n_ftp(t,i)=zu;
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         target_wage=(1-bpw)*U(n_ftp(t,i))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
%                         n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                     end
%                 else %Firm does not hire
%                     firm_status(t,i)=1;
%                     m_ftp(t,i)=0;
%                     n_ftp(t,i)=0;
%                     m_fwage(t,i)=0.0;
%                     n_fwage(t,i)=0.0;
%                 end
%             %Match employed
%             elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<=cdf_emp_dist(2))
%                 [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                 % Hire employed manager?
%                 if (h_e_m(am,zm,a_ftp(t,i))>0)
%                     %Where does it go?
%                     if (p_e_m(a_ftp(t,i),zm)>0) %Put as manager
%                         firm_status(t,i)=2;
%                         hire_m(t,i)=1;
%                         m_ftp(t,i)=zm;
%                         n_ftp(t,i)=0;
%                         target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
%                         m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
%                         n_fwage(t,i)=0.0;
%                     else %Put as non manager
%                         firm_status(t,i)=3;
%                         hire_n(t,i)=1;
%                         n_ftp(t,i)=zm;
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
%                         n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                     end
%                 else %Firm does not hire
%                     firm_status(t,i)=1;
%                     m_ftp(t,i)=0;
%                     n_ftp(t,i)=0;
%                     m_fwage(t,i)=0.0;
%                     n_fwage(t,i)=0.0;
%                 end
%             %Match employed non manager
%             elseif (r.state(t,i)>cdf_emp_dist(2) && r.state(t,i)<=cdf_emp_dist(3))
%                 [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                 % Hire employed non manager?
%                 if (h_e_nm(an,zn,a_ftp(t,i))>0)
%                     %Where does it go?
%                     if (p_e_m(a_ftp(i),zn)>0) %Put as manager
%                         firm_status(t,i)=2;
%                         hire_m(t,i)=1;
%                         m_ftp(t,i)=zn;
%                         n_ftp(t,i)=0;
%                         n_fwage(t,i)=0.0;
%                         target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
%                         m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
%                     else %Put as non manager
%                         firm_status(t,i)=3;
%                         hire_n(t,i)=1;
%                         n_ftp(t,i)=zn;
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
%                         n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                     end
%                 else %Firm does not hire
%                     firm_status(t,i)=1;
%                     m_ftp(t,i)=0;
%                     n_ftp(t,i)=0;
%                     m_fwage(t,i)=0.0;
%                     n_fwage(t,i)=0.0;
%                 end
%             %Match with team
%             else 
%                 [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                 % Hire the manager?
%                 if (h_e_t_m(at,zt,nt,a_ftp(t,i))>0)
%                     %Where does it go?
%                     if (p_e_m(a_ftp(t,i),zt)>0) %Put as manager
%                         firm_status(t,i)=2;
%                         hire_m(t,i)=1;
%                         m_ftp(t,i)=zt;
%                         n_ftp(t,i)=0;
%                         target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
%                         m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
%                         n_fwage(t,i)=0.0;
%                     else %Put as non manager
%                         firm_status(t,i)=3;
%                         hire_n(t,i)=1;
%                         n_ftp(t,i)=zt;
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
%                         n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                     end
%                 % Hire the non manager?
%                 elseif (h_e_t_nm(at,zt,nt,a_ftp(t,i))>0)
%                     %Where does it go?
%                     if (p_e_m(a_ftp(t,i),nt)>0) %Put as manager
%                         firm_status(t,i)=2;
%                         hire_m(t,i)=1;
%                         m_ftp(t,i)=nt;
%                         n_ftp(t,i)=0;
%                         target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
%                         m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
%                         n_fwage(t,i)=0.0;
%                     else %Put as non manager
%                         firm_status(t,i)=3;
%                         hire_n(t,i)=1;
%                         n_ftp(t,i)=nt;
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
%                         n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                     end
%                 else %Firm does not hire
%                     firm_status(t,i)=1;
%                     m_ftp(t,i)=0;
%                     n_ftp(t,i)=0;
%                     m_fwage(t,i)=0.0;
%                     n_fwage(t,i)=0.0;
%                 end
%             end
%             else %Firm does not match
%                   firm_status(t,i)=1;
%                   m_ftp(t,i)=0;
%                   n_ftp(t,i)=0;
%                   m_fwage(t,i)=0.0;
%                   n_fwage(t,i)=0.0;
%             end
%         end 
        
%         % fprintf('SM for status 2')
%         if (firm_status(t-1,i)==2) %Firm with manager
%             %Check productivity shock
%             if (r.ashock(t,i)<down_a(a_ftp(t-1,i))) %Firm go down in a
%                 a_ftp(t,i)=a_ftp(t-1,i)-1;
%             elseif (r.ashock(t,i)< down_a(a_ftp(t-1,i))+stay_a(a_ftp(t-1,i)) && r.ashock(t,i)>=down_a(a_ftp(t-1,i))) %Firm stay in a
%                 a_ftp(t,i)=a_ftp(t-1,i);
%             else %Firm go up in a
%                 a_ftp(t,i)=a_ftp(t-1,i)+1;
%             end
%             %Update the states do we can use below
%             m_ftp(t,i)=m_ftp(t-1,i);
%             n_ftp(t,i)=n_ftp(t-1,i);
%             m_fwage(t,i)=m_fwage(t-1,i);
%             n_fwage(t,i)=n_fwage(t-1,i);
%             % Check if worker leaves the firm
%             if (r.event(t,i)<=del+death)
%                 firm_status(t,i)=1;
%                 m_ftp(t,i)=0;
%                 n_ftp(t,i)=0;
%                 m_fwage(t,i)=0.0;
%                 n_fwage(t,i)=0.0;
%             elseif (r.event(t,i)> del+death && r.event(t,i)<= del+death+prob_matching) %Match with someone
%                 %Match unemp
%                 if (r.state(t,i)<=cdf_emp_dist(1))
%                     zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
%                     % Hire unemp?
%                     if (h_m_u(zu,a_ftp(t,i),m_ftp(t-1,i))>0)
%                         %Where does it go?
%                         if (p_m_m_u(a_ftp(t,i),m_ftp(t-1,i),zu)>0) %Put as manager firing current
%                             firm_status(t,i)=2;
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*U(zu)+ bpw*(Vmh(a_ftp(t,i),zu)+U(m_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),zu ),target_wage,wgrid);
%                             n_fwage(t,i)=0.0;
%                             m_ftp(t,i)=zu;
%                             n_ftp(t,i)=0;
%                         elseif (p_m_m_d(a_ftp(t,i),m_ftp(t-1,i),zu)>0) %Put as manager with demotion
%                             firm_status(t,i)=4; %becomes a team
%                             hire_m(t,i)=1;
%                             prom_sm(t,i)=-1; %Demote the manager
%                             target_wage_m=(1-bpw)*U(m_ftp(t,i))+bpw*(Vth(a_ftp(t,i),zu,m_ftp(t,i))-cost_d - Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zu,m_ftp(t,i)),target_wage_m,wgrid);
%                             n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
%                             m_ftp(t,i)=zu;
%                             n_ftp(t,i)=m_ftp(t-1,i); %Demote the manager
%                         else %Put as non manager
%                             firm_status(t,i)=4; %becomes a team
%                             hire_n(t,i)=1;
%                             m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
%                             target_wage=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zu)- Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zu),target_wage,wgrid);
%                             m_ftp(t,i)=m_ftp(t-1,i); %Manager stays
%                             n_ftp(t,i)=zu;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=2;
%                         m_ftp(t,i)=m_ftp(t-1,i);
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                         n_fwage(t,i)=0.0;
%                     end
%                 %Match employed manager
%                 elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<=cdf_emp_dist(2))
%                     [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                     % Hire employed manager?
%                     if (h_m_m(am,zm,a_ftp(t,i),m_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_m_m_u(a_ftp(t,i),m_ftp(t-1,i),zm)>0) %Put as manager firing current
%                             firm_status(t,i)=2;
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vmh(a_ftp(t,i),zm)+U(m_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),zm),target_wage,wgrid);
%                             n_fwage(t,i)=0.0;
%                             m_ftp(t,i)=zm;
%                             n_ftp(t,i)=0;
%                         elseif (p_m_m_d(a_ftp(t,i),m_ftp(t-1,i),zm)>0) %Put as manager with demotion
%                             firm_status(t,i)=4; %becomes a team
%                             hire_m(t,i)=1;
%                             prom_sm(t,i)=-1; %Demote the manager
%                             target_wage_m=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),zm,m_ftp(t,i))-cost_d - Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zm,m_ftp(t,i)),target_wage_m,wgrid);
%                             n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
%                             m_ftp(t,i)=zm;
%                             n_ftp(t,i)=m_ftp(t-1,i); %Demote the manager
%                         else %Put as non manager
%                             firm_status(t,i)=4; %becomes a team
%                             hire_n(t,i)=1;
%                             m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
%                             target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zm)- Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zm),target_wage,wgrid);
%                             m_ftp(t,i)=m_ftp(t-1,i); %Manager stays
%                             n_ftp(t,i)=zm;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=2;
%                         m_ftp(t,i)=m_ftp(t-1,i);
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                         n_fwage(t,i)=0.0;
%                     end
%                 %Match employed non manager
%                 elseif (r.state(t,i)>cdf_emp_dist(2) && r.state(t,i)<=cdf_emp_dist(3))
%                     [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                     % Hire employed non manager?
%                     if (h_m_nm(an,zn,a_ftp(t,i),m_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_m_m_u(a_ftp(t,i),m_ftp(t-1,i),zn)>0) %Put as manager firing current
%                             firm_status(t,i)=2;
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vmh(a_ftp(t,i),zn)+U(m_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),zn),target_wage,wgrid);
%                             n_fwage(t,i)=0.0;
%                             m_ftp(t,i)=zn;
%                             n_ftp(t,i)=0;
%                         elseif (p_m_m_d(a_ftp(t,i),m_ftp(t-1,i),zn)>0) %Put as manager with demotion
%                             firm_status(t,i)=4; %becomes a team
%                             hire_m(t,i)=1;
%                             prom_sm(t,i)=-1; %Demote the manager
%                             target_wage_m=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),zn,m_ftp(t,i))-cost_d - Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zn,m_ftp(t,i)),target_wage_m,wgrid);
%                             n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
%                             m_ftp(t,i)=zn;
%                             n_ftp(t,i)=m_ftp(t-1,i); %Demote the manager
%                         else %Put as non manager
%                             firm_status(t,i)=4; %becomes a team
%                             hire_n(t,i)=1;
%                             m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
%                             target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zn)- Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zn),target_wage,wgrid);
%                             m_ftp(t,i)=m_ftp(t-1,i); %Manager stays
%                             n_ftp(t,i)=zn;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=2;
%                         m_ftp(t,i)=m_ftp(t-1,i);
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                         n_fwage(t,i)=0.0;
%                     end
%                 %Match with team
%                 else 
%                     [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                     % Hire the manager?
%                     if (h_m_t_m(at,zt,nt,a_ftp(t,i),m_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_m_m_u(a_ftp(t,i),m_ftp(t-1,i),zt)>0) %Put as manager firing current
%                             firm_status(t,i)=2;
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vmh(a_ftp(t,i),zt)+U(m_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),zt),target_wage,wgrid);
%                             n_fwage(t,i)=0.0;
%                             m_ftp(t,i)=zt;
%                             n_ftp(t,i)=0;
%                         elseif (p_m_m_d(a_ftp(t,i),m_ftp(t-1,i),zt)>0) %Put as manager with demotion
%                             firm_status(t,i)=4; %becomes a team
%                             hire_m(t,i)=1;
%                             prom_sm(t,i)=-1; %Demote the manager
%                             target_wage_m=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),zt,m_ftp(t,i))-cost_d - Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,m_ftp(t,i)),target_wage_m,wgrid);
%                             n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
%                             m_ftp(t,i)=zt;
%                             n_ftp(t,i)=m_ftp(t-1,i); %Demote the manager
%                         else %Put as non manager
%                             firm_status(t,i)=4; %becomes a team
%                             hire_n(t,i)=1;
%                             m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
%                             target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt)) + bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zt)- Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zt),target_wage,wgrid);
%                             n_ftp(t,i)=zt;
%                         end
%                     % Hire the non manager?
%                     elseif (h_m_t_nm(at,zt,nt,a_ftp(t,i),m_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_m_m_u(a_ftp(t,i),m_ftp(t-1,i),nt)>0) %Put as manager firing current
%                             firm_status(t,i)=2;
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vmh(a_ftp(t,i),nt)+U(m_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),nt),target_wage,wgrid);
%                             n_fwage(t,i)=0.0;
%                             m_ftp(t,i)=nt;
%                             n_ftp(t,i)=0;
%                         elseif (p_m_m_d(a_ftp(t,i),m_ftp(t-1,i),nt)>0) %Put as manager with demotion
%                             firm_status(t,i)=4; %becomes a team
%                             hire_m(t,i)=1;
%                             prom_sm(t,i)=-1; %Demote the manager
%                             target_wage_m=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),nt,m_ftp(t,i))-cost_d - Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),nt,m_ftp(t,i)),target_wage_m,wgrid);
%                             n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
%                             m_ftp(t,i)=nt;
%                             n_ftp(t,i)=m_ftp(t-1,i); %Demote the manager
%                         else %Put as non manager
%                             firm_status(t,i)=4; %becomes a team
%                             hire_n(t,i)=1;
%                             m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
%                             target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),nt)- Vmh(a_ftp(t,i),m_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),nt),target_wage,wgrid);
%                             m_ftp(t,i)=m_ftp(t-1,i); %Manager stays
%                             n_ftp(t,i)=nt;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=2;
%                         m_ftp(t,i)=m_ftp(t-1,i);
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                         n_fwage(t,i)=0.0;
%                     end
%                 end
%             elseif (r.event(t,i)>del+death+prob_matching && r.event(t,i)<=del+death+prob_matching+ lam) % Firm is met by another
%                 if (r.event(t,i)<cdf_firm_dist(1)) %Met By vacant firm
%                     av=draw_CDF_1d(cdf_e_edist,r.type(t,i));
%                     %Does the poacher hire?
%                     if (h_e_m(a_ftp(t,i),m_ftp(t,i),av)>0)
%                         firm_status(t,i)=1;
%                         m_ftp(t,i)=0;
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=0.0;
%                     else %Firm does not hire
%                         firm_status(t,i)=2;
%                         %Check if wage goes up
%                         target_wage=max(Vmh(av,m_ftp(t,i)),Vnh(av,n_ftp(t,i)))-Veh(av); %Marginal of the poaching firm
%                         current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wmh(:,a_ftp(t,i),m_ftp(t,i)));
%                         if (target_wage>current_wage)
%                             m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
%                         else
%                             m_fwage(t,i)=m_fwage(t-1,i);
%                         end
%                         n_ftp(t,i)=0;
%                         n_fwage(t,i)=0.0;
%                     end
%                 elseif (r.event(t,i)>=cdf_firm_dist(1) && r.event(t,i)<cdf_firm_dist(2)) %Met by firm with manager
%                     [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                     %Does the poacher hire?
%                     if (h_m_m(a_ftp(t,i),m_ftp(t,i),am,zm)>0)
%                         firm_status(t,i)=1;
%                         m_ftp(t,i)=0;
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=0.0;
%                     else %Firm does not hire
%                         firm_status(t,i)=2;
%                         %Check if wage goes up
%                         nu=[Vmh(am,m_ftp(t,i))+U(zm),Vth(am,zm,m_ftp(t,i)), Vth(am,m_ftp(t,i),zm)-cost_d];
%                         target_wage_m=max(nu)-Vmh(am,zm); %Marginal of the poaching firm
%                         current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wmh(:,a_ftp(t,i),m_ftp(t,i)));
%                         if (target_wage_m>current_wage)
%                             m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage_m,wgrid);
%                         else
%                             m_fwage(t,i)=m_fwage(t-1,i);
%                         end
%                         n_ftp(t,i)=0;
%                         n_fwage(t,i)=0.0;
%                     end
%                 elseif (r.event(t,i)>=cdf_firm_dist(2) && r.event(t,i)<cdf_firm_dist(3)) %Met by firm with non manager
%                     [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                     %Does the poacher hire?
%                     if (h_nm_m(a_ftp(t,i),m_ftp(t,i),an,zn)>0)
%                         firm_status(t,i)=1;
%                         m_ftp(t,i)=0;
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=0.0;
%                     else %Firm does not hire
%                         firm_status(t,i)=2;
%                         %Check if wage goes up
%                         nu=[Vnh(an,m_ftp(t,i))+U(zn),Vth(an,zn,m_ftp(t,i))-cost_p, Vth(an,m_ftp(t,i),zn)];
%                         target_wage_m=max(nu)-Vnh(an,zn); %Marginal of the poaching firm
%                         current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wmh(:,a_ftp(t,i),m_ftp(t,i)));
%                         if (target_wage_m>current_wage)
%                             m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage_m,wgrid);
%                         else
%                             m_fwage(t,i)=m_fwage(t-1,i);
%                         end
%                         n_ftp(t,i)=0;
%                         n_fwage(t,i)=0.0;
%                     end
%                 else %Met by team
%                     [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                     %Does the poacher hire?
%                     if (h_t_m(a_ftp(t,i),m_ftp(t,i),at,zt,nt)>0)
%                         firm_status(t,i)=1;
%                         m_ftp(t,i)=0;
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=0.0;
%                     else %Firm does not hire
%                         firm_status(t,i)=2;
%                         %Check if wage goes up
%                         nu=[Vth(at,zt,m_ftp(t,i))+U(nt),Vth(at,m_ftp(t,i),zt)-cost_d+U(nt), Vth(at,nt,m_ftp(t,i))+U(zt)-cost_p, Vth(at,m_ftp(t,i),nt)+U(zt)];
%                         target_wage_m=max(nu)-Vth(at,zt,nt) %Marginal of the poaching firm
%                         current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wmh(:,a_ftp(t,i),m_ftp(t,i)));
%                         if (target_wage_m>current_wage)
%                             m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage_m,wgrid);
%                         else
%                             m_fwage(t,i)=m_fwage(t-1,i);
%                         end
%                         n_ftp(t,i)=0;
%                         n_fwage(t,i)=0.0;
%                     end
%                 end
%             else %Nothing happens to the firm at SM
%                   firm_status(t,i)=2;
%                   m_ftp(t,i)=m_ftp(t-1,i);
%                   n_ftp(t,i)=0;
%                   m_fwage(t,i)=m_fwage(t-1,i);
%                   n_fwage(t,i)=0.0;
%             end
%             % %Reallocation post SM
%             % if (d_m(a_ftp(t,i),m_ftp(t,i))>0) %Fire the manager
%             %     firm_status(t,i)=1;
%             %     m_ftp(t,i)=0;
%             %     n_ftp(t,i)=0;
%             %     m_fwage(t,i)=0.0;
%             %     n_fwage(t,i)=0.0;
%             % elseif (r_m(a_ftp(t,i),m_ftp(t,i))>0) %reallocate the manager
%             %     firm_status(t,i)=3;
%             %     n_ftp(t,i)=m_ftp(t,i);
%             %     n_fwage(t,i)=m_fwage(t,i); %Wage does not change upon reallocation
%             %     m_ftp(t,i)=0;
%             %     m_fwage(t,i)=0.0;
%             % else %Manager stays
%             %     firm_status(t,i)=2;
%             %     %Does the wage goes down?
%             %     current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wmh(:,a_ftp(t,i),m_ftp(t,i)));
%             %     target_wage=Vmh(a_ftp(t,i),m_ftp(t,i))-Veh(a_ftp(t,i));
%             %     if (target_wage<current_wage)
%             %         m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
%             %     else
%             %         m_fwage(t,i)=m_fwage(t-1,i);
%             %     end
%             %     n_ftp(t,i)=0;
%             %     n_fwage(t,i)=0.0;
%             % end
%         end
        
%         % fprintf('SM for status 3')
%         if (firm_status(t-1,i)==3) %Firm with non manager
%             %Check productivity shock
%             if (r.ashock(t,i)<down_a(a_ftp(t-1,i))) %Firm go down in a
%                 a_ftp(t,i)=a_ftp(t-1,i)-1;
%             elseif (r.ashock(t,i)< down_a(a_ftp(t-1,i))+stay_a(a_ftp(t-1,i)) && r.ashock(t,i)>=down_a(a_ftp(t-1,i))) %Firm stay in a
%                 a_ftp(t,i)=a_ftp(t-1,i);
%             else %Firm go up in a
%                 a_ftp(t,i)=a_ftp(t-1,i)+1;
%             end
%             %Check Learning shock on the non manager
%             if (r.qshock(t,i)<down_q(a_ftp(t-1,i), n_ftp(t-1,i))) %Firm go down in q
%                 n_ftp(t,i)=n_ftp(t-1,i)-1;    
%             elseif (r.qshock(t,i)< down_q(a_ftp(t-1,i), n_ftp(t-1,i))+stay_q(a_ftp(t-1,i), n_ftp(t-1,i)) && r.qshock(t,i)>=down_q(a_ftp(t-1,i), n_ftp(t-1,i))) %Firm stay in q
%                 n_ftp(t,i)=n_ftp(t-1,i);
%             else %Firm go up in q
%                 n_ftp(t,i)=n_ftp(t-1,i)+1;
%             end
%             %Update the states do we can use below
%             m_ftp(t,i)=m_ftp(t-1,i);
%             m_fwage(t,i)=m_fwage(t-1,i);
%             n_fwage(t,i)=n_fwage(t-1,i);
%             % Check if worker leaves the firm
%             if (r.event(t,i)<=del+death)
%                 firm_status(t,i)=1;
%                 m_ftp(t,i)=0;
%                 n_ftp(t,i)=0;
%                 m_fwage(t,i)=0.0;
%                 n_fwage(t,i)=0.0;
%             elseif (r.event(t,i)> del+death && r.event(t,i)<= del+death+prob_matching) %Match with someone
%                 %Match unemp
%                 if (r.state(t,i)<=cdf_emp_dist(1)) zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
%                     % Hire unemp?
%                     if (h_nm_u(zu,a_ftp(t,i),n_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_n_n_u(a_ftp(t,i),n_ftp(t,i),zu)>0) %Put as non manager firing current
%                             firm_status(t,i)=3;
%                             hire_n(t,i)=1;
%                             target_wage=(1-bpw)*U(zu)+bpw*(Vnh(a_ftp(t,i),zu)+U(n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),zu),target_wage,wgrid);
%                             m_fwage(t,i)=0.0;
%                             m_ftp(t,i)=0;
%                             n_ftp(t,i)=zu;
%                         elseif (p_n_n_p(a_ftp(t,i),n_ftp(t,i),zu)>0) %Put as non manager with promotion
%                             firm_status(t,i)=4; %becomes a team
%                             hire_n(t,i)=1;
%                             prom_sm(t,i)=1; %Promote the non manager 
%                             target_wage_n=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zu)-cost_p - Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zu),target_wage_n,wgrid);
%                             m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
%                             m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
%                             n_ftp(t,i)=zu;
%                         else %Put as manager
%                             firm_status(t,i)=4; %becomes a team
%                             hire_m(t,i)=1;
%                             n_fwage(t,i)=n_fwage(t-1,i); %Wage of the non manager do not change now
%                             target_wage=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),zu,n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zu,n_ftp(t,i)),target_wage,wgrid);
%                             m_ftp(t,i)=zu;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=3;
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     end
%                 %Match employed manager
%                 elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<=cdf_emp_dist(2))
%                     %Match employed manager
%                     [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                     % Hire employed manager?
%                     if (h_nm_m(am,zm,a_ftp(t,i),n_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_n_n_u(a_ftp(t,i),n_ftp(t,i),zm)>0) %Put as non manager firing current
%                             firm_status(t,i)=3;
%                             hire_n(t,i)=1;
%                             target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vnh(a_ftp(t,i),zm)+U(n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),zm),target_wage,wgrid);
%                             m_fwage(t,i)=0.0;
%                             m_ftp(t,i)=0;
%                             n_ftp(t,i)=zm;
%                         elseif (p_n_n_p(a_ftp(t,i),n_ftp(t,i),zm)>0) %Put as non manager with promotion
%                             firm_status(t,i)=4; %becomes a team
%                             hire_n(t,i)=1;
%                             prom_sm(t,i)=1; %Promote the non manager
%                             target_wage_n=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zm)-cost_p - Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zm),target_wage_n,wgrid);
%                             m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
%                             m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
%                             n_ftp(t,i)=zm;
%                         else %Put as manager
%                             firm_status(t,i)=4; %becomes a team
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),zm,n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zm,n_ftp(t,i)),target_wage,wgrid);
%                             n_fwage(t,i)=n_fwage(t-1,i); %Wage of the non manager do not change now
%                             m_ftp(t,i)=zm;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=3;
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     end
%                 %Match employed non manager
%                 elseif (r.state(t,i)>cdf_emp_dist(2) && r.state(t,i)<=cdf_emp_dist(3))
%                     [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                     % Hire employed non manager?
%                     if (h_nm_nm(an,zn,a_ftp(t,i),n_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_n_n_u(a_ftp(t,i),n_ftp(t,i),zn)>0) %Put as non manager firing current
%                             firm_status(t,i)=3;
%                             hire_n(t,i)=1;
%                             target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vnh(a_ftp(t,i),zn)+U(n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),zn),target_wage,wgrid);
%                             m_fwage(t,i)=0.0;
%                             m_ftp(t,i)=0;
%                             n_ftp(t,i)=zn;
%                         elseif (p_n_n_p(a_ftp(t,i),n_ftp(t,i),zn)>0) %Put as non manager with promotion
%                             firm_status(t,i)=4; %becomes a team
%                             hire_n(t,i)=1;
%                             prom_sm(t,i)=1; %Promote the non manager
%                             target_wage_n=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zn)-cost_p - Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zn),target_wage_n,wgrid);
%                             m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
%                             m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
%                             n_ftp(t,i)=zn;
%                         else %Put as manager
%                             firm_status(t,i)=4; %becomes a team
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),zn,n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zn,n_ftp(t,i)),target_wage,wgrid);
%                             n_fwage(t,i)=n_fwage(t-1,i); %Wage of the non
%                             m_ftp(t,i)=zn;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=3;
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     end
%                 %Match with team
%                 else 
%                     [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                     % Hire the manager?
%                     if (h_nm_t_m(at,zt,nt,a_ftp(t,i),n_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_n_n_u(a_ftp(t,i),n_ftp(t,i),zt)>0) %Put as non manager firing current
%                             firm_status(t,i)=3;
%                             hire_n(t,i)=1;
%                             target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vnh(a_ftp(t,i),zt)+U(n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),zt),target_wage,wgrid);
%                             m_fwage(t,i)=0.0;
%                             m_ftp(t,i)=0;
%                             n_ftp(t,i)=zt;
%                         elseif (p_n_n_p(a_ftp(t,i),n_ftp(t,i),zt)>0) %Put as non manager with promotion
%                             firm_status(t,i)=4; %becomes a team
%                             hire_n(t,i)=1;
%                             prom_sm(t,i)=1; %Promote the non manager
%                             target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zt)-cost_p - Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zt),target_wage_n,wgrid);
%                             m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
%                             m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
%                             n_ftp(t,i)=zt;
%                         else %Put as manager
%                             firm_status(t,i)=4; %becomes a team
%                             hire_m(t,i)=1;
%                             n_fwage(t,i)=n_fwage(t-1,i); %Wage of the non manager do not change now
%                             target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),zt,n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,n_ftp(t,i)),target_wage,wgrid);
%                             m_ftp(t,i)=zt;
%                         end
%                     % Hire the non manager?
%                     elseif (h_nm_t_nm(at,zt,nt,a_ftp(t,i),n_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_n_n_u(a_ftp(t,i),n_ftp(t,i),nt)>0) %Put as non manager firing current
%                             firm_status(t,i)=3;
%                             hire_n(t,i)=1;
%                             target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vnh(a_ftp(t,i),nt)+U(n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),nt),target_wage,wgrid);
%                             m_fwage(t,i)=0.0;
%                             m_ftp(t,i)=0;
%                             n_ftp(t,i)=nt;
%                         elseif (p_n_n_p(a_ftp(t,i),n_ftp(t,i),nt)>0) %Put as non manager with promotion
%                             firm_status(t,i)=4; %becomes a team
%                             hire_n(t,i)=1;
%                             prom_sm(t,i)=1; %Promote the non manager
%                             target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),nt)-cost_p - Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),nt),target_wage_n,wgrid);
%                             m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
%                             m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
%                             n_ftp(t,i)=nt;
%                         else %Put as manager
%                             firm_status(t,i)=4; %becomes a team
%                             hire_m(t,i)=1;
%                             n_fwage(t,i)=n_fwage(t-1,i); %Wage of the non manager do not change now
%                             target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),nt,n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),nt,n_ftp(t,i)),target_wage,wgrid);
%                             m_ftp(t,i)=nt;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=3;
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     end
%                 end
%             elseif (r.event(t,i)>del+death+prob_matching && r.event(t,i)<=del+death+prob_matching+ lam) % Firm is met by another
%                 if (r.event(t,i)<cdf_firm_dist(1)) %Met By vacant firm
%                     av=draw_CDF_1d(cdf_e_edist,r.type(t,i));
%                     %Does the poacher hire?
%                     if (h_e_nm(a_ftp(t,i),n_ftp(t,i),av)>0)
%                         firm_status(t,i)=1;
%                         m_ftp(t,i)=0;
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=0.0;
%                     else %Firm does not hire
%                         firm_status(t,i)=3;
%                         %Check if wage goes up
%                         target_wage=max(Vnh(av,n_ftp(t,i)),Vmh(av,n_ftp(t,i)))-Veh(av); %Marginal of the poaching firm
%                         current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wnh(:,a_ftp(t,i),n_ftp(t,i)));
%                         if (target_wage>current_wage)
%                             n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                         else
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                         end
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                     end
%                 elseif (r.event(t,i)>=cdf_firm_dist(1) && r.event(t,i)<cdf_firm_dist(2)) %Met by firm with manager
%                     [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                     %Does the poacher hire?
%                     if (h_m_nm(a_ftp(t,i),n_ftp(t,i),am,zm)>0)
%                         firm_status(t,i)=1;
%                         m_ftp(t,i)=0;
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=0.0;
%                     else %Firm does not hire
%                         firm_status(t,i)=3;
%                         %Check if wage goes up
%                         nu=[Vmh(am,n_ftp(t,i))+U(zm),Vth(am,zm,n_ftp(t,i)), Vth(am,n_ftp(t,i),zm)-cost_d];
%                         target_wage_m=max(nu)-Vmh(am,zm); %Marginal of the poaching firm
%                         current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wnh(:,a_ftp(t,i),n_ftp(t,i)));
%                         if (target_wage_m>current_wage)
%                             n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage_m,wgrid);
%                         else
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                         end
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                     end
%                 elseif (r.event(t,i)>=cdf_firm_dist(2) && r.event(t,i)<cdf_firm_dist(3)) %Met by firm with non manager
%                     [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                     %Does the poacher hire?
%                     if (h_nm_nm(a_ftp(t,i),n_ftp(t,i),an,zn)>0)
%                         firm_status(t,i)=1;
%                         m_ftp(t,i)=0;
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=0.0;
%                     else %Firm does not hire
%                         firm_status(t,i)=3;
%                         %Check if wage goes up
%                         nu=[Vnh(an,n_ftp(t,i))+U(zn),Vth(an,zn,n_ftp(t,i))-cost_p, Vth(an,n_ftp(t,i),zn)];
%                         target_wage_m=max(nu)-Vnh(an,zn); %Marginal of the poaching firm
%                         current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wnh(:,a_ftp(t,i),n_ftp(t,i)));
%                         if (target_wage_m>current_wage)
%                             n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage_m,wgrid);
%                         else
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                         end
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                     end
%                 else %Met by team
%                     [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                     %Does the poacher hire?
%                     if (h_t_nm(a_ftp(t,i),n_ftp(t,i),at,zt,nt)>0)
%                         firm_status(t,i)=1;
%                         m_ftp(t,i)=0;
%                         n_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=0.0;
%                     else %Firm does not hire
%                         firm_status(t,i)=3;
%                         %Check if wage goes up
%                         nu=[Vth(at,zt,n_ftp(t,i))+U(nt),Vth(at,n_ftp(t,i),zt)-cost_d+U(nt), Vth(at,nt,n_ftp(t,i))+U(zt)-cost_p, Vth(at,n_ftp(t,i),nt)+U(zt)];
%                         target_wage_m=max(nu)-Vth(at,zt,nt) %Marginal of the poaching firm
%                         current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wnh(:,a_ftp(t,i),n_ftp(t,i)));
%                         if (target_wage_m>current_wage)
%                             n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage_m,wgrid);
%                         else
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                         end
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                     end
%                 end
%             else %Nothing happens to the firm at SM
%                   firm_status(t,i)=3;
%                   m_ftp(t,i)=0;
%                   n_ftp(t,i)=n_ftp(t-1,i);
%                   m_fwage(t,i)=0.0;
%                   n_fwage(t,i)=n_fwage(t-1,i);
%             end
%             % %Reallocation post SM
%             % if (d_n(a_ftp(t,i),n_ftp(t,i))>0) %Fire the non manager
%             %     firm_status(t,i)=1;
%             %     m_ftp(t,i)=0;
%             %     n_ftp(t,i)=0;
%             %     m_fwage(t,i)=0.0;
%             %     n_fwage(t,i)=0.0;
%             % elseif (r_n(a_ftp(t,i),n_ftp(t,i))>0) %reallocate the non manager
%             %     firm_status(t,i)=2; %Firm with manager
%             %     m_ftp(t,i)=n_ftp(t,i);
%             %     m_fwage(t,i)=n_fwage(t,i); %Wage does not change upon reallocation
%             %     n_ftp(t,i)=0;
%             %     n_fwage(t,i)=0.0;
%             % else %Non manager stays
%             %     firm_status(t,i)=3;
%             %     %Does the wage goes down?
%             %     current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wnh(:,a_ftp(t,i),n_ftp(t,i)));
%             %     target_wage=Vnh(a_ftp(t,i),n_ftp(t,i))-Veh(a_ftp(t,i));
%             %     if (target_wage<current_wage)
%             %         n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%             %     else
%             %         n_fwage(t,i)=n_fwage(t-1,i);
%             %     end
%             %     m_ftp(t,i)=0;
%             %     m_fwage(t,i)=0.0;
%             % end
%         end
                        
%         % fprintf('SM for status 4')
%         if (firm_status(t-1,i)==4) %Firm with team
%             %Check productivity shock
%             if (r.ashock(t,i)<down_a(a_ftp(t-1,i))) %Firm go down in a
%                 a_ftp(t,i)=a_ftp(t-1,i)-1;
%             elseif (r.ashock(t,i)< down_a(a_ftp(t-1,i))+stay_a(a_ftp(t-1,i)) && r.ashock(t,i)>=down_a(a_ftp(t-1,i))) %Firm stay in a
%                 a_ftp(t,i)=a_ftp(t-1,i);
%             else %Firm go up in a
%                 a_ftp(t,i)=a_ftp(t-1,i)+1;
%             end
%             %Check Learning shock on the non manager
%             if (r.qshock(t,i)<down_q(a_ftp(t-1,i), n_ftp(t-1,i))) %Firm go down in q
%                 n_ftp(t,i)=n_ftp(t-1,i)-1;
%             elseif (r.qshock(t,i)< down_q(a_ftp(t-1,i), n_ftp(t-1,i)) + stay_q(a_ftp(t-1,i),n_ftp(t-1,i)) && r.qshock(t,i)>=down_q(a_ftp(t-1,i), n_ftp(t-1,i))) %Firm stay in q
%                 n_ftp(t,i)=n_ftp(t-1,i);
%             else %Firm go up in q
%                 n_ftp(t,i)=n_ftp(t-1,i)+1;
%             end
%             %Update the states do we can use below
%             m_ftp(t,i)=m_ftp(t-1,i);
%             m_fwage(t,i)=m_fwage(t-1,i);
%             n_fwage(t,i)=n_fwage(t-1,i);
%             % Check if manager leaves the firm
%             if (r.event(t,i)<=del+death)
%                 firm_status(t,i)=3;
%                 m_ftp(t,i)=0;
%                 n_ftp(t,i)=n_ftp(t-1,i);
%                 m_fwage(t,i)=0.0;
%                 n_fwage(t,i)=n_fwage(t-1,i);
%             %check if non manager leaves the firm
%             elseif (r.event(t,i)> del+death && r.event(t,i)<= del+death+del+death) 
%                 firm_status(t,i)=2;
%                 m_ftp(t,i)=m_ftp(t-1,i);
%                 n_ftp(t,i)=0;
%                 m_fwage(t,i)=m_fwage(t-1,i);
%                 n_fwage(t,i)=0.0;
%             elseif (r.event(t,i)> del+death+del+death && r.event(t,i)<= del+death+del+death+prob_matching) %Match with someone
%                 %Match unemp
%                 if (r.state(t,i)<=cdf_emp_dist(1)) 
%                     zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
%                     % Hire unemp?
%                     if (h_t_u(zu,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_t_m_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zu)>0) %Put as manager firing current
%                             firm_status(t,i)=4;
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),zu,n_ftp(t,i))+U(m_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zu,n_ftp(t,i)),target_wage,wgrid);
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                             m_ftp(t,i)=zu;
%                         elseif (p_t_m_d(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zu)>0) %Put as manager with demotion
%                             firm_status(t,i)=4; %still a team
%                             hire_m(t,i)=1;
%                             prom_sm(t,i)=-1; %Demote the manager
%                             target_wage_m=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),zu,m_ftp(t,i))+U(n_ftp(t,i))-cost_d- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zu,m_ftp(t,i)),target_wage_m,wgrid);
%                             n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
%                             n_ftp(t,i)=m_ftp(t,i); %Demote the manager
%                             m_ftp(t,i)=zu;
%                         elseif (p_t_n_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zu)>0) %Put as non manager firing current
%                             firm_status(t,i)=4; %still a team
%                             hire_n(t,i)=1;
%                             target_wage_n=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zu)+U(n_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zu),target_wage_n,wgrid);
%                             m_fwage(t,i)=m_fwage(t-1,i); %Wage of the non manager do not change now
%                             n_ftp(t,i)=zu;
%                         else %Put as non manager with promotion
%                             firm_status(t,i)=4; %still a team
%                             hire_n(t,i)=1;
%                             prom_sm(t,i)=1; %Promote the non manager
%                             target_wage_n=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zu)-cost_p+ U(m_ftp(t,i)) -Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zu),target_wage_n,wgrid);
%                             m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
%                             m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
%                             n_ftp(t,i)=zu;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=4;
%                         m_ftp(t,i)=m_ftp(t-1,i);
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     end
%                 %Match employed manager
%                 elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<=cdf_emp_dist(2))
%                     %Match employed manager
%                     [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                     % Hire employed manager?
%                     if (h_t_m(am,zm,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_t_m_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zm)>0) %Put as manager firing current
%                             firm_status(t,i)=4;
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),zm,n_ftp(t,i))+U(m_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zm,n_ftp(t,i)),target_wage,wgrid);
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                             m_ftp(t,i)=zm;
%                         elseif (p_t_m_d(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zm)>0) %Put as manager with demotion
%                             firm_status(t,i)=4; %still a team
%                             hire_m(t,i)=1;
%                             prom_sm(t,i)=-1; %Demote the manager
%                             target_wage_m=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),zm,m_ftp(t,i))+U(n_ftp(t,i))-cost_d- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zm,m_ftp(t,i)),target_wage_m,wgrid);
%                             n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
%                             n_ftp(t,i)=m_ftp(t,i); %Demote the manager
%                             m_ftp(t,i)=zm;
%                         elseif (p_t_n_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zm)>0) %Put as non manager firing current
%                             firm_status(t,i)=4; %still a team
%                             hire_n(t,i)=1;
%                             target_wage_n=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zm)+U(n_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zm),target_wage_n,wgrid);
%                             m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
%                             n_ftp(t,i)=zm;
%                         else %Put as non manager with promotion
%                             firm_status(t,i)=4; %still a team
%                             hire_n(t,i)=1;
%                             prom_sm(t,i)=1; %Promote the non manager
%                             target_wage_n=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zm)-cost_p + U(m_ftp(t,i))-Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zm),target_wage_n,wgrid);
%                             m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
%                             m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
%                             n_ftp(t,i)=zm;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=4;
%                         m_ftp(t,i)=m_ftp(t-1,i);
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     end
%                 %Match employed non manager
%                 elseif (r.state(t,i)>cdf_emp_dist(2) && r.state(t,i)<=cdf_emp_dist(3))
%                     [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                     % Hire employed non manager?
%                     if (h_t_nm(an,zn,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_t_m_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zn)>0) %Put as manager firing current
%                             firm_status(t,i)=4;
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),zn,n_ftp(t,i))+U(m_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zn,n_ftp(t,i)),target_wage,wgrid);
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                             m_ftp(t,i)=zn;
%                         elseif (p_t_m_d(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zn)>0) %Put as manager with demotion
%                             firm_status(t,i)=4; %still a team
%                             hire_m(t,i)=1;
%                             prom_sm(t,i)=-1; %Demote the manager
%                             target_wage_m=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),zn,m_ftp(t,i))+U(n_ftp(t,i))-cost_d- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zn,m_ftp(t,i)),target_wage_m,wgrid);
%                             n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
%                             n_ftp(t,i)=m_ftp(t,i); %Demote the manager
%                             m_ftp(t,i)=zn;
%                         elseif (p_t_n_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zn)>0) %Put as non manager firing current
%                             firm_status(t,i)=4; %still a team
%                             hire_n(t,i)=1;
%                             target_wage_n=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zn)+U(n_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zn),target_wage_n,wgrid);
%                             m_fwage(t,i)=m_fwage(t-1,i); %Wage of the non manager do not change now
%                             n_ftp(t,i)=zn;
%                         else %Put as non manager with promotion
%                             firm_status(t,i)=4; %still a team
%                             hire_n(t,i)=1;
%                             prom_sm(t,i)=1; %Promote the non manager
%                             target_wage_n=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zn)-cost_p + U(m_ftp(t,i))-Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zn),target_wage_n,wgrid);
%                             m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
%                             m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
%                             n_ftp(t,i)=zn;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=4;
%                         m_ftp(t,i)=m_ftp(t-1,i);
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     end
%                 %Match with team
%                 else 
%                     [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                     % Hire the manager?
%                     if (h_t_t_m(at,zt,nt,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_t_m_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zt)>0) %Put as manager firing current
%                             firm_status(t,i)=4;
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),zt,n_ftp(t,i))+U(m_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,n_ftp(t,i)),target_wage,wgrid);
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                             m_ftp(t,i)=zt;
%                         elseif (p_t_m_d(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zt)>0) %Put as manager with demotion
%                             firm_status(t,i)=4; %still a team
%                             hire_m(t,i)=1;
%                             prom_sm(t,i)=-1; %Demote the manager
%                             target_wage_m=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),zt,m_ftp(t,i))+U(n_ftp(t,i))-cost_d- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,m_ftp(t,i)),target_wage_m,wgrid);
%                             n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
%                             n_ftp(t,i)=m_ftp(t,i); %Demote the manager
%                             m_ftp(t,i)=zt;
%                         elseif (p_t_n_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zt)>0) %Put as non manager firing current
%                             firm_status(t,i)=4; %still a team
%                             hire_n(t,i)=1;
%                             target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zt)+U(n_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zt),target_wage_n,wgrid);
%                             m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
%                             n_ftp(t,i)=zt;
%                         else %Put as non manager with promotion
%                             firm_status(t,i)=4; %still a team
%                             hire_n(t,i)=1;
%                             prom_sm(t,i)=1; %Promote the non manager
%                             target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zt)-cost_p + U(m_ftp(t,i))-Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zt),target_wage_n,wgrid);
%                             m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
%                             m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
%                             n_ftp(t,i)=zt;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=4;
%                         m_ftp(t,i)=m_ftp(t-1,i);
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     end
        
%                     % Hire the non manager?
%                     if (h_t_t_nm(at,zt,nt,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0)
%                         %Where does it go?
%                         if (p_t_m_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),nt)>0) %Put as manager firing current
%                             firm_status(t,i)=4;
%                             hire_m(t,i)=1;
%                             target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),zt,n_ftp(t,i))+U(m_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,n_ftp(t,i)),target_wage,wgrid);
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                             m_ftp(t,i)=zt;
%                         elseif (p_t_m_d(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),nt)>0) %Put as manager with demotion
%                             firm_status(t,i)=4; %still a team
%                             hire_m(t,i)=1;
%                             prom_sm(t,i)=-1; %Demote the manager
%                             target_wage_m=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),zt,m_ftp(t,i))+U(n_ftp(t,i))-cost_d- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,m_ftp(t,i)),target_wage_m,wgrid);
%                             n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
%                             n_ftp(t,i)=m_ftp(t,i); %Demote the manager
%                             m_ftp(t,i)=zt;
%                         elseif (p_t_n_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),nt)>0) %Put as non manager firing current
%                             firm_status(t,i)=4; %still a team
%                             hire_n(t,i)=1;
%                             target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zt)+U(n_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zt),target_wage_n,wgrid);
%                             m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
%                             n_ftp(t,i)=zt;
%                         else %Put as non manager with promotion
%                             firm_status(t,i)=4; %still a team
%                             hire_n(t,i)=1;
%                             prom_sm(t,i)=1; %Promote the non manager
%                             target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zt)-cost_p + U(m_ftp(t,i))-Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zt),target_wage_n,wgrid);
%                             m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
%                             m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
%                             n_ftp(t,i)=zt;
%                         end
%                     else %Firm does not hire
%                         firm_status(t,i)=4;
%                         m_ftp(t,i)=m_ftp(t-1,i);
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     end
%                 end
%             elseif (r.event(t,i)>del+death+prob_matching && r.event(t,i)<=del+death+prob_matching+ lam) % Firm is met by another
%                 if (r.event(t,i)<cdf_firm_dist(1)) %Met By vacant firm
%                     av=draw_CDF_1d(cdf_e_edist,r.type(t,i));
%                     %Does the poacher hire the manager?
%                     if (h_e_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),av)>0)
%                         firm_status(t,i)=3; %becomes a firm with non manager
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     %Does the poacher hire the non manager?    
%                     elseif (h_e_t_nm(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),av)>0)
%                         firm_status(t,i)=2; %becomes a firm with manager
%                         n_ftp(t,i)=0;
%                         n_fwage(t,i)=0.0;
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                     else %Firm does not hire
%                         firm_status(t,i)=4;
%                         %Check if wage goes up
%                         %Who had the biggest gain from trade? This will determine which, if any, wage goes up
%                         gain_m=max(Vmh(av,m_ftp(t,i)),Vnh(av,m_ftp(t,i)))-Veh(av) - Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vnh(a_ftp(t,i),n_ftp(t,i));
%                         gain_n=max(Vmh(av,n_ftp(t,i)),Vnh(av,n_ftp(t,i)))-Veh(av) - Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vmh(a_ftp(t,i),m_ftp(t,i));
%                         if (gain_m>=gain_n) %Manager might get a raise
%                             target_wage=max(Vmh(av,m_ftp(t,i)),Vnh(av,m_ftp(t,i)))-Veh(av); %Marginal of the poaching firm
%                             current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             if (target_wage>current_wage)
%                                 m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                             else
%                                 m_fwage(t,i)=m_fwage(t-1,i);
%                             end
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                         else %Non manager might get a raise
%                             target_wage=max(Vmh(av,n_ftp(t,i)),Vnh(av,n_ftp(t,i)))-Veh(av); %Marginal of the poaching firm
%                             current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             if (target_wage>current_wage)
%                                 n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                             else
%                                 n_fwage(t,i)=n_fwage(t-1,i);
%                             end
%                             m_fwage(t,i)=m_fwage(t-1,i);
%                         end
%                     end
%                 elseif (r.event(t,i)>=cdf_firm_dist(1) && r.event(t,i)<cdf_firm_dist(2)) %Met by firm with manager
%                     [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                     %Does the poacher hire the manager?
%                     if (h_m_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),am,zm)>0)
%                         firm_status(t,i)=3; %becomes a firm with non manager
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     %Does the poacher hire the non manager?    
%                     elseif (h_m_t_nm(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),am,zm)>0)
%                         firm_status(t,i)=2; %becomes a firm with manager
%                         n_ftp(t,i)=0;
%                         n_fwage(t,i)=0.0;
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                     else %Firm does not hire
%                         firm_status(t,i)=4;
%                         %Check if wage goes up
%                         %Who had the biggest gain from trade? This will determine which, if any, wage goes up
%                         nu_m=[Vmh(am,m_ftp(t,i))+U(zm),Vth(am,zm,m_ftp(t,i)), Vth(am,m_ftp(t,i),zm)-cost_d];
%                         nu_n=[Vmh(am,n_ftp(t,i))+U(zm),Vth(am,zm,n_ftp(t,i)), Vth(am,n_ftp(t,i),zm)];
%                         gain_m=max(nu_m)-Vmh(am,zm)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vnh(a_ftp(t,i),n_ftp(t,i));
%                         gain_n=max(nu_n)-Vnh(am,zm)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vmh(a_ftp(t,i),m_ftp(t,i));
%                         if (gain_m>=gain_n) %Manager might get a raise
%                             target_wage=max(nu_m)-Vmh(am,zm); %Marginal of the poaching firm
%                             current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             if (target_wage>current_wage)
%                                 m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                             else
%                                 m_fwage(t,i)=m_fwage(t-1,i);
%                             end
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                         else %Non manager might get a raise
%                             target_wage=max(nu_n)-Vnh(am,zm); %Marginal of the poaching firm
%                             current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             if (target_wage>current_wage)
%                                 n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                             else
%                                 n_fwage(t,i)=n_fwage(t-1,i);
%                             end
%                             m_fwage(t,i)=m_fwage(t-1,i);
%                         end
%                     end
%                 elseif (r.event(t,i)>=cdf_firm_dist(2) && r.event(t,i)<cdf_firm_dist(3)) %Met by firm with non manager
%                     [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                     %Does the poacher hire the manager?
%                     if (h_nm_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),an,zn)>0)
%                         firm_status(t,i)=3; %becomes a firm with non manager
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     %Does the poacher hire the non manager?    
%                     elseif (h_nm_t_nm(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),an,zn)>0)
%                         firm_status(t,i)=2; %becomes a firm with manager
%                         n_ftp(t,i)=0;
%                         n_fwage(t,i)=0.0;
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                     else %Firm does not hire
%                         firm_status(t,i)=4;
%                         %Check if wage goes up
%                         %Who had the biggest gain from trade? This will determine which, if any, wage goes up
%                         nu_m=[Vnh(an,m_ftp(t,i))+U(zn),Vth(an,zn,m_ftp(t,i)), Vth(an,m_ftp(t,i),zn)-cost_d];
%                         nu_n=[Vnh(an,n_ftp(t,i))+U(zn),Vth(an,zn,n_ftp(t,i))-cost_p, Vth(an,n_ftp(t,i),zn)];
%                         gain_m=max(nu_m)-Vnh(an,zn)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vnh(a_ftp(t,i),n_ftp(t,i));
%                         gain_n=max(nu_n)-Vnh(an,zn)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vmh(a_ftp(t,i),m_ftp(t,i));
%                         if (gain_m>=gain_n) %Manager might get a raise
%                             target_wage=max(nu_m)-Vnh(an,zn); %Marginal of the poaching firm
%                             current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             if (target_wage>current_wage)
%                                 m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                             else
%                                 m_fwage(t,i)=m_fwage(t-1,i);
%                             end
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                         else %Non manager might get a raise
%                             target_wage=max(nu_n)-Vnh(an,zn); %Marginal of the poaching firm
%                             current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             if (target_wage>current_wage)
%                                 n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                             else
%                                 n_fwage(t,i)=n_fwage(t-1,i);
%                             end
%                             m_fwage(t,i)=m_fwage(t-1,i);
%                         end
%                     end
%                 else %Met by firm with team
%                     [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                     %Does the poacher hire the manager?
%                     if (h_t_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),at,zt,nt)>0)
%                         firm_status(t,i)=3; %becomes a firm with non manager
%                         m_ftp(t,i)=0;
%                         m_fwage(t,i)=0.0;
%                         n_fwage(t,i)=n_fwage(t-1,i);
%                     %Does the poacher hire the non manager?
%                     elseif (h_t_t_nm(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),at,zt,nt)>0)
%                         firm_status(t,i)=2; %becomes a firm with manager
%                         n_ftp(t,i)=0;
%                         n_fwage(t,i)=0.0;
%                         m_fwage(t,i)=m_fwage(t-1,i);
%                     else %Firm does not hire
%                         firm_status(t,i)=4;
%                         %Check if wage goes up
%                         %Who had the biggest gain from trade? This will determine which, if any, wage goes up
%                         nu_m=[Vth(at,zt,m_ftp(t,i))+U(nt),Vth(at,m_ftp(t,i),zt)-cost_d+U(nt), Vth(at,m_ftp(t,i),nt)+U(zt), Vth(at,nt,m_ftp(t,i))+U(zt)-cost_p];
%                         nu_n=[Vth(at,zt,n_ftp(t,i))+U(nt),Vth(at,n_ftp(t,i),zt)-cost_d+U(nt), Vth(at,n_ftp(t,i),nt)+U(zt), Vth(at,nt,n_ftp(t,i))+U(zt)-cost_p];
%                         gain_m=max(nu_m)-Vth(at,zt,nt)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vnh(a_ftp(t,i),n_ftp(t,i));
%                         gain_n=max(nu_n)-Vth(at,zt,nt)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vmh(a_ftp(t,i),m_ftp(t,i));
%                         if (gain_m>=gain_n) %Manager might get a raise
%                             target_wage=max(nu_m)-Vth(at,zt,nt); %Marginal of the poaching firm
%                             current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             if (target_wage>current_wage)
%                                 m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                             else
%                                 m_fwage(t,i)=m_fwage(t-1,i);
%                             end
%                             n_fwage(t,i)=n_fwage(t-1,i);
%                         else %Non manager might get a raise
%                             target_wage=max(nu_n)-Vth(at,zt,nt); %Marginal of the poaching firm
%                             current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                             if (target_wage>current_wage)
%                                 n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                             else
%                                 n_fwage(t,i)=n_fwage(t-1,i);
%                             end
%                             m_fwage(t,i)=m_fwage(t-1,i);
%                         end
%                     end
%                 end
%             else %Nothing happens at SM
%                 firm_status(t,i)=4;
%                 m_ftp(t,i)=m_ftp(t-1,i);
%                 m_fwage(t,i)=m_fwage(t-1,i);
%                 n_fwage(t,i)=n_fwage(t-1,i);
%             end
%         end
        
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %% Reallocation post SM
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % fprintf('R for status 2')
%         if (firm_status(t,i)==2) %Gets to the end of SM with a manager
%             if (d_m(a_ftp(t,i),m_ftp(t,i))>0)  %Fire the manager
%                 firm_status(t,i)=1;
%                 fire_m(t,i)=1;
%                 m_ftp(t,i)=0;
%                 n_ftp(t,i)=0;
%                 m_fwage(t,i)=0.0;
%                 n_fwage(t,i)=0.0;
%             elseif (r_m(a_ftp(t,i),m_ftp(t,i))>0) %reallocate the manager
%                 firm_status(t,i)=3;
%                 prom_sm(t,i)=-1; %Demote the manager
%                 n_ftp(t,i)=m_ftp(t,i);
%                 n_fwage(t,i)=m_fwage(t,i); %Wage does not change upon reallocation
%                 m_ftp(t,i)=0;
%                 m_fwage(t,i)=0.0;
%             else %Manager stays
%                 firm_status(t,i)=2;
%                 %Does the wage goes down?
%                 current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wm(:,a_ftp(t,i),m_ftp(t,i)));
%                 target_wage=Vm(a_ftp(t,i),m_ftp(t,i))-Ve(a_ftp(t,i));
%                 if (target_wage<current_wage)
%                     m_fwage(t,i)= InterpolateWage(Wm(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
%                 end
%                 n_ftp(t,i)=0;
%                 n_fwage(t,i)=0.0;
%             end
%         end

%         % fprintf('R for status 3')
%         if (firm_status(t,i)==3) %Gets to the end of SM with a non manager
%             %Reallocation post SM
%             if (d_n(a_ftp(t,i),n_ftp(t,i))>0) %Fire the non manager
%                 firm_status(t,i)=1;
%                 fire_n(t,i)=1;
%                 m_ftp(t,i)=0;
%                 n_ftp(t,i)=0;
%                 m_fwage(t,i)=0.0;
%                 n_fwage(t,i)=0.0;
%             elseif (r_n(a_ftp(t,i),n_ftp(t,i))>0) %reallocate the non manager
%                 firm_status(t,i)=2; %Firm with manager
%                 prom_sm(t,i)=1; %Promote the non manager
%                 m_ftp(t,i)=n_ftp(t,i);
%                 m_fwage(t,i)=n_fwage(t,i); %Wage does not change upon reallocation
%                 n_ftp(t,i)=0;
%                 n_fwage(t,i)=0.0;
%             else %Non manager stays
%                 firm_status(t,i)=3;
%                 %Does the wage goes down?
%                 current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wn(:,a_ftp(t,i),n_ftp(t,i)));
%                 target_wage=Vn(a_ftp(t,i),n_ftp(t,i))-Ve(a_ftp(t,i));
%                 if (target_wage<current_wage)
%                     n_fwage(t,i)= InterpolateWage(Wn(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
%                 end
%                 m_ftp(t,i)=0;
%                 m_fwage(t,i)=0.0;
%             end
%         end
        
%         % fprintf('R for status 4')
%         if (firm_status(t,i)==4) % Gets to the end of Sm with a team 
%             if (d_t_b(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0) %Fire the team
%                 firm_status(t,i)=1;
%                 fire_m(t,i)=1;
%                 fire_n(t,i)=1;
%                 m_ftp(t,i)=0;
%                 n_ftp(t,i)=0;
%                 m_fwage(t,i)=0.0;
%                 n_fwage(t,i)=0.0;
%             elseif (r_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))*r_t_n(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0) %reallocate the team
%                 %Lets be careful here to ensure the swap
%                 firm_status(t,i)=4;
%                 prom_realloc(t,i)=1; %Promote the non manager
%                 temp_m=m_ftp(t,i);
%                 temp_n=n_ftp(t,i);
%                 temp_mw=m_fwage(t,i);
%                 temp_nw=n_fwage(t,i);
%                 m_ftp(t,i)=temp_n;
%                 n_ftp(t,i)=temp_m;
%                 m_fwage(t,i)=temp_nw;
%                 n_fwage(t,i)=temp_mw;
%                 clear temp_m temp_n temp_mw temp_nw;
%             elseif (d_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))*(1-r_t_n(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)))>0) %Fire the manager do not reallocate nm
%                 firm_status(t,i)=3;
%                 fire_m(t,i)=1;
%                 m_ftp(t,i)=0;
%                 m_fwage(t,i)=0.0;
%             elseif (d_t_n(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))*(1-r_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)))>0) %Fire the non manager do not reallocate m
%                 firm_status(t,i)=2;
%                 fire_n(t,i)=1;
%                 n_ftp(t,i)=0;
%                 n_fwage(t,i)=0.0;
%             elseif (d_t_n(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))*r_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0) %Fire the non manager reallocate m
%                 firm_status(t,i)=3;
%                 fire_n(t,i)=1;
%                 prom_realloc(t,i)=-1; %Demote the manager
%                 n_ftp(t,i)=m_ftp(t,i);
%                 n_fwage(t,i)=m_fwage(t,i); %Wage does not change upon reallocation
%                 m_ftp(t,i)=0;
%                 m_fwage(t,i)=0.0;
%             elseif (d_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))*r_t_n(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0) %Fire the manager reallocate nm
%                 firm_status(t,i)=2;
%                 fire_m(t,i)=1;
%                 prom_realloc(t,i)=1; %Promote the non manager
%                 m_ftp(t,i)=n_ftp(t,i);
%                 m_fwage(t,i)=n_fwage(t,i); %Wage does not change upon reallocation
%                 n_ftp(t,i)=0;
%                 n_fwage(t,i)=0.0;
%             else %Team stays
%                 %Does the wage goes down?
%                 %For the manager
%                 current_wage_m=InterpolateWage(wgrid, m_fwage(t,i),Wtm(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                 target_wage_m=Vt(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))-Vn(a_ftp(t,i),n_ftp(t,i));
%                 if (target_wage_m<current_wage_m)
%                     m_fwage(t,i)= InterpolateWage(Wtm(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage_m,wgrid);
%                 end
%                 %For the non manager
%                 current_wage_n=InterpolateWage(wgrid, n_fwage(t,i),Wtn(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
%                 target_wage_n=Vt(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))-Vm(a_ftp(t,i),m_ftp(t,i));
%                 if (target_wage_n<current_wage_n)
%                     n_fwage(t,i)= InterpolateWage(Wtn(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage_n,wgrid);
%                 end
%             end
%         end
%     end
% end
% toc
        


% %% Simulation for Workers
% n_months=10;
% n_workers=30;
% max_age=12*200;
% for i=1:12*200
%     cdf_age_dist(i) = 1.0  - (1.0 - death)^i;
% end
% cdf_age_dist(max_age) = 1.0; %age distribution (for inital draw)


% %Qunatiles of wage dist to be used in the simulation for cowarkers wage
% nquantiles=10;
% % quantiles=zeros(nquantiles,tpts,tpts);

% %From the firm simulation, get cdf distributions of wage for each type of firm
% %Firm with man only (a,z)
% quantiles_m=zeros(nquantiles,ats,tpts);
% for a=1:ats
%     for z=1:tpts
%         I=find(firm_status==2 & a_ftp==a & m_ftp==z);
%         if isempty(I)==0
%             quantiles_m(:,a,z) =quantile(m_fwage(I),nquantiles);
%         end
%     end
% end

% %quantiles_n=zeros(nquantiles,ats,tpts);
% for a=1:ats
%     for z=1:tpts
%         I=find(firm_status==3 & a_ftp==a & n_ftp==z);
%         if isempty(I)==0
%             quantiles_n(:,a,z) =quantile(n_fwage(I),nquantiles);
%         end
%     end
% end

% %Quantiles of manager in teams 
% quantiles_tm=zeros(nquantiles,ats,tpts,tpts);
% quantiles_tn=zeros(nquantiles,ats,tpts,tpts);
% for a=1:ats
%     for z=1:tpts
%         for n=1:tpts
%             Im=find(firm_status==4 & a_ftp==a & m_ftp==z & n_ftp==n);
%             if isempty(I)==0
%                 quantiles_tm(:,a,z,n) =quantile(m_fwage(Im),nquantiles);
%             end
%             In=find(firm_status==4 & a_ftp==a & m_ftp==z & n_ftp==n);
%             if isempty(I)==0
%                 quantiles_tn(:,a,z,n) =quantile(n_fwage(In),nquantiles);
%             end
%         end
%     end
% end

% cdf_type_birth=measure_to_cdf(typebirth);
        
% %Random matrices
% rng(1,"twister");

% r.state=rand(n_months,n_workers);
% r.type=rand(n_months,n_workers);
% r.event=rand(n_months,n_workers);
% r.type=rand(n_months,n_workers);
% r.ashock=rand(n_months,n_workers);
% r.qshock=rand(n_months,n_workers);
% r.ushock=rand(n_months,n_workers);
% r.death=rand(n_months,n_workers);
% r.tie_break=rand(n_months,n_workers);

% %Store the results
% worker_status=zeros(n_months,n_workers); %(1=unemployed, 2=employed Manager alone, 3=employed Non Manager alone, 4=employed in a team as a manager, 5=employed in a team as a non manager)
% self_z=zeros(n_months,n_workers); %Workers own productivity
% co_z=zeros(n_months,n_workers); %Productivity of the co-worker
% a_wtp=zeros(n_months,n_workers); %Firm type of the worker 
% self_wage=zeros(n_months,n_workers); %Wage of the worker
% co_wage=zeros(n_months,n_workers); %Wage of the co-worker
% self_age=zeros(n_months,n_workers); %Age of the worker
% self_tenure=zeros(n_months,n_workers); %Tenure of the worker
% co_tenure=zeros(n_months,n_workers); %Tenure of the co-worker
% hire_m=zeros(n_months,n_workers); %Hired as the manager
% hire_n=zeros(n_months,n_workers); %Hired as the non manager
% fire_m=zeros(n_months,n_workers); %Fired as the manager
% fire_n=zeros(n_months,n_workers); %Fired as the non manager
% prom_sm=zeros(n_months,n_workers); %Promoted or demoted as the manager at the SM
% prom_realloc=zeros(n_months,n_workers); %Promoted or demoted as the manager at the reallocation
% fire_m_realloc=zeros(n_months,n_workers); %Fired as the manager at the reallocation
% fire_n_realloc=zeros(n_months,n_workers); %Fired as the non manager at the reallocation


% tic
% %Initial conditions
% t=1;
% for i=1:n_workers
%     if (r.state(t,i) < cdf_emp_dist(1)) %Unemployed
%         worker_status(t,i)=1;
%         zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
%         self_z(t,i)=zu;
%         self_wage(t,i)=0.0;
%         self_age(t,i)=draw_CDF_1d(cdf_age_dist,r.type(t,i));
%         self_tenure(t,i)=1;
%         co_z(t,i)=0.0;
%         co_wage(t,i)=0.0;
%         co_tenure(t,i)=0;
%         a_wtp(t,i)=0;
%     elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<cdf_emp_dist(2)) %Employed as a manager alone
%         [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%         worker_status(t,i)=2;
%         self_z(t,i)=zm;
%         self_wage(t,i)=InterpolateWage(Wm(:,am,zm),Vm(am,zm)-Ve(am),wgrid);
%         self_age(t,i)=draw_CDF_1d(cdf_age_dist,r.type(t,i));
%         self_tenure(t,i)=1;
%         a_wtp(t,i)=am;
%         co_z(t,i)=0.0;
%         co_wage(t,i)=0.0;
%         co_tenure(t,i)=0;
%     elseif (r.state(t,i)>=cdf_emp_dist(2) && r.state(t,i)<cdf_emp_dist(3)) %Employed as a non manager alone
%         [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%         worker_status(t,i)=3;
%         self_z(t,i)=zn;
%         self_wage(t,i)=InterpolateWage(Wn(:,an,zn),Vn(an,zn)-Ve(an),wgrid);
%         self_age(t,i)=draw_CDF_1d(cdf_age_dist,r.type(t,i));
%         self_tenure(t,i)=1;
%         a_wtp(t,i)=an;
%         co_z(t,i)=0.0;
%         co_wage(t,i)=0.0;
%         co_tenure(t,i)=0;
%     else %Employed in a team
%         [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%         %Flip a coin to see if the worker is the manager or the non manager
%         if (r.event(t,i)<0.5) %As a manager
%             worker_status(t,i)=4;
%             self_z(t,i)=zt;
%             self_wage(t,i)=InterpolateWage(Wtm(:,at,zt,nt),Vt(at,zt,nt)-Vn(at,nt),wgrid);
%             self_age(t,i)=draw_CDF_1d(cdf_age_dist,r.type(t,i));
%             self_tenure(t,i)=1;
%             a_wtp(t,i)=at;
%             co_z(t,i)=nt;
%             co_wage(t,i)=InterpolateWage(Wtn(:,at,zt,nt),Vt(at,zt,nt)-Vm(at,zt),wgrid);
%             co_tenure(t,i)=1;
%         else %As a non manager
%             worker_status(t,i)=5;
%             self_z(t,i)=nt;
%             self_wage(t,i)=InterpolateWage(Wtn(:,at,zt,nt),Vt(at,zt,nt)-Vm(at,zt),wgrid);
%             self_age(t,i)=draw_CDF_1d(cdf_age_dist,r.type(t,i));
%             self_tenure(t,i)=1;
%             a_wtp(t,i)=at;
%             co_z(t,i)=zt;
%             co_wage(t,i)=InterpolateWage(Wtm(:,at,zt,nt),Vt(at,zt,nt)-Vn(at,nt),wgrid);
%             co_tenure(t,i)=1;
%         end
%     end
% end 

% %For other months, lets follow the amazing journey of the workers
% for i=1:n_workers
%     for t=2:n_months
%         %% SM for workers
%         % fprintf('SM Status 1')
%         if (worker_status(t-1,i)==1) %unemployed workers
%             %Check if they die
%             if (r.death(t,i)<death | self_age(t-1,i) >= max_age)
%                 worker_status(t,i)=1;
%                 %new draw but directly on unemployed
%                 self_z(t,i)=draw_CDF_1d(cdf_type_birth,r.type(t,i));
%                 self_wage(t,i)=0.0;
%                 self_age(t,i)=0;
%                 self_tenure(t,i)=0;
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 a_wtp(t,i)=0;
%             else %Does not die
%                 self_age(t,i)=self_age(t-1,i)+1;
%                 %Check if lose human capital
%                 if (r.ushock(t,i)<down_u(self_z(t-1,i)))
%                     self_z(t,i)=self_z(t-1,i)-1;
%                 % elseif (r.ushock(t,i)>up_u(self_z(t-1,i)))
%                 %     self_z(t,i)=self_z(t-1,i)+1;
%                 else
%                     self_z(t,i)=self_z(t-1,i);
%                 end
%                 %Events 
%                 if (r.event(t,i)<lamu) %Meet another firm 
%                     if (r.state(t,i)<cdf_firm_dist(1))%meet vacant firm
%                         av=draw_CDF_1d(cdf_e_edist,r.type(t,i));
%                         %Does get hired?
%                         if (h_e_u(self_z(t,i),av)>0)
%                             %Where does it go?
%                             if (p_e_m(av,self_z(t,i))>0) %Hired as a manager
%                                 worker_status(t,i)=2;
%                                 a_wtp(t,i)=av;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vmh(av,self_z(t,i))-Veh(av));
%                                 self_wage(t,i)= InterpolateWage(Wmh(:,av,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             else %Hired as non-manager
%                                 worker_status(t,i)=3;
%                                 a_wtp(t,i)=av;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vnh(av,self_z(t,i))-Veh(av));
%                                 self_wage(t,i)= InterpolateWage(Wnh(:,av,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=1;
%                             self_wage(t,i)=0.0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %Unemployed for one more month
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             a_wtp(t,i)=0;
%                         end
%                     elseif (r.state(t,i)>=cdf_firm_dist(1) && r.state(t,i)<cdf_firm_dist(2)) %Meet firm with manager
%                         [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_m_u(self_z(t,i),am,zm)>0)
%                             %Where does it go?
%                             if (p_m_m_u(am,zm,self_z(t,i))>0) %Hired as a manager firing curent one
%                                 worker_status(t,i)=2;
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vmh(am,self_z(t,i))+ U(zm) -Vmh(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wmh(:,am,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             elseif (p_m_m_d(am,zm,self_z(t,i))>0) %Hired as a manager demoting the current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vth(am,self_z(t,i),zm)-cost_d-Vm(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,am,zm,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                                 co_z(t,i)=zm;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_m(k,am,zm);
%                                 co_tenure(t,i)=1;
%                             else %Hired as a non manager
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vth(am,zm,self_z(t,i))-Vm(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,am,zm,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=zm;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_m(k,am,zm);
%                                 co_tenure(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=1;
%                             self_wage(t,i)=0.0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %Unemployed for one more month
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             a_wtp(t,i)=0;
%                         end
%                     elseif (r.state(t,i)>=cdf_firm_dist(2) && r.state(t,i)<cdf_firm_dist(3)) %Meet firm with non manager
%                         [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_nm_u(self_z(t,i),an,zn)>0)
%                             %Where does it go?
%                             if (p_n_n_u(an,zn,self_z(t,i))>0) %Hired as a non manager firing curent one
%                                 worker_status(t,i)=3;
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vnh(an,self_z(t,i))+U(zn)-Vnh(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wnh(:,an,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             elseif (p_n_n_p(an,zn,self_z(t,i))>0) %Hired as a non manager promoting the current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vth(an,zn,self_z(t,i))-cost_p-Vn(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,an,zn,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                                 co_z(t,i)=zn;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_n(k,an,zn);
%                                 co_tenure(t,i)=1;
%                             else %Hired as a manager
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vth(an,self_z(t,i),zn)-Vn(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,an,zn,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=zn;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_n(k,an,zn);
%                                 co_tenure(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=1;
%                             self_wage(t,i)=0.0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %Unemployed for one more month
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             a_wtp(t,i)=0;
%                         end
%                     else %Meet firm with team
%                         [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_t_u(self_z(t,i),at,zt,nt)>0)
%                             %Where does it go?
%                             if (p_t_m_u(at,zt,nt,self_z(t,i))>0) %Hired as a manager firing curent one
%                                 worker_status(t,i)=4;
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vth(at,self_z(t,i),nt)+U(zt)-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,at,self_z(t,i),nt),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=nt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tn(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             elseif (p_t_m_d(at,zt,nt,self_z(t,i))>0) %Hired as a manager demoting the current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vth(at,self_z(t,i),zt)+U(nt)-cost_d-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,at,self_z(t,i),zt),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                                 co_z(t,i)=zt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tm(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             elseif (p_t_n_u(at,zt,nt,self_z(t,i))>0) %Hired as a non manager firing curent one
%                                 worker_status(t,i)=5;
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vth(at,zt,self_z(t,i))+U(nt)-Vth(at,nt,zt));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,at,zt,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=zt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tm(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             else %Hired as a non manager promoting the current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*U(self_z(t,i))+bpw*(Vth(at,nt,self_z(t,i))+U(zt)-cost_p-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,at,nt,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                                 co_z(t,i)=nt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tn(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=1;
%                             self_wage(t,i)=0.0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %Unemployed for one more month
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             a_wtp(t,i)=0;
%                         end
%                     end
%                 else %Does not meet a firm
%                     worker_status(t,i)=1;
%                     self_wage(t,i)=0.0;
%                     self_tenure(t,i)=self_tenure(t-1,i)+1; %Unemployed for one more month
%                     co_z(t,i)=0.0;
%                     co_wage(t,i)=0.0;
%                     co_tenure(t,i)=0;
%                     a_wtp(t,i)=0;
%                 end
%             end
%         end
%         % fprintf('SM Status 2')
%         if (worker_status(t-1,i)==2) %Manager alone
%             %Check productivity shock
%             if (r.ashock(t,i)<down_a(a_wtp(t-1,i))) %Firm go down in a
%                 a_wtp(t,i)=a_wtp(t-1,i)-1;
%             elseif (r.ashock(t,i)< down_a(a_wtp(t-1,i))+stay_a(a_wtp(t-1,i)) && r.ashock(t,i)>=down_a(a_wtp(t-1,i))) %Firm stay in a
%                 a_wtp(t,i)=a_wtp(t-1,i);
%             else %Firm go up in a
%                 a_wtp(t,i)=a_wtp(t-1,i)+1;
%             end
%             %Check if they die
%             if (r.death(t,i)<death) 
%                 worker_status(t,i)=1;
%                 %new draw but directly on unemployed
%                 self_z(t,i)=draw_CDF_1d(cdf_type_birth,r.type(t,i));
%                 self_wage(t,i)=0.0;
%                 self_age(t,i)=0;
%                 self_tenure(t,i)=0;
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 a_wtp(t,i)=0;
%             else %Does not die
%                 self_age(t,i)=self_age(t-1,i)+1;
%                 %No learnig because it is a manager
%                 self_z(t,i)=self_z(t-1,i);
%                 self_wage(t,i)=self_wage(t-1,i); %Wage does not change
%                 %Events
%                 if (r.event(t,i)<=del) %Leave the firm exogenoulsy
%                     worker_status(t,i)=1;
%                     self_wage(t,i)=0.0;
%                     self_tenure(t,i)=0;
%                     co_z(t,i)=0.0;
%                     co_wage(t,i)=0.0;
%                     co_tenure(t,i)=0;
%                     a_wtp(t,i)=0;
%                 elseif (r.event(t,i)>del && r.event(t,i)<=del+prob_matching) %Firm meets another worker
%                     if (r.state(t,i)<cdf_emp_dist(1))%meet unemployed
%                         zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
%                         %Does get hired?
%                         if (h_m_u(zu,a_wtp(t,i),self_z(t,i))>0)
%                             %Where does it go?
%                             if (p_m_m_u(a_wtp(t,i),self_z(t,i),zu)>0) %Hired as a manager firing curent one (me)
%                                 worker_status(t,i)=1; 
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             elseif (p_m_m_d(a_wtp(t,i),self_z(t,i),zu)>0) %Hired as a manager demoting the current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zu;
%                                 target_wage= (1-bpw)*U(zu)+bpw*(Vth(a_wtp(t,i),zu,self_z(t,i))-cost_d-Vmh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zu,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                             else %Hired as a non manager
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zu;
%                                 target_wage= (1-bpw)*U(zu)+bpw*(Vth(a_wtp(t,i),self_z(t,i),zu)-Vmh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zu),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=2;    
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             a_wtp(t,i)=0;
%                         end
%                     elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<cdf_emp_dist(2)) %meet firm with manager
%                         [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_m_m(am,zm,a_wtp(t,i),self_z(t,i))>0)
%                             %Where does it go?
%                             if (p_m_m_u(a_wtp(t,i),self_z(t,i),zm)>0) %Hired as a manager firing curent one (me)
%                                 worker_status(t,i)=1; 
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             elseif (p_m_m_d(a_wtp(t,i),self_z(t,i),zm)>0) %Hired as a manager demoting the current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zm;
%                                 target_wage= (1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_wtp(t,i),zm,self_z(t,i))-cost_d-Vmh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zm,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                             else %Hired as a non manager
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zm;
%                                 target_wage= (1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zm)- Vmh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zm),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=2;
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             a_wtp(t,i)=0;
%                         end
%                     elseif (r.state(t,i)>=cdf_emp_dist(2) && r.state(t,i)<cdf_emp_dist(3)) %meet firm with non manager
%                         [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_m_nm(an,zn,a_wtp(t,i),self_z(t,i))>0)
%                             %Where does it go?
%                             if (p_m_m_u(a_wtp(t,i),self_z(t,i),zn)>0) %Hired as a manager firing curent one (me)
%                                 worker_status(t,i)=1; 
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             elseif (p_m_m_d(a_wtp(t,i),self_z(t,i),zn)>0) %Hired as a manager demoting the current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zn;
%                                 target_wage= (1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_wtp(t,i),zn,self_z(t,i))-cost_d-Vmh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zn,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                             else %Hired as a non manager
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zn;
%                                 target_wage= (1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zn)- Vmh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zn),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=2;
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             a_wtp(t,i)=0;
%                         end
%                     else %meet firm with team
%                         [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                         %Does hire the maneger? 
%                         if (h_m_t_m(at,zt,nt,a_wtp(t,i),self_z(t,i)>0))
%                             %Where does it go?
%                             if (p_m_m_u(a_wtp(t,i),self_z(t,i),zt)>0) %Hired as a manager firing curent one (me)
%                                 worker_status(t,i)=1; 
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             elseif (p_m_m_d(a_wtp(t,i),self_z(t,i),zt)>0) %Hired as a manager demoting the current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zt;
%                                 target_wage= (1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_wtp(t,i),zt,self_z(t,i))-cost_d-Vmh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zt,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                             else %Hired as a non manager
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zt;
%                                 target_wage= (1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zt)- Vmh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zt),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                             end
%                         %Does hire the non maneger?    
%                         elseif (h_m_t_nm(at,zt,nt,a_wtp(t,i),self_z(t,i)>0))
%                             %Where does it go?
%                             if (p_m_m_u(a_wtp(t,i),self_z(t,i),nt)>0) %Hired as a manager firing curent one (me)
%                                 worker_status(t,i)=1; 
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             elseif (p_m_m_d(a_wtp(t,i),self_z(t,i),nt)>0) %Hired as a manager demoting the current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=nt;
%                                 target_wage= (1-bpw)*(Vth(at,nt,zt)-Vmh(at,zt))+bpw*(Vth(a_wtp(t,i),nt,self_z(t,i))-cost_d-Vmh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),nt,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                             else %Hired as a non manager
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=nt;
%                                 target_wage= (1-bpw)*(Vth(at,nt,zt)-Vmh(at,zt))+bpw*(Vth(a_wtp(t,i),self_z(t,i),nt)- Vmh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),nt),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=2;
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             a_wtp(t,i)=0;
%                         end
%                     end
%                 elseif (r.event(t,i)>del+prob_matching && r.event(t,i)<=del+prob_matching+lam) %Workers meets another firm 
%                     if (r.state(t,i)<=cdf_firm_dist(1))%meet vacant firm
%                         av=draw_CDF_1d(cdf_e_edist,r.type(t,i));
%                         %Does get hired?
%                         if (h_e_m(a_wtp(t,i),self_z(t,i),av)>0)
%                             %Where does it go?
%                             if (p_e_m(av,self_z(t,i))>0) %Hired as a manager
%                                 worker_status(t,i)=2;
%                                 a_wtp(t,i)=av;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vmh(av,self_z(t,i))-Veh(av));
%                                 self_wage(t,i)= InterpolateWage(Wmh(:,av,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             else %Hired as non-manager
%                                 worker_status(t,i)=3;
%                                 a_wtp(t,i)=av;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vnh(av,self_z(t,i))-Veh(av));
%                                 self_wage(t,i)= InterpolateWage(Wnh(:,av,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=2;
%                             %Does the wage increase?
%                             current_wage=InterpolateWage(wgrid,self_wage(t,i),Wmh(:,a_wtp(t,i),self_z(t,i)));
%                             target_wage= max(Vmh(av,self_z(t,i)),Vnh(av,self_z(t,i)))-Veh(av);
%                             if (current_wage<target_wage)
%                                 self_wage(t,i)=InterpolateWage(Wmh(:,a_wtp(t,i),self_z(t,i)),target_wage,wgrid);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     elseif (r.state(t,i)>=cdf_firm_dist(1) && r.state(t,i)<cdf_firm_dist(2)) %meet firm with manager
%                         [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_m_m(a_wtp(t,i),self_z(t,i),am,zm)>0)
%                             %Where does it go?
%                             if (p_m_m_u(am,zm,self_z(t,i))>0) %Hired as a manager firing curent one
%                                 worker_status(t,i)=2; %manager alone
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vmh(am,self_z(t,i))-Veh(am));
%                                 self_wage(t,i)= InterpolateWage(Wmh(:,am,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             elseif (p_m_m_d(am,zm,self_z(t,i))>0) %Hired as a manager demoting the current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(am,self_z(t,i),zm)-cost_d-Vm(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,am,self_z(t,i),zm),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                                 co_z(t,i)=zm;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_m(k,am,zm);
%                                 co_tenure(t,i)=1;
%                             else %Hired as a non manager
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(am,zm,self_z(t,i))-Vm(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,am,zm,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=zm;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_m(k,am,zm);
%                                 co_tenure(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=2;
%                             %Does the wage increase?
%                             current_wage=InterpolateWage(wgrid,self_wage(t,i),Wmh(:,a_wtp(t,i),self_z(t,i)));
%                             nu=[Vth(am,zm,self_z(t,i)),Vmh(am,self_z(t,i))+U(zm),Vth(am,self_z(t,i),zm)-cost_d]
%                             target_wage= max(nu)-Vmh(am,zm);
%                             if (current_wage<target_wage)
%                                 self_wage(t,i)=InterpolateWage(Wmh(:,a_twp(t,i),self_z(t,i)),target_wage,wgrid);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     elseif (r.state(t,i)>=cdf_firm_dist(2) && r.state(t,i)<cdf_firm_dist(3)) %meet firm with non manager
%                         [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_nm_m(a_wtp(t,i),self_z(t,i),an,zn)>0)
%                             %Where does it go?
%                             if (p_n_n_u(an,zn,self_z(t,i))>0) %Hired as a non manager firing curent one
%                                 worker_status(t,i)=3;
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vnh(an,self_z(t,i))+U(zn)-Vnh(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wnh(:,an,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             elseif (p_n_n_p(an,zn,self_z(t,i))>0) %Hired as a non manager promoting the current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(an,zn,self_z(t,i))-cost_p-Vn(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,an,zn,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                                 co_z(t,i)=zn;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_n(k,an,zn);
%                                 co_tenure(t,i)=1;
%                             else %Hired as a manager
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(an,self_z(t,i),zn)-Vn(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,an,zn,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=zn;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_n(k,an,zn);
%                                 co_tenure(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=2;
%                             %Does the wage increase?
%                             current_wage=InterpolateWage(wgrid,self_wage(t,i),Wmh(:,a_wtp(t,i),self_z(t,i)));
%                             nu=[Vnh(an,self_z(t,i))+U(zn),Vth(an,zn,self_z(t,i))-cost_p,Vth(an,self_z(t,i),zn)]
%                             target_wage= max(nu)-Vmh(an,zn);
%                             if (current_wage<target_wage)
%                                 self_wage(t,i)=InterpolateWage(Wmh(:,a_twp(t,i),self_z(t,i)),target_wage,wgrid);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     else %meet firm with team
%                         [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_t_m(a_wtp(t,i),self_z(t,i),at,zt,nt)>0)
%                             %Where does it go?
%                             if (p_t_m_u(at,zt,nt,self_z(t,i))>0) %Hired as a manager firing curent one
%                                 worker_status(t,i)=4;
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(at,self_z(t,i),nt)+U(zt)-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,at,self_z(t,i),nt),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=nt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tn(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             elseif (p_t_m_d(at,zt,nt,self_z(t,i))>0) %Hired as a manager demoting the current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(at,self_z(t,i),zt)+U(nt)-cost_d-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,at,self_z(t,i),zt),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                                 co_z(t,i)=zt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tm(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             elseif (p_t_n_u(at,zt,nt,self_z(t,i))>0) %Hired as a non manager firing curent one
%                                 worker_status(t,i)=5;
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(at,zt,self_z(t,i))+U(nt)-Vth(at,nt,zt));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,at,zt,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=zt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tm(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             else %Hired as a non manager promoting the current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vmh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(at,nt,self_z(t,i))+U(zt)-cost_p-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,at,nt,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                                 co_z(t,i)=nt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tn(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=2;
%                             %Does the wage increase?
%                             current_wage=InterpolateWage(wgrid,self_wage(t,i),Wmh(:,a_wtp(t,i),self_z(t,i)));
%                             nu=[Vth(at,zt,self_z(t,i))+U(nt),Vth(at,nt,self_z(t,i))-cost_p+U(zt),Vth(at,self_z(t,i),zt)-cost_d+U(nt),Vth(at,self_z(t,i),nt)+U(zt)]
%                             target_wage= max(nu)-Vth(at,zt,nt);
%                             if (current_wage<target_wage)
%                                 self_wage(t,i)=InterpolateWage(Wmh(:,a_twp(t,i),self_z(t,i)),target_wage,wgrid);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     end
%                 else %Does not meet a firm
%                     worker_status(t,i)=2;
%                     self_wage(t,i)=self_wage(t-1,i); %Wage does not change
%                     self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                     co_z(t,i)=0.0;
%                     co_wage(t,i)=0.0;
%                     co_tenure(t,i)=0;
%                 end
%             end
%         end
%         % fprintf('SM Status 3')
%         if (worker_status(t-1,i)==3) %Non Manager alone
%             %Check productivity shock
%             if (r.ashock(t,i)<down_a(a_wtp(t-1,i))) %Firm go down in a
%                 a_wtp(t,i)=a_wtp(t-1,i)-1;
%             elseif (r.ashock(t,i)< down_a(a_wtp(t-1,i))+stay_a(a_wtp(t-1,i)) && r.ashock(t,i)>=down_a(a_wtp(t-1,i))) %Firm stay in a
%                 a_wtp(t,i)=a_wtp(t-1,i);
%             else %Firm go up in a
%                 a_wtp(t,i)=a_wtp(t-1,i)+1;
%             end
%             %Check if they die
%             if (r.death(t,i)<death) 
%                 worker_status(t,i)=1;
%                 %new draw but directly on unemployed
%                 self_z(t,i)=draw_CDF_1d(cdf_type_birth,r.type(t,i));
%                 self_wage(t,i)=0.0;
%                 self_age(t,i)=0;
%                 self_tenure(t,i)=0;
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 a_wtp(t,i)=0;
%             else %Does not die
%                 self_age(t,i)=self_age(t-1,i)+1;
%                 %Check if learns
%                 if (r.qshock(t,i)<down_q(a_wtp(t-1,i),self_z(t-1,i))) %Firm go down in q
%                     self_z(t,i)=self_z(t-1,i)-1;    
%                 elseif (r.qshock(t,i)< down_q(a_wtp(t-1,i),self_z(t-1,i))+stay_q(a_wtp(t-1,i),self_z(t-1,i)) && r.qshock(t,i)>=down_q(a_wtp(t-1,i),self_z(t-1,i)))%Firm stay in q
%                     self_z(t,i)=self_z(t-1,i);
%                 else %Firm go up in q
%                     self_z(t,i)=self_z(t-1,i)+1;
%                 end
%                 % if (r.qshock(t,i) < down_q(self_z(t-1,i),a_wtp(t-1,i)))
%                 %     self_z(t,i)=self_z(t-1,i)-1;
%                 % elseif (r.qshock(t,i) > up_q(self_z(t-1,i),a_wtp(t-1,i)))
%                 %     self_z(t,i)=self_z(t-1,i)+1;
%                 % else
%                 %     self_z(t,i)=self_z(t-1,i);
%                 % end
%                 self_wage(t,i)=self_wage(t-1,i); %Wage does not change
%                 %Events
%                 if (r.event(t,i)<=del) %Leave the firm exogenoulsy
%                     worker_status(t,i)=1;
%                     self_wage(t,i)=0.0;
%                     self_tenure(t,i)=0;
%                     co_z(t,i)=0.0;
%                     co_wage(t,i)=0.0;
%                     co_tenure(t,i)=0;
%                     a_wtp(t,i)=0;
%                 elseif (r.event(t,i)>del && r.event(t,i)<=del+prob_matching) %Firm meets another worker
%                     if (r.state(t,i)<=cdf_emp_dist(1))%meet unemployed
%                         zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
%                         %Does get hired?
%                         if (h_nm_u(zu,a_wtp(t,i),self_z(t,i))>0)
%                             %Where does it go?
%                             if (p_n_n_u(a_wtp(t,i),self_z(t,i),zu)>0) %Hired as a non manager firing curent one (me)
%                                 worker_status(t,i)=1; 
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_n(t,i)=1;
%                             elseif (p_n_n_p(a_wtp(t,i),self_z(t,i),zu)>0) %Hired as a non manager promoting the current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zu;
%                                 target_wage= (1-bpw)*U(zu)+bpw*(Vth(a_wtp(t,i),self_z(t,i),zu)-cost_p-Vnh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zu),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                             else %Hired as a manager
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zu;
%                                 target_wage= (1-bpw)*U(zu)+bpw*(Vth(a_wtp(t,i),zu,self_z(t,i))-Vnh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zu,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=3; %non manager alone
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<cdf_emp_dist(2)) %meet firm with manager
%                         [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_nm_m(am,zm,a_wtp(t,i),self_z(t,i)>0))
%                             %Where does it go?
%                             if (p_n_n_u(a_wtp(t,i),self_z(t,i),zm)>0) %Hired as a non manager firing curent one (me)
%                                 worker_status(t,i)=1; 
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_n(t,i)=1;
%                             elseif (p_n_n_p(a_wtp(t,i),self_z(t,i),zm)>0) %Hired as a non manager promoting the current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zm;
%                                 target_wage= (1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zm)-cost_p-Vnh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zm),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                             else %Hired as a manager
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zm;
%                                 target_wage= (1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_wtp(t,i),zm,self_z(t,i))-Vnh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zm,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=3; %non manager alone
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     elseif (r.state(t,i)>=cdf_emp_dist(2) && r.state(t,i)<cdf_emp_dist(3)) %meet firm with non manager
%                         [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_nm_nm(an,zn,a_wtp(t,i),self_z(t,i))>0)
%                             %Where does it go?
%                             if (p_n_n_u(a_wtp(t,i),self_z(t,i),zn)>0) %Hired as a non manager firing curent one (me)
%                                 worker_status(t,i)=1; 
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_n(t,i)=1;
%                             elseif (p_n_n_p(a_wtp(t,i),self_z(t,i),zn)>0) %Hired as a non manager promoting the current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zn;
%                                 target_wage= (1-bpw)*(Vmh(an,zn)-Veh(an))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zn)-cost_p-Vnh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zn),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                             else %Hired as a manager
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zn;
%                                 target_wage= (1-bpw)*(Vmh(an,zn)-Veh(an))+bpw*(Vth(a_wtp(t,i),zn,self_z(t,i))-Vnh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zn,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=3; %non manager alone
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     else %meet firm with team
%                         [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                         %Hire the maneger?
%                         if (h_nm_t_m(at,zt,nt,a_wtp(t,i),self_z(t,i)>0))
%                             %Where does it go?
%                             if (p_n_n_u(a_wtp(t,i),self_z(t,i),zt)>0) %Hired as a non manager firing curent one (me)
%                                 worker_status(t,i)=1; 
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_n(t,i)=1
%                             elseif (p_n_n_p(a_wtp(t,i),self_z(t,i),zt)>0) %Hired as a non manager promoting the current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zt;
%                                 target_wage= (1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zt)-cost_p-Vnh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zt),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                             else %Hired as a manager
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zt;
%                                 target_wage= (1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_wtp(t,i),zt,self_z(t,i))-Vnh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zt,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                             end
%                         %Does hire the non maneger?    
%                         elseif (h_nm_t_nm(at,zt,nt,a_wtp(t,i),self_z(t,i)>0))
%                             %where does it go?
%                             if (p_n_n_u(a_wtp(t,i),self_z(t,i),nt)>0) %Hired as a non manager firing curent one (me)
%                                 worker_status(t,i)=1; 
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_n(t,i)=1;
%                             elseif (p_n_n_p(a_wtp(t,i),self_z(t,i),nt)>0) %Hired as a non manager promoting the current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=nt;
%                                 target_wage= (1-bpw)*(Vth(at,nt,zt)-Vmh(at,zt))+bpw*(Vth(a_wtp(t,i),self_z(t,i),nt)-cost_p-Vnh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),nt),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                             else %Hired as a manager
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=nt;
%                                 target_wage= (1-bpw)*(Vth(at,nt,zt)-Vmh(at,zt))+bpw*(Vth(a_wtp(t,i),nt,self_z(t,i))-Vnh(a_wtp(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),nt,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=3; %non manager alone
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     end
%                 elseif (r.event(t,i)>del+prob_matching && r.event(t,i)<=del+prob_matching+lam) %Workers meets another firm
%                     if (r.state(t,i)<=cdf_firm_dist(1))%meet vacant firm
%                         av=draw_CDF_1d(cdf_e_edist,r.type(t,i));
%                         %Does get hired?
%                         if (h_e_nm(a_wtp(t,i),self_z(t,i),av)>0)
%                             %Where does it go?
%                             if (p_e_m(a,self_z(t,i))>0) %Hired as a manager
%                                 worker_status(t,i)=2;
%                                 a_wtp(t,i)=av;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vmh(av,self_z(t,i))-Veh(av));
%                                 self_wage(t,i)= InterpolateWage(Wmh(:,av,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             else %Hired as non-manager
%                                 worker_status(t,i)=3;
%                                 a_wtp(t,i)=av;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vnh(av,self_z(t,i))-Veh(av));
%                                 self_wage(t,i)= InterpolateWage(Wnh(:,av,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=3; %non manager alone
%                             %Does the wage increase?
%                             current_wage=InterpolateWage(wgrid,self_wage(t,i),Wnh(:,a_wtp(t,i),self_z(t,i)));
%                             target_wage= max(Vnh(av,self_z(t,i)),Vmh(av,self_z(t,i)))-Veh(av);
%                             if (current_wage<target_wage)
%                                 self_wage(t,i)=InterpolateWage(Wnh(:,a_wtp(t,i),self_z(t,i)),target_wage,wgrid);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     elseif (r.state(t,i)>=cdf_firm_dist(1) && r.state(t,i)<cdf_firm_dist(2)) %meet firm with manager
%                         [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_m_nm(a_wtp(t,i),self_z(t,i),am,zm)>0)
%                             %Where does it go?
%                             if (p_m_m_u(am,zm,self_z(t,i))>0) %Hired as a manager firing curent one
%                                 worker_status(t,i)=2; %manager alone
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vmh(am,self_z(t,i))+U(zm)-Vmh(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wmh(:,am,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             elseif (p_m_m_d(am,zm,self_z(t,i))>0) %Hired as a manager demoting the current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(am,self_z(t,i),zm)-cost_d-Vm(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,am,self_z(t,i),zm),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                                 co_z(t,i)=zm;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_m(k,am,zm);
%                                 co_tenure(t,i)=1;
%                             else %Hired as a non manager
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(am,zm,self_z(t,i))-Vm(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,am,zm,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=zm;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_m(k,am,zm);
%                                 co_tenure(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=3; %non manager alone
%                             %Does the wage increase?
%                             current_wage=InterpolateWage(wgrid,self_wage(t,i),Wmh(:,a_wtp(t,i),self_z(t,i)));
%                             nu=[Vmh(am,self_z(t,i))+U(zm),Vth(am,zm,self_z(t,i)),Vth(am,self_z(t,i),zm)-cost_d];
%                             target_wage= max(nu)-Vmh(am,zm);
%                             if (current_wage<target_wage)
%                                 self_wage(t,i)=InterpolateWage(Wmh(:,a_twp(t,i),self_z(t,i)),target_wage,wgrid);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     elseif (r.state(t,i)>=cdf_firm_dist(2) && r.state(t,i)<cdf_firm_dist(3)) %meet firm with non manager
%                         [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_nm_nm(a_wtp(t,i),self_z(t,i),an,zn)>0)
%                             %Where does it go?
%                             if (p_n_n_u(an,zn,self_z(t,i))>0) %Hired as a non manager firing curent one
%                                 worker_status(t,i)=3; %non manager alone
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vnh(an,self_z(t,i))+U(zn)-Vnh(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wnh(:,an,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             elseif (p_n_n_p(an,zn,self_z(t,i))>0) %Hired as a non manager promoting the current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(an,zn,self_z(t,i))-cost_p-Vn(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,an,zn,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                                 co_z(t,i)=zn;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_n(k,an,zn);
%                                 co_tenure(t,i)=1;
%                             else %Hired as a manager
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(an,self_z(t,i),zn)-Vn(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,an,zn,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=zn;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_n(k,an,zn);
%                                 co_tenure(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=3; %non manager alone
%                             %Does the wage increase?
%                             current_wage=InterpolateWage(wgrid,self_wage(t,i),Wnh(:,a_wtp(t,i),self_z(t,i)));
%                             nu=[Vnh(an,self_z(t,i))+U(zn),Vth(an,self_z(t,i),zn),Vth(an,zn,self_z(t,i))-cost_p];
%                             target_wage= max(nu)-Vnh(an,zn);
%                             if (current_wage<target_wage)
%                                 self_wage(t,i)=InterpolateWage(Wnh(:,a_twp(t,i),self_z(t,i)),target_wage,wgrid);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     else %meet firm with team
%                         [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_t_nm(a_wtp(t,i),self_z(t,i),at,zt,nt)>0)
%                             %Where does it go?
%                             if (p_t_m_u(at,zt,nt,self_z(t,i))>0) %Hired as a manager firing curent one
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(at,self_z(t,i),nt)+U(zt)-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,at,self_z(t,i),nt),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=nt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tn(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             elseif (p_t_m_d(at,zt,nt,self_z(t,i))>0) %Hired as a manager demoting the current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(at,self_z(t,i),zt)+U(nt)-cost_d-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,at,self_z(t,i),zt),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                                 co_z(t,i)=zt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tm(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             elseif (p_t_n_u(at,zt,nt,self_z(t,i))>0) %Hired as a non manager firing curent one
%                                 worker_status(t,i)=5;
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(at,zt,self_z(t,i))+U(nt)-Vth(at,nt,zt));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,at,zt,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=zt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tm(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             else %Hired as a non manager promoting the current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vnh(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i)))+bpw*(Vth(at,nt,self_z(t,i))+U(zt)-cost_p-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,at,nt,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                                 co_z(t,i)=nt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tn(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=3; %non manager alone
%                             %Does the wage increase?
%                             current_wage=InterpolateWage(wgrid,self_wage(t,i),Wnh(:,a_wtp(t,i),self_z(t,i)));
%                             nu=[Vth(at,zt,self_z(t,i))+U(nt),Vth(at,nt,self_z(t,i))-cost_p+U(zt),Vth(at,self_z(t,i),zt)-cost_d+U(nt),Vth(at,self_z(t,i),nt)+U(zt)]
%                             target_wage= max(nu)-Vth(at,zt,nt);
%                             if (current_wage<target_wage)
%                                 self_wage(t,i)=InterpolateWage(Wnh(:,a_twp(t,i),self_z(t,i)),target_wage,wgrid);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                         end
%                     end
%                 else %Does not meet a firm
%                     worker_status(t,i)=3; %non manager alone
%                     self_wage(t,i)=self_wage(t-1,i); %Wage does not change
%                     self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                     co_z(t,i)=0.0;
%                     co_wage(t,i)=0.0;
%                     co_tenure(t,i)=0;
%                 end
%             end
%         end
%         % fprintf('SM Status 4')
%         if (worker_status(t-1,i)==4) %Manager in a team
%             %Check productivity shock
%             if (r.ashock(t,i)<down_a(a_wtp(t-1,i))) %Firm go down in a
%                 a_wtp(t,i)=a_wtp(t-1,i)-1;
%             elseif (r.ashock(t,i)< down_a(a_wtp(t-1,i))+stay_a(a_wtp(t-1,i)) && r.ashock(t,i)>=down_a(a_wtp(t-1,i))) %Firm stay in a
%                 a_wtp(t,i)=a_wtp(t-1,i);
%             else %Firm go up in a
%                 a_wtp(t,i)=a_wtp(t-1,i)+1;
%             end
%             %Check if co-workes learns
%             if (r.qshock(t,i)<down_q(a_wtp(t-1,i),co_z(t-1,i))) %Firm go down in q
%                 co_z(t,i)=co_z(t-1,i)-1;    
%             elseif (r.qshock(t,i)< down_q(a_wtp(t-1,i),co_z(t-1,i))+stay_q(a_wtp(t-1,i),co_z(t-1,i)) && r.qshock(t,i)>=down_q(a_wtp(t-1,i),co_z(t-1,i)))%Firm stay in q
%                 co_z(t,i)=co_z(t-1,i);
%             else %Firm go up in q
%                 co_z(t,i)=co_z(t-1,i)+1;
%             end
%             %Check if it dies
%             if (r.death(t,i)<death)
%                 worker_status(t,i)=1; %Unemployed
%                 %new draw but directly on unemployed
%                 self_z(t,i)=draw_CDF_1d(cdf_type_birth,r.type(t,i));
%                 self_wage(t,i)=0.0;
%                 self_age(t,i)=0;
%                 self_tenure(t,i)=0;
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 a_wtp(t,i)=0;
%             else %Does not die
%                 self_age(t,i)=self_age(t-1,i)+1;
%                 self_z(t,i)=self_z(t-1,i);
%                 self_wage(t,i)=self_wage(t-1,i);
%                 %Events
%                 if (r.event(t,i)<=del) %Leave the firm exogenoulsy
%                     worker_status(t,i)=1;
%                     self_wage(t,i)=0.0;
%                     self_tenure(t,i)=0;
%                     co_z(t,i)=0.0;
%                     co_wage(t,i)=0.0;
%                     co_tenure(t,i)=0;
%                     a_wtp(t,i)=0;
%                 elseif (r.event(t,i)>del && r.event(t,i)<=del+del+death) %Coworker leaves or die
%                     worker_status(t,i)=2; %Manager alone
%                     self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                     self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                     co_z(t,i)=0.0;
%                     co_wage(t,i)=0.0;
%                     co_tenure(t,i)=0;
%                 elseif (r.event(t,i)>del+del+death && r.event(t,i)<=del+del+death+prob_matching) %Firm meets another worker
%                     if (r.state(t,i)<=cdf_emp_dist(1))%meet unemployed
%                         zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
%                         %Does get hired?
%                         if (h_t_u(zu,a_wtp(t,i),self_z(t,i),co_z(t,i))>0)
%                             %Where does it go?
%                             if (p_t_m_u(a_wtp(t,i),self_z(t,i),co_z(t,i),zu)>0 | p_t_n_p(a_wtp(t,i),self_z(t,i),co_z(t,i),zu)>0) %Hired as a manager firing curent one (me) or as non manager with promotion
%                                 worker_status(t,i)=1; %Uneomployed
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             elseif (p_t_n_u(a_wtp(t,i),self_z(t,i),co_z(t,i),zu)>0) %Hired as a non manager firing curent one 
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zu;
%                                 target_wage= (1-bpw)*U(zu)+bpw*(Vth(a_wtp(t,i),self_z(t,i),zu)+U(co_z(t,i))-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zu),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                             else %Hired as a manager with demotion
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zu;
%                                 target_wage= (1-bpw)*U(zu)+bpw*(Vth(a_wtp(t,i),zu,self_z(t,i))-cost_d-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zu,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=4; %manager alone
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<cdf_emp_dist(2)) %meet firm with manager
%                         [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_t_m(am,zm,a_wtp(t,i),self_z(t,i),co_z(t,i))>0)
%                             %Where does it go?
%                             if (p_t_m_u(a_wtp(t,i),self_z(t,i),co_z(t,i),zm)>0 | p_t_n_p(a_wtp(t,i),self_z(t,i),co_z(t,i),zm)>0) %Hired as a manager firing curent one (me) or as non manager with promotion
%                                 worker_status(t,i)=1; %Uneomployed
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             elseif (p_t_n_u(a_wtp(t,i),self_z(t,i),co_z(t,i),zm)>0) %Hired as a non manager firing curent one 
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zm;
%                                 target_wage= (1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zm)+U(co_z(t,i))-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zm),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                             else %Hired as a manager with demotion
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zm;
%                                 target_wage= (1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_wtp(t,i),zm,self_z(t,i))-cost_d+U(co_z(t,i))-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zm,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=4; %manager alone
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     elseif (r.state(t,i)>=cdf_emp_dist(2) && r.state(t,i)<cdf_emp_dist(3)) %meet firm with non manager
%                         [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_t_nm(an,zn,a_wtp(t,i),self_z(t,i),co_z(t,i))>0)
%                             %Where does it go?
%                             if (p_t_m_u(a_wtp(t,i),self_z(t,i),co_z(t,i),zn)>0 | p_t_n_p(a_wtp(t,i),self_z(t,i),co_z(t,i),zn)>0) %Hired as a manager firing curent one (me) or as non manager with promotion
%                                 worker_status(t,i)=1; %Uneomployed
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             elseif (p_t_n_u(a_wtp(t,i),self_z(t,i),co_z(t,i),zn)>0) %Hired as a non manager firing curent one 
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zn;
%                                 target_wage= (1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zn)+U(co_z(t,i))-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zn),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                             else %Hired as a manager with demotion
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zn;
%                                 target_wage= (1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_wtp(t,i),zn,self_z(t,i))-cost_d+U(co_z(t,i))-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i)));   
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zn,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=4; %manager alone
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     else %meet firm with team
%                         [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                         %Hire the manager?
%                         if (h_t_t_m(at,zt,nt,a_wtp(t,i),self_z(t,i),co_z(t,i))>0)
%                             %where does it go?
%                             if (p_t_m_u(a_wtp(t,i),self_z(t,i),co_z(t,i),zt) >0 |p_t_n_p(a_wtp(t,i),self_z(t,i),co_z(t,i),zt)>0) %Hired as a manager firing curent one (me) or as non manager with promotion
%                                 worker_status(t,i)=1; %Uneomployed
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             elseif (p_t_n_u(a_wtp(t,i),self_z(t,i),co_z(t,i),zt)>0) %Hired as a non manager firing curent one
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zt;
%                                 target_wage= (1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zt)+U(co_z(t,i))-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zt),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                             else %Hired as a manager with demotion
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zt;
%                                 target_wage= (1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_wtp(t,i),zt,self_z(t,i))-cost_d+U(co_z(t,i))-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zt,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                             end
%                         % Hire the non manager?    
%                         elseif (h_t_t_nm(at,zt,nt,a_wtp(t,i),self_z(t,i),co_z(t,i))>0) 
%                             %Where does it go?
%                             if (p_t_m_u(a_wtp(t,i),self_z(t,i),co_z(t,i),nt)>0 | p_t_n_p(a_wtp(t,i),self_z(t,i),co_z(t,i),nt)>0) %Hired as a manager firing curent one (me) or as non manager with promotion
%                                 worker_status(t,i)=1; %Uneomployed
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             elseif (p_t_n_u(a_wtp(t,i),self_z(t,i),co_z(t,i),nt)>0) %Hired as a non manager firing curent one
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=nt;
%                                 target_wage= (1-bpw)*(Vth(at,nt,zt)-Vmh(at,zt))+bpw*(Vth(a_wtp(t,i),self_z(t,i),nt)+U(co_z(t,i))-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),nt),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                             else %Hired as a manager with demotion
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=nt;
%                                 target_wage= (1-bpw)*(Vth(at,nt,zt)-Vmh(at,zt))+bpw*(Vth(a_wtp(t,i),nt,self_z(t,i))-cost_d+U(co_z(t,i))-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),nt,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=4; %manager alone
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     end
%                 elseif (r.event(t,i)>del+del+death+prob_matching && r.event(t,i)<=del+del+death+prob_matching+lam) %Worker meets another firm
%                     if (r.state(t,i)<=cdf_firm_dist(1)) %meet vacant firm
%                         av=draw_CDF_1d(cdf_e_edist,r.type(t,i));
%                         %Does get hired?
%                         if (h_e_t_m(a_wtp(t,i),self_z(t,i),co_z(t,i),av)>0)
%                             %Where does it go?
%                             if (p_e_m(av,self_z(t,i))>0) %Hired as a manager
%                                 worker_status(t,i)=2; %manager alone
%                                 a_wtp(t,i)=av;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vmh(av,self_z(t,i))-Veh(av));
%                                 self_wage(t,i)= InterpolateWage(Wmh(:,av,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             else %Hired as a non manager
%                                 worker_status(t,i)=3; %non manager alone
%                                 a_wtp(t,i)=av;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vnh(av,self_z(t,i))-Veh(av));
%                                 self_wage(t,i)= InterpolateWage(Wnh(:,av,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             end
%                         elseif (h_e_t_nm(a_wtp(t,i),self_z(t,i),co_z(t,i),av)>0) %My cowrkers gets hired
%                             %Dont need to know whete it goes just that I lose it
%                             worker_status(t,i)=2; %manager in a team 
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             self_wage(t,i)=self_wage(t-1,i); %Wage does not change
%                         else %NO one gets hired
%                             worker_status(t,i)=4; %manager in a team
%                             %Does the increase? Only if  m was the one tried to get hired
%                             nu_m=[Vmh(av,self_z(t,i)),Vnh(av,self_z(t,i))];
%                             nu_n=[Vnh(av,cot_z(t,i)),Vmh(av,co_z(t,i))];
%                             gain_m=max(nu_m)-Ve(av)+Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i));
%                             gain_n=max(nu_n)-Ve(av)+Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vmh(a_wtp(t,i),self_z(t,i));
%                             target_wage_m=(1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(max(nu_m)-Ve(av));
%                             target_wage_n=(1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vmh(a_wtp(t,i),self_z(t,i)))+bpw*(max(nu_n)-Ve(av));
%                             current_wage_m=InterpolateWage(wgrid,self_wage(t,i),Wtmh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                             current_wage_n=InterpolateWage(wgrid,co_wage(t,i),Wtnh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                             if (gain_m>gain_n && current_wage_m<target_wage_m)
%                                 self_wage(t,i)=InterpolateWage(Wtmh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)),target_wage_m,wgrid);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             elseif (gain_m<gain_n && current_wage_n<target_wage_n)
%                                 co_wage(t,i)=InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)),target_wage_n,wgrid);
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     elseif (r.state(t,i)>cdf_firm_dist(1) && r.state(t,i)<=cdf_firm_dist(2)) %meet firm with manager
%                         [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_m_t_m(a_wtp(t,i),self_z(t,i),co_z(t,i),am,zm)>0)
%                             %Where does it go?
%                             if (p_m_m_u(am,zm,self_z(t,i))>0) %Hired as a manager firing curent one
%                                 worker_status(t,i)=2; %manager alone
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vmh(am,self_z(t,i))+U(zm)-Vmh(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wmh(:,am,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             elseif (p_m_m_d(am,zm,self_z(t,i))>0) %Hired as manager demoting current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(am,self_z(t,i),zm)-cost_d-Vmh(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,am,self_z(t,i),zm),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                                 co_z(t,i)=zm;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_m(k,am,zm);
%                                 co_tenure(t,i)=1;
%                             else %Hired as a non manager
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(am,zm,self_z(t,i))-Vmh(am,zm)); 
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,am,zm,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=zm;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_m(k,am,zm);
%                                 co_tenure(t,i)=1;
%                             end
%                         % My coworker gets hired    
%                         elseif (h_m_t_nm(a_wtp(t,i),self_z(t,i),co_z(t,i),am,zm)>0)
%                             %Dont need to know whete it goes just that I lose it 
%                             worker_status(t,i)=2; %manager alone
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             self_wage(t,i)=self_wage(t-1,i); %Wage does not change
%                         else %NO one gets hired
%                             worker_status(t,i)=4; %manager in a team
%                             %Does the increase? Only if  m was the one tried to get hired
%                             nu_m=[Vmh(am,self_z(t,i))+U(zm),Vth(am,zm,self_z(t,i)),Vth(am,self_z(t,i),zm)-cost_d];
%                             nu_n=[Vmh(am,cot_z(t,i))+U(zm),Vth(am,zm,cot_z(t,i)),Vth(am,cot_z(t,i),zm)-cost_d];
%                             gain_m=max(nu_m)-Vmh(am,zm)-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))+Vnh(a_wtp(t,i),co_z(t,i));
%                             gain_n=max(nu_n)-Vmh(am,zm)-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))+Vmh(a_wtp(t,i),self_z(t,i));
%                             target_wage_m=(1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(max(nu_m)-Vmh(am,zm));
%                             target_wage_n=(1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vmh(a_wtp(t,i),self_z(t,i)))+bpw*(max(nu_n)-Vmh(am,zm));
%                             current_wage_m=InterpolateWage(wgrid,self_wage(t,i),Wtmh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                             current_wage_n=InterpolateWage(wgrid,co_wage(t,i),Wtnh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                             if (gain_m>gain_n && current_wage_m<target_wage_m)
%                                 self_wage(t,i)=InterpolateWage(Wtmh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)),target_wage_m,wgrid);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             elseif (gain_m<gain_n && current_wage_n<target_wage_n)
%                                 co_wage(t,i)=InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)),target_wage_n,wgrid);
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     elseif (r.state(t,i)>cdf_firm_dist(2) && r.state(t,i)<=cdf_firm_dist(3)) %meet firm with non manager
%                         [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_nm_t_m(a_wtp(t,i),self_z(t,i),co_z(t,i),an,zn)>0)
%                             %Where does it go?
%                             if (p_n_n_u(an,zn,self_z(t,i))>0) %Hired as a non manager firing the curent one
%                                 worker_status(t,i)=3; %non manager alone
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vnh(an,self_z(t,i))+U(zn)-Vnh(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wnh(:,an,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             elseif (p_n_n_p(an,zn,self_z(t,i))>0) %Hired as non manager promoting current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(an,zn,self_z(t,i))-cost_p-Vnh(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,an,zn,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                                 co_z(t,i)=zn;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_n(k,an,zn);
%                                 co_tenure(t,i)=1;
%                             else %Hired as a manager
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(an,self_z(t,i),zn)-Vnh(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,an,self_z(t,i),zn),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=zn;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_n(k,an,zn);
%                                 co_tenure(t,i)=1;
%                             end
%                         % My coworker gets hired
%                         elseif (h_nm_t_nm(a_wtp(t,i),self_z(t,i),co_z(t,i),an,zn)>0)
%                             %Does not need to know where it goes just that I lose it
%                             worker_status(t,i)=2; %manager alone
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             self_wage(t,i)=self_wage(t-1,i); %Wage does not change
%                         else %NO one gets hired
%                             worker_status(t,i)=4; %manager in a team
%                             %Does the wage increase? Only if  m was the one tried to get hired
%                             nu_m=[Vnh(an,self_z(t,i))+U(zn),Vth(an,zn,self_z(t,i))-cost_p,Vth(an,self_z(t,i),zn)];
%                             nu_n=[Vnh(an,cot_z(t,i))+U(zn),Vth(an,zn,cot_z(t,i))-cost_p,Vth(an,cot_z(t,i),zn)];
%                             gain_m=max(nu_m)-Vnh(an,zn)-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))+Vnh(a_wtp(t,i),co_z(t,i));
%                             gain_n=max(nu_n)-Vnh(an,zn)-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))+Vmh(a_wtp(t,i),self_z(t,i));
%                             target_wage_m=(1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(max(nu_m)-Vnh(an,zn));
%                             target_wage_n=(1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vmh(a_wtp(t,i),self_z(t,i)))+bpw*(max(nu_n)-Vnh(an,zn));
%                             current_wage_m=InterpolateWage(wgrid,self_wage(t,i),Wtmh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                             current_wage_n=InterpolateWage(wgrid,co_wage(t,i),Wtnh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                             if (gain_m>gain_n && current_wage_m<target_wage_m)
%                                 self_wage(t,i)=InterpolateWage(Wtmh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)),target_wage_m,wgrid);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             elseif (gain_m<gain_n && current_wage_n<target_wage_n)
%                                 co_wage(t,i)=InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)),target_wage_n,wgrid);
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     else %meet firm with team
%                         [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                         %Hire the manager?
%                         if (h_t_t_m(a_wtp(t,i),self_z(t,i),co_z(t,i),at,zt,nt)>0) 
%                             %where does it go?
%                             if (p_t_m_u(at,zt,nt,self_z(t,i))>0) %Hired as a manager firing curent one
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(at,self_z(t,i),nt)+U(zt)-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,at,self_z(t,i),nt),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=nt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tn(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             elseif (p_t_m_d(at,zt,nt,self_z(t,i))>0) %Hired as manager demoting current one
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(at,self_z(t,i),zt)-cost_d+U(nt)-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,at,self_z(t,i),zt),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                                 co_z(t,i)=zt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tm(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             elseif (p_t_n_u(at,zt,nt,self_z(t,i))>0) %Hire as non manager firing current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(at,zt,self_z(t,i))+U(nt)-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,at,zt,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=zt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tm(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             else %Hired as non manager promoting current one
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(at,nt,self_z(t,i))-cost_p+U(zt)-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,at,nt,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                                 co_z(t,i)=nt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tn(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             end
%                         % Hire the non manager?
%                         elseif (h_t_t_nm(a_wtp(t,i),self_z(t,i),co_z(t,i),at,zt,nt)>0) 
%                             %Don't need to know where it goes just that I lose it
%                             worker_status(t,i)=2; %manager alone
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             self_wage(t,i)=self_wage(t-1,i); %Wage does not change
%                         else %NO one gets hired
%                             worker_status(t,i)=4; %manager in a team
%                             %Does the wage increase? Only if  m was the one tried to get hired
%                             nu_m=[Vth(at,self_z(t,i),nt)+U(zt),Vth(at,zt,self_z(t,i))+U(nt),Vth(at,self_z(t,i),zt)-cost_d+U(nt),Vth(at,nz,self_z(t,i))-cost_p+U(zt)];
%                             nu_n=[Vth(at,cot_z(t,i),nt)+U(zt),Vth(at,zt,cot_z(t,i))+U(nt),Vth(at,cot_z(t,i),zt)-cost_d+U(nt),Vth(at,nz,cot_z(t,i))-cost_p+U(zt)];
%                             gain_m=max(nu_m)-Vth(at,zt,nt)-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))+Vnh(a_wtp(t,i),co_z(t,i));
%                             gain_n=max(nu_n)-Vth(at,zt,nt)-Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))+Vmh(a_wtp(t,i),self_z(t,i));
%                             target_wage_m=(1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),co_z(t,i)))+bpw*(max(nu_m)-Vth(at,zt,nt));
%                             target_wage_n=(1-bpw)*(Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vmh(a_wtp(t,i),self_z(t,i)))+bpw*(max(nu_n)-Vth(at,zt,nt));
%                             current_wage_m=InterpolateWage(wgrid,self_wage(t,i),Wtmh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                             current_wage_n=InterpolateWage(wgrid,co_wage(t,i),Wtnh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                             if (gain_m>gain_n && current_wage_m<target_wage_m)
%                                 self_wage(t,i)=InterpolateWage(Wtmh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)),target_wage_m,wgrid);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             elseif (gain_m<gain_n && current_wage_n<target_wage_n)
%                                 co_wage(t,i)=InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)),target_wage_n,wgrid);
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     end
%                 else %Does not meet anyone
%                     worker_status(t,i)=4; %manager in a team
%                     self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                     co_wage(t,i)=co_wage(t-1,i); %does not chenge now
%                     self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                     co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                 end
%             end
%         end
%         % fprintf('SM Status 5')
%         if (worker_status(t-1,i)==5) %Non Manager in a team
%             %Check productivity shock
%             if (r.ashock(t,i)<down_a(a_wtp(t-1,i))) %Firm go down in a
%                 a_wtp(t,i)=a_wtp(t-1,i)-1;
%             elseif (r.ashock(t,i)< down_a(a_wtp(t-1,i))+stay_a(a_wtp(t-1,i)) && r.ashock(t,i)>=down_a(a_wtp(t-1,i))) %Firm stay in a
%                 a_wtp(t,i)=a_wtp(t-1,i);
%             else %Firm go up in a
%                 a_wtp(t,i)=a_wtp(t-1,i)+1;
%             end
%             %Check if learns
%             if (r.qshock(t,i)<down_q(a_wtp(t-1,i),self_z(t-1,i))) %Firm go down in q
%                 self_z(t,i)=self_z(t-1,i)-1;    
%             elseif (r.qshock(t,i)< down_q(a_wtp(t-1,i),self_z(t-1,i))+stay_q(a_wtp(t-1,i),self_z(t-1,i)) && r.qshock(t,i)>=down_q(a_wtp(t-1,i),self_z(t-1,i)))%Firm stay in q
%                 self_z(t,i)=self_z(t-1,i);
%             else %Firm go up in q
%                 self_z(t,i)=self_z(t-1,i)+1;
%             end
%             %Coworker stays the same 
%             co_z(t,i)=co_z(t-1,i);
%             co_wage(t,i)=co_wage(t-1,i);
%             %Check if I die
%             if (r.death(t,i)<death)
%                 worker_status(t,i)=1; %Unemployed        
%                 %new draw but directly on unemployed
%                 self_z(t,i)=draw_CDF_1d(cdf_type_birth,r.type(t,i));
%                 self_wage(t,i)=0.0;
%                 self_age(t,i)=0;
%                 self_tenure(t,i)=0;
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 a_wtp(t,i)=0;
%             else %Dont die 
%                 self_age(t,i)=self_age(t-1,i)+1;
%                 %Events
%                 if (r.event(t,i)<=del) %Leave the firm exogenoulsy
%                     worker_status(t,i)=1; %Unemployed
%                     %new draw but directly on unemployed
%                     self_z(t,i)=draw_CDF_1d(cdf_type_birth,r.type(t,i));
%                     self_wage(t,i)=0.0;
%                     self_tenure(t,i)=0;
%                     co_z(t,i)=0.0;
%                     co_wage(t,i)=0.0;
%                     co_tenure(t,i)=0;
%                     a_wtp(t,i)=0; 
%                 elseif (r.event(t,i)>del && r.event(t,i)<=del+del+death) %Coworker leaves or die
%                     worker_status(t,i)=3; %Non manager alone
%                     self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                     self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                     co_z(t,i)=0.0;
%                     co_wage(t,i)=0.0;
%                     co_tenure(t,i)=0;
%                 elseif (r.event(t,i)>del+del+death && r.event(t,i)<=del+del+death+prob_matching) %Firm meets another worker
%                     if (r.state(t,i)<=cdf_emp_dist(1))%meet unemployed
%                         zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
%                         %Does get hired?
%                         if (h_t_u(zu,a_wtp(t,i),co_z(t,i),self_z(t,i))>0)
%                             %Where does it go?
%                             if (p_t_m_u(a_wtp(t,i),co_z(t,i),self_z(t,i),zu)>0) %Hired as a manager firing current one (my cowroker)
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zu;
%                                 target_wage= (1-bpw)*U(zu)+bpw*(Vth(a_wtp(t,i),zu,self_z(t,i))+U(co_z(t,i))-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zu,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                             elseif (p_t_m_d(a_wtp(t,i),co_z(t,i),self_z(t,i),zu)>0 | p_t_n_u(a_wtp(t,i),co_z(t,i),self_z(t,i),zu)>0) %Hired as a manager with demotion or as non manager with fire
%                                 %In any case I am unemployed
%                                 worker_status(t,i)=1; %Unemployed
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             else %Hired as a non manager with promotion of the current one (me)
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zu;
%                                 target_wage= (1-bpw)*U(zu)+bpw*(Vth(a_wtp(t,i),self_z(t,i),zu)+U(co_z(t,i))-cost_p-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zu),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=5; %non manager in a team
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                             co_wage(t,i)=co_wage(t-1,i); %does not chenge now
%                         end
%                     elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<cdf_emp_dist(2)) %meet employed manager
%                         [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_t_m(am,zm,a_wtp(t,i),co_z(t,i),self_z(t,i))>0)
%                             %Where does it go?
%                             if (p_t_m_u(a_wtp(t,i),co_z(t,i),self_z(t,i),zm)>0) %Hired as a manager firing current one (my cowroker)
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zm;
%                                 target_wage= (1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_wtp(t,i),zm,self_z(t,i))+U(co_z(t,i))-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),self_z(t,i),zm),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                             elseif (p_t_m_d(a_wtp(t,i),co_z(t,i),self_z(t,i),zm)>0 | p_t_n_u(a_wtp(t,i),co_z(t,i),self_z(t,i),zm)>0) %Hired as a manager with demotion or as non manager with fire
%                                 %In any case I am unemployed
%                                 worker_status(t,i)=1; %Unemployed
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             else %Hired as a non manager with promotion of the current one (me)
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zm;
%                                 target_wage= (1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zm)+U(co_z(t,i))-cost_p-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zm),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=5; %non manager in a team
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     elseif (r.state(t,i)>=cdf_emp_dist(2) && r.state(t,i)<cdf_emp_dist(3)) %meet employed non manager
%                         [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                         %Does get hired?
%                         if (h_t_nm(an,zn,a_wtp(t,i),co_z(t,i),self_z(t,i))>0)
%                             %Where does it go?
%                             if (p_t_m_u(a_wtp(t,i),co_z(t,i),self_z(t,i),zn)>0) %Hired as a manager firing current one (my cowroker)
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zn;
%                                 target_wage= (1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_wtp(t,i),zn,self_z(t,i))+U(co_z(t,i))-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zn,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                             elseif (p_t_m_d(a_wtp(t,i),co_z(t,i),self_z(t,i),zn)>0 | p_t_n_u(a_wtp(t,i),co_z(t,i),self_z(t,i),zn)>0) %Hired as a manager with demotion or as non manager with fire
%                                 %In any case I am unemployed
%                                 worker_status(t,i)=1; %Unemployed
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             else %Hired as a non manager with promotion of the current one (me)
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zn;
%                                 target_wage= (1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zn)+U(co_z(t,i))-cost_p-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zn),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                             end
%                         else %Does not get hired
%                             worker_status(t,i)=5; %non manager in a team
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     else %meet employed team
%                         [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                         %Hire the manager?
%                         if (h_t_t_m(a_wtp(t,i),co_z(t,i),self_z(t,i),at,zt,nt)>0) 
%                             %where does it go?
%                             if (p_t_m_u(a_wtp(t,i),co_z(t,i),self_z(t,i),zt)>0) %Hired as a manager firing current one (my cowroker)
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zt;
%                                 target_wage= (1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_wtp(t,i),zt,self_z(t,i))+U(co_z(t,i))-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),zt,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                             elseif (p_t_m_d(a_wtp(t,i),co_z(t,i),self_z(t,i),zt)>0 | p_t_n_u(a_wtp(t,i),co_z(t,i),self_z(t,i),zt)>0) %Hired as a manager with demotion or as non manager with fire
%                                 %In any case I am unemployed
%                                 worker_status(t,i)=1; %Unemployed
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             else %Hired as a non manager with promotion of the current one (me) 
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=zt;
%                                 target_wage= (1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_wtp(t,i),self_z(t,i),zt)+U(co_z(t,i))-cost_p-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),zt),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                             end
%                         % Hire the non manager?
%                         elseif (h_t_t_nm(a_wtp(t,i),co_z(t,i),self_z(t,i),at,zt,nt)>0) 
%                             %Where does it go?
%                             if (p_t_m_u(a_wtp(t,i),co_z(t,i),self_z(t,i),nt)>0) %Hired as a manager firing current one (my cowroker)
%                                 worker_status(t,i)=5; %non manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=nt;
%                                 target_wage= (1-bpw)*(Vth(at,znt,nt)-Vmh(at,nt))+bpw*(Vth(a_wtp(t,i),nt,self_z(t,i))+U(co_z(t,i))-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtmh(:,a_wtp(t,i),nt,self_z(t,i)),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                             elseif (p_t_m_d(a_wtp(t,i),co_z(t,i),self_z(t,i),nt)>0 | p_t_n_u(a_wtp(t,i),co_z(t,i),self_z(t,i),nt)>0) %Hired as a manager with demotion or as non manager with fire
%                                 %In any case I am unemployed
%                                 worker_status(t,i)=1; %Unemployed
%                                 self_wage(t,i)=0.0;
%                                 self_tenure(t,i)=0;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                                 a_wtp(t,i)=0;
%                                 fire_m(t,i)=1;
%                             else %Hired as a non manager with promotion of the current one (me)
%                                 worker_status(t,i)=4; %manager in a team
%                                 self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                                 co_z(t,i)=nt;
%                                 target_wage= (1-bpw)*(Vth(at,nt,zt)-Vmh(at,zt))+bpw*(Vth(a_wtp(t,i),self_z(t,i),nt)+U(co_z(t,i))-cost_p-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                                 co_wage(t,i)= InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),nt),target_wage,wgrid);
%                                 co_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                             end
%                         else %Does not hire anyone
%                             worker_status(t,i)=5; %non manager in a team
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     end
%                 elseif (r.event(t,i)>del+del+death+prob_matching && r.event(t,i)<=del+del+death+prob_matching+lam) %Worker meets another firm
%                     if (r.state(t,i)<=cdf_firm_dist(1)) %meet Vacant firm
%                         av=draw_CDF_1d(cdf_e_edist,r.type(t,i));
%                         %Does Manager get hired? (My coworker)
%                         if (h_e_t_m(a_wtp(t,i),co_z(t,i),self_z(t,i),av)>0)
%                             worker_status(t,i)=3; %Non manager alone
%                             %Does not matter where it goes all that matters is that I lose it
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                         %Does non manager get hired? (Me)
%                         elseif (h_e_t_nm(a_wtp(t,i),co_z(t,i),self_z(t,i),av)>0)
%                             %Where does it go?
%                             if (p_e_m(av,self_z(t,i))>0) %Hire as manager
%                                 worker_status(t,i)=2; %manager alone
%                                 a_wtp(t,i)=av;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vmh(av,self_z(t,i))-Veh(av));
%                                 self_wage(t,i)= InterpolateWage(Wmh(:,av,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             else %Hire as non manager
%                                 worker_status(t,i)=3;
%                                 a_wtp(t,i)=av;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vnh(av,self_z(t,i))-Veh(av));
%                                 self_wage(t,i)= InterpolateWage(Wnh(:,av,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             end
%                         else %No one gets hired
%                             worker_status(t,i)=5; %non manager in a team
%                             %Does the wage increase? Only if  nm was the one tried to get hired
%                             nu_m=[Vmh(av,co_z(t,i)),Vnh(av,co_z(t,i))]
%                             nu_n=[Vmh(av,self_z(t,i)),Vnh(av,self_z(t,i))]
%                             gain_m=max(nu_m)-Veh(av)-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))+Vnh(a_wtp(t,i),self_z(t,i));
%                             gain_n=max(nu_n)-Veh(av)-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))+Vmh(a_wtp(t,i),co_z(t,i));
%                             target_wage_m=(1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vnh(a_wtp(t,i),self_z(t,i)))+bpw*(max(nu_m)-Veh(av));
%                             target_wage_n=(1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(max(nu_n)-Veh(av));
%                             current_wage_m=InterpolateWage(wgrid,co_wage(t,i),Wtmh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                             current_wage_n=InterpolateWage(wgrid,self_wage(t,i),Wtnh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                             if (gain_m>gain_n && current_wage_m<target_wage_m)
%                                 co_wage(t,i)=InterpolateWage(Wtmh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)),target_wage_m,wgrid);
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             elseif (gain_m<gain_n && current_wage_n<target_wage_n)
%                                 self_wage(t,i)=InterpolateWage(Wtnh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)),target_wage_n,wgrid);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     elseif (r.state(t,i)>cdf_firm_dist(1) && r.state(t,i)<=cdf_firm_dist(2)) %meet firm with manager alone
%                         [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
%                         %Does it hire the manager? (My coworker)
%                         if (h_m_t_m(a_wtp(t,i),co_z(t,i),self_z(t,i),am,zm)>0)
%                             %Does not matter where it goes all that matters is that I lose it
%                             worker_status(t,i)=3; %Non manager alone
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                         %Does it hire the non manager? (Me)    
%                         elseif (h_m_t_nm(a_wtp(t,i),co_z(t,i),self_z(t,i),am,zm)>0)
%                             %Where does it go?
%                             if (p_m_m_u(am,zm,self_z(t,i))>0) %Hire as manager firing current 
%                                 worker_status(t,i)=2; %manager alone
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vmh(am,self_z(t,i))-Vmh(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wmh(:,am,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             elseif (p_m_m_d(am,zm,self_z(t,i))>0) %Hire as manager demoting current
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(am,self_z(t,i),zm)-cost_d-Vmh(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,am,self_z(t,i),zm),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                                 co_z(t,i)=zm;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_m(k,am,zm);
%                                 co_tenure(t,i)=1;
%                             else %Hire as non manager 
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=am;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(am,zm,self_z(t,i))-Vmh(am,zm));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,am,zm,self_z(t,i)),target_wage,wgrid);   
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=zm;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_m(k,am,zm);
%                                 co_tenure(t,i)=1;
%                             end
%                         else %NO one gets hired
%                             worker_status(t,i)=5; %non manager in a team
%                             %Does the wage increase? Only if  nm was the one tried to get hired
%                             nu_m=[Vmh(am,cot_z(t,i))+U(zm),Vth(am,zm,co_z(t,i)),Vth(am,co_z(t,i),zm)-cost_d]
%                             nu_n=[Vmh(am,self_z(t,i))+U(zm),Vth(am,zm,co_z(t,i)),Vth(am,co_z(t,i),zm)-cost_d]
%                             gain_m=max(nu_m)-Vmh(am,zm)-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))+Vnh(a_wtp(t,i),self_z(t,i));
%                             gain_n=max(nu_n)-Vmh(am,zm)-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))+Vmh(a_wtp(t,i),co_z(t,i));
%                             target_wage_m=(1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vnh(a_wtp(t,i),self_z(t,i)))+bpw*(max(nu_m)-Vmh(am,zm));
%                             target_wage_n=(1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(max(nu_n)-Vmh(am,zm));
%                             current_wage_m=InterpolateWage(wgrid,co_wage(t,i),Wtmh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                             current_wage_n=InterpolateWage(wgrid,self_wage(t,i),Wtnh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                             if (gain_m>gain_n && current_wage_m<target_wage_m)
%                                 co_wage(t,i)=InterpolateWage(Wtmh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)),target_wage_m,wgrid);
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             elseif (gain_m<gain_n && current_wage_n<target_wage_n)
%                                 self_wage(t,i)=InterpolateWage(Wtnh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)),target_wage_n,wgrid);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     elseif (r.state(t,i)>cdf_firm_dist(2) && r.state(t,i)<=cdf_firm_dist(3)) %meet firm with non manager alone
%                         [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
%                         %Does it hire the manager? (My coworker)
%                         if (h_nm_t_m(a_wtp(t,i),co_z(t,i),self_z(t,i),an,zn)>0)
%                             %Does not matter where it goes all that matters is that I lose it
%                             worker_status(t,i)=3; %Non manager alone
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                         %Does it hire the non manager? (Me) 
%                         elseif (h_nm_t_nm(a_wtp(t,i),co_z(t,i),self_z(t,i),an,zn)>0)
%                             %Where does it go?
%                             if (p_n_n_u(an,zn,self_z(t,i))>0) %Hire as non manager firing current
%                                 worker_status(t,i)=3; %non manager alone
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vnh(an,self_z(t,i))+U(zn)-Vnh(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wnh(:,an,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=0.0;
%                                 co_wage(t,i)=0.0;
%                                 co_tenure(t,i)=0;
%                             elseif (p_n_n_p(an,zn,self_z(t,i))>0) %Hire as non manager promoting current
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(an,zn,self_z(t,i))-cost_p-Vnh(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,an,zn,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                                 co_z(t,i)=zn;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_n(k,an,zn);
%                                 co_tenure(t,i)=1;
%                             else %Hire as manager
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=an;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(an,self_z(t,i),zn)-Vnh(an,zn));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,an,self_z(t,i),zn),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=zn;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_n(k,an,zn);
%                                 co_tenure(t,i)=1;
%                             end
%                         else %NO one gets hired
%                             worker_status(t,i)=5; %non manager in a team
%                             %Does the wage increase? Only if  nm was the one tried to get hired
%                             nu_m=[Vnh(an,cot_z(t,i))+U(zn),Vth(an,zn,co_z(t,i))-cost_p,Vth(an,co_z(t,i),zn)]
%                             nu_n=[Vnh(an,self_z(t,i))+U(zn),Vth(an,zn,self_z(t,i))-cost_p,Vth(an,self_z(t,i),zn)]
%                             gain_m=max(nu_m)-Vnh(an,zn)-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))+Vnh(a_wtp(t,i),self_z(t,i));
%                             gain_n=max(nu_n)-Vnh(an,zn)-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))+Vmh(a_wtp(t,i),co_z(t,i));
%                             target_wage_m=(1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vnh(a_wtp(t,i),self_z(t,i)))+bpw*(max(nu_m)-Vnh(an,zn));
%                             target_wage_n=(1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(max(nu_n)-Vnh(an,zn));
%                             current_wage_m=InterpolateWage(wgrid,co_wage(t,i),Wtmh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                             current_wage_n=InterpolateWage(wgrid,self_wage(t,i),Wtnh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                             if (gain_m>gain_n && current_wage_m<target_wage_m)
%                                 co_wage(t,i)=InterpolateWage(Wtmh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)),target_wage_m,wgrid);
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             elseif (gain_m<gain_n && current_wage_n<target_wage_n)
%                                 self_wage(t,i)=InterpolateWage(Wtnh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)),target_wage_n,wgrid);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     else %Meet team
%                         [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
%                         %Does it hire the manager? (My coworker)
%                         if (h_t_t_m(a_wtp(t,i),co_z(t,i),self_z(t,i),at,zt,nt)>0)
%                             %Does not matter where it goes all that matters is that I lose it
%                             worker_status(t,i)=3;
%                             co_z(t,i)=0.0;
%                             co_wage(t,i)=0.0;
%                             co_tenure(t,i)=0;
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                         %Does it hire the non manager? (Me) 
%                         elseif (h_t_t_nm(a_wtp(t,i),co_z(t,i),self_z(t,i),at,zt,nt)>0)
%                             %Where does it go?
%                             if (p_t_n_u(at,zt,nt,self_z(t,i))>0) %Hire as non manager firing current
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(at,zt,self_z(t,i))+U(nt)-Vth(at,zt,nz));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,at,zt,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 co_z(t,i)=zt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tm(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             elseif (p_t_n_p(at,zt,nt,self_z(t,i))>0) %Hire as non manager promoting current
%                                 worker_status(t,i)=5; %non manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(at,nt,self_z(t,i))-cost_p+U(zt)-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtnh(:,at,nt,self_z(t,i)),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_n(t,i)=1;
%                                 prom_sm(t,i)=1;
%                                 co_z(t,i)=nt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tn(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             elseif  (p_t_m_u(at,zt,nt,self_z(t,i))>0) %Hire as manager firing current
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(at,self_z(t,i),nt)+U(zt)-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,at,self_z(t,i),nt),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 co_z(t,i)=nt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tn(k,at,zt,nt);
%                                 co_tenure(t,i)=1;
%                             else %Hire as manager demoting current
%                                 worker_status(t,i)=4; %manager in a team
%                                 a_wtp(t,i)=at;
%                                 target_wage= (1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(Vth(at,self_z(t,i),zt)-cost_d+U(nt)-Vth(at,zt,nt));
%                                 self_wage(t,i)= InterpolateWage(Wtmh(:,at,self_z(t,i),zt),target_wage,wgrid);
%                                 self_tenure(t,i)=1;
%                                 hire_m(t,i)=1;
%                                 prom_sm(t,i)=-1;
%                                 co_z(t,i)=zt;
%                                 k=1+floor(r.tie_break(t,i)*nquantiles);
%                                 co_wage(t,i)=quantiles_tm(k,at,zt,nt); 
%                                 co_tenure(t,i)=1;
%                             end
%                         else %NO one gets hired
%                             worker_status(t,i)=5; %non manager in a team
%                             %Does the wage increase? Only if  nm was the one tried to get hired
%                             nu_m=[Vth(at,zt,co_z(t,i))+U(nt),Vth(at,nt,co_z(t,i))-cost_p+U(zt),Vth(at,co_z(t,i),nt)+U(zt),Vth(at,co_z(t,i),zt)-cost_d+U(nt)]
%                             nu_n=[Vth(at,zt,self_z(t,i))+U(nt),Vth(at,nt,self_z(t,i))-cost_p+U(zt),Vth(at,self_z(t,i),nt)+U(zt),Vth(at,self_z(t,i),zt)-cost_d+U(nt)]
%                             gain_m=max(nu_m)-Vth(at,zt,nt)-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))+Vnh(a_wtp(t,i),self_z(t,i));
%                             gain_n=max(nu_n)-Vth(at,zt,nt)-Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))+Vmh(a_wtp(t,i),co_z(t,i));
%                             target_wage_m=(1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vnh(a_wtp(t,i),self_z(t,i)))+bpw*(max(nu_m)-Vth(at,zt,nt));
%                             target_wage_n=(1-bpw)*(Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i)))+bpw*(max(nu_n)-Vth(at,zt,nt));
%                             current_wage_m=InterpolateWage(wgrid,co_wage(t,i),Wtmh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                             current_wage_n=InterpolateWage(wgrid,self_wage(t,i),Wtnh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                             if (gain_m>gain_n && current_wage_m<target_wage_m)
%                                 co_wage(t,i)=InterpolateWage(Wtmh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)),target_wage_m,wgrid);
%                                 self_wage(t,i)=self_wage(t-1,i);
%                             elseif (gain_m<gain_n && current_wage_n<target_wage_n)
%                                 self_wage(t,i)=InterpolateWage(Wtnh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)),target_wage_n,wgrid);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             else
%                                 self_wage(t,i)=self_wage(t-1,i);
%                                 co_wage(t,i)=co_wage(t-1,i);
%                             end
%                             self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                             co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                         end
%                     end
%                 else %Does not meet any firm
%                     worker_status(t,i)=5; %non manager in a team
%                     self_wage(t,i)=self_wage(t-1,i); %does not chenge now
%                     self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                     co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                     co_wage(t,i)=co_wage(t-1,i); %does not chenge now
%                 end
%             end
%         end



%         %% Reallocation stage for workers
%         % fprintf('Reallocation 2')
%         if (worker_status(t,i)==2) %Manager alone
%             %Gets fired?
%             if (d_m(a_wtp(t,i),self_z(t,i))>0)
%                 worker_status(t,i)=1; %Unemployed
%                 self_wage(t,i)=0.0;
%                 self_tenure(t,i)=0;
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 a_wtp(t,i)=0;
%                 fire_m(t,i)=1;
%             %Gets reallocated?
%             elseif (r_m(a_wtp(t,i),self_z(t,i))>0)
%                 worker_status(t,i)=3; %non manager alone
%                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 prom_realloc(t,i)=-1;
%             else %Does not change
%                 %Does the wage decrease?
%                 current_wage=InterpolateWage(wgrid,self_wage(t,i),Wmh(:,a_wtp(t,i),self_z(t,i)));
%                 target_wage= Vm(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i));
%                 if (current_wage>target_wage)
%                     self_wage(t,i)=InterpolateWage(Wmh(:,a_wtp(t,i),self_z(t,i)),target_wage,wgrid);
%                 end
%             end
%         end
%         % fprintf('Reallocation 3')
%         if (worker_status(t,i)==3) %Non manager alone
%             %Gets fired?
%             if (d_n(a_wtp(t,i),self_z(t,i))>0)
%                 worker_status(t,i)=1; %Unemployed
%                 self_wage(t,i)=0.0;
%                 self_tenure(t,i)=0;
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 a_wtp(t,i)=0;
%                 fire_n(t,i)=1;
%             %Gets reallocated?
%             elseif (r_n(a_wtp(t,i),self_z(t,i))>0)
%                 worker_status(t,i)=2; %manager alone
%                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 prom_realloc(t,i)=1;
%             else %Does not change
%                 %Does the wage decrease?
%                 current_wage=InterpolateWage(wgrid,self_wage(t,i),Wnh(:,a_wtp(t,i),self_z(t,i)));
%                 target_wage= Vn(a_wtp(t,i),self_z(t,i))-Veh(a_wtp(t,i));
%                 if (current_wage>target_wage)
%                     self_wage(t,i)=InterpolateWage(Wnh(:,a_wtp(t,i),self_z(t,i)),target_wage,wgrid);
%                 end
%             end
%         end
%         % fprintf('Reallocation 4')
%         if (worker_status(t,i)==4) %Manager in a team
%             if (d_t_m(a_wtp(t,i),self_z(t,i),co_z(t,i))>0 || d_t_b(a_wtp(t,i),self_z(t,i),co_z(t,i))>0) %I get fired
%                 worker_status(t,i)=1; %Unemployed
%                 self_wage(t,i)=0.0;
%                 self_tenure(t,i)=0;
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 a_wtp(t,i)=0;
%                 fire_m(t,i)=1;
%             elseif (d_t_n(a_wtp(t,i),self_z(t,i),co_z(t,i))>0 && r_t_m(a_wtp(t,i),self_z(t,i),co_z(t,i))>0) %Fire nm and reallocate me
%                 worker_status(t,i)=3; %non manager alone
%                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 prom_realloc(t,i)=-1;
%                 fire_n(t,i)=1;
%             elseif (d_t_n(a_wtp(t,i),self_z(t,i),co_z(t,i))>0 && r_t_m(a_wtp(t,i),self_z(t,i),co_z(t,i))<=0) %Fire nm and do not reallocate me
%                 worker_status(t,i)=2; %manager alone
%                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 prom_realloc(t,i)=1;
%                 fire_n(t,i)=1;
%             elseif (r_t_m(a_wtp(t,i),self_z(t,i),co_z(t,i))>0 && r_t_n(a_wtp(t,i),self_z(t,i),co_z(t,i))>0) %reallocate both
%                 worker_status(t,i)=5; %non manager in a team
%                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                 co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                 prom_realloc(t,i)=1;
%             else %Does not change
%                 %Does the wage decrease?
%                 current_wage_m=InterpolateWage(wgrid,self_wage(t,i),Wtmh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                 current_wage_n=InterpolateWage(wgrid,co_wage(t,i),Wtnh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)));
%                 target_wage_m=Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vnh(a_wtp(t,i),self_z(t,i));
%                 target_wage_n=Vth(a_wtp(t,i),self_z(t,i),co_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i));
%                 if (current_wage_m>target_wage_m)
%                     self_wage(t,i)=InterpolateWage(Wtmh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)),target_wage_m,wgrid);
%                 end
%                 if (current_wage_n>target_wage_n)
%                     co_wage(t,i)=InterpolateWage(Wtnh(:,a_wtp(t,i),self_z(t,i),co_z(t,i)),target_wage_n,wgrid);
%                 end
%             end
%         end
%         % fprintf('Reallocation 5')
%         if (worker_status(t,i)==5) %Non manager in a team
%             if (d_t_n(a_wtp(t,i),co_z(t,i),self_z(t,i))>0 || d_t_b(a_wtp(t,i),co_z(t,i),self_z(t,i))>0) %I get fired
%                 worker_status(t,i)=1; %Unemployed
%                 self_wage(t,i)=0.0;
%                 self_tenure(t,i)=0;
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 a_wtp(t,i)=0;
%                 fire_m(t,i)=1;
%             elseif (d_t_m(a_wtp(t,i),co_z(t,i),self_z(t,i))>0 && r_t_n(a_wtp(t,i),co_z(t,i),self_z(t,i))>0) %Fire m and realloc me
%                 worker_status(t,i)=2; %manager alone
%                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 fire_m(t,i)=1;
%                 prom_realloc(t,i)=1;
%             elseif (d_t_m(a_wtp(t,i),co_z(t,i),self_z(t,i))>0 && r_t_n(a_wtp(t,i),co_z(t,i),self_z(t,i))<=0) %Fire m and do not realloc me
%                 worker_status(t,i)=3; %non manager alone
%                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                 co_z(t,i)=0.0;
%                 co_wage(t,i)=0.0;
%                 co_tenure(t,i)=0;
%                 fire_m(t,i)=1;
%             elseif (r_t_n(a_wtp(t,i),co_z(t,i),self_z(t,i))>0 && r_t_m(a_wtp(t,i),co_z(t,i),self_z(t,i))>0) %reallocate both
%                 worker_status(t,i)=4; %manager in a team
%                 self_tenure(t,i)=self_tenure(t-1,i)+1; %One more month in the firm
%                 co_tenure(t,i)=co_tenure(t-1,i)+1; %One more month with the same co worker
%                 prom_realloc(t,i)=1; % At this stage promotion means both reallocated
%             else %Does not change
%                 %Does the wage decrease?
%                 current_wage_m=InterpolateWage(wgrid,co_wage(t,i),Wtmh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                 current_wage_n=InterpolateWage(wgrid,self_wage(t,i),Wtnh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)));
%                 target_wage_m=Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vnh(a_wtp(t,i),self_z(t,i));
%                 target_wage_n=Vth(a_wtp(t,i),co_z(t,i),self_z(t,i))-Vmh(a_wtp(t,i),co_z(t,i));
%                 if (current_wage_m>target_wage_m)
%                     co_wage(t,i)=InterpolateWage(Wtmh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)),target_wage_m,wgrid);
%                 end
%                 if (current_wage_n>target_wage_n)
%                     self_wage(t,i)=InterpolateWage(Wtnh(:,a_wtp(t,i),co_z(t,i),self_z(t,i)),target_wage_n,wgrid);
%                 end
%             end
%         end
%     end
% end
% toc


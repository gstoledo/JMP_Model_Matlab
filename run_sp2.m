%% Specification 2
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
aup=0.2;
adown=0.1;
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
        fprintf('Sp3 Joint Converged in %d iterations\n',it_joint)
    end

end
%End timers
toc
    
%Check sum of the eplus final distributions
nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");

save('model_solutions_sp2'+location,'Ve','Vm','Vn','Vt','U','Veh','Vmh','Vnh','Vth','eplus_udist','eplus_edist','eplus_mdist','eplus_ndist','eplus_tdist','nplus','popplus');


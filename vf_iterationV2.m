%Function to perform the value function iteration
%Correcting for the u issue
function [Ve, Vm, Vn, Vt, U, Veh, Vmh, Vnh, Vth, Vetl, Vmtl, Vntl,Vttl,Utl] = vf_iterationV2(e_edist,e_mdist,e_ndist,e_tdist,e_udist,Veini,Vmini,Vnini,Vtini,Uini,ats,tpts,cost_d,cost_p,tr,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed)
    %Use initial guesses for value functions
    Veup = Veini;
    Vmup = Vmini;
    Vnup = Vnini;
    Vtup = Vtini;
    Uup = Uini;
    
    %Mass of U
    %u=sum(e_udist); %mass of unemployed
    
    %Iteration parameters
    diff=100;
    diffmax=1e-6;
     
    it=0; 
    itmax=100000;
    
    %speed=1; %Speed of convergence
     
    while diff>diffmax & it<itmax
     it=it+1;
     if it==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Update stage
        Ve=Veup;
        Vm=Vmup;
        Vn=Vnup;
        Vt=Vtup;
        U=Uup;
     end
    
     if it>1
        Ve=(1-speed)*Ve + speed*Veup;
        Vm= (1-speed)*Vm + speed*Vmup;
        Vn= (1-speed)*Vn + speed*Vnup;
        Vt= (1-speed)*Vt + speed*Vtup;
        U= (1-speed)*U + speed*Uup;
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Fire/Promotion VFs
        %Empty firm
        Veh=Ve;
        
        %Firm with manager (a,z)
        Vm_alloc(:,:,1)=Vm; %Keep as it is
        Vm_alloc(:,:,2)=Vn-cost_d; %Demote z
        Vm_alloc(:,:,3)=repmat(Ve',1, tpts)+repmat(U,ats,1);%Fire z
        Vmh=max(Vm_alloc,[],3);
        
        
        %Firm with non-manager (a,q) 
        Vn_alloc(:,:,1)=Vn; %Keep as it is
        Vn_alloc(:,:,2)=Vm-cost_p; %Promote q
        Vn_alloc(:,:,3)=repmat(Ve',1, tpts)+repmat(U,ats,1); %Fire q
        Vnh=max(Vn_alloc,[],3);
        
        %Firm with team (a,z,q)
        Vt_alloc(:,:,:,1)=Vt; %Keep as it is
        Vt_alloc(:,:,:,2)=permute(repmat(Vm,1,1,tpts),[1 3 2])+repmat(U,ats,1,tpts)-cost_p; %Fire z and promote q
        Vt_alloc(:,:,:,3)=repmat(Vn,1,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]) - cost_d; %Fire q and demote z
        Vt_alloc(:,:,:,4)=repmat(Vm,1,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]); %Fire q 
        Vt_alloc(:,:,:,5)=permute(repmat(Vn,1,1,tpts),[1 3 2])+repmat(U,ats,1,tpts); %Fire z
        Vt_alloc(:,:,:,6)=permute(Vt, [1 3 2])-cost_d-cost_p; %Swap z and q with promotion and demotion costs
        Vt_alloc(:,:,:,7)=permute(repmat(Ve',1,tpts,tpts),[1 3 2])+repmat(U,ats,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]); %Fire z and q
        Vth=max(Vt_alloc,[],4);
        
        %Unemployed
        Uh=U;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%Empty firm continuation value Ve(a)
        %Gains from trade 
        [gt_ef_meet_u, gt_ef_meet_m, gt_ef_meet_nm, gt_ef_meet_t_m, gt_ef_meet_t_nm, gt_ef_meet_t]=ef_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,tr);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Firm with manager continuation value Vm(a,z)
        %Gains from trade
        [gt_em_meet_u, gt_em_meet_m, gt_em_meet_nm, gt_em_meet_t_m, gt_em_meet_t_nm, gt_em_meet_t]=mf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,tr,cost_d);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Firm with no manager continuation value Vn(a,q)
        %Gains from trade
        [gt_en_meet_u, gt_en_meet_m, gt_en_meet_nm, gt_en_meet_t_m, gt_en_meet_t_nm, gt_en_meet_t]=nf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,tr,cost_p);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Firm with team continuation value Vt(a,z,q)
        %Gains from trade
        [gt_tf_meet_u, gt_tf_meet_m, gt_tf_meet_nm, gt_tf_meet_t_m, gt_tf_meet_t_nm, gt_tf_meet_t]=tf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,tr,cost_p,cost_d);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Expected Gains of poaching 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Expected gains from trade for empty firm
        %Empty firm (a) poaching from unemplyment
        Egt_e_u=e_udist*gt_ef_meet_u; %From unemployment to empty firm
        
        %Empty firm (a) poaching from firm with manager
        Egt_e_m=reshape(sum(e_mdist.*gt_ef_meet_m, [1 2]),[1 ats]); %From manager to empty firm
        
        %Empty firm (a) poaching from firm with non manager
        Egt_e_n=reshape(sum(e_ndist.*gt_ef_meet_nm, [1 2]), [1 ats]);
        
        %Empty firm (a) poaching from firm with team
        Egt_e_t=reshape(sum(e_tdist.*gt_ef_meet_t, [1 2 3]), [1 ats]);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Expected Gains of poaching for firm with manager
        %Firm (a,z) Poaching from unemployment
        Egt_m_u=reshape(permute(sum(e_udist.*permute(gt_em_meet_u, [2 1 3]),2),[2 1 3]),ats,tpts); %From unemployment to firm with manager
        
        %Firm (a,z) Poaching from firm with manager
        Egt_m_m= reshape(sum(e_mdist.*gt_em_meet_m, [1 2]), ats, tpts); %From firm with manager manager to firm with manager
        
        %Firm (a,z) Poaching from firm with no manager
        Egt_m_n=reshape(sum(e_ndist.*gt_em_meet_nm, [1 2]), ats,tpts); %From firm with no manager to firm with manager
        
        %Firm (a,z) Poaching from firm with team
        Egt_m_t=reshape(sum(e_tdist.*gt_em_meet_t, [1 2 3]), ats, tpts); %From firm with team to firm with manager
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Expected Gains of poaching for firm with no manager
        %Firm (a,q) Poaching from unemployment
        Egt_n_u=reshape(permute(sum(e_udist.*permute(gt_en_meet_u, [2 1 3]),2),[2 1 3]),ats,tpts); %From unemployment to firm with no manager
        
        %Firm (a,q) Poaching from firm with manager
        Egt_n_m= reshape(sum(e_mdist.*gt_en_meet_m, [1 2]), ats, tpts); %From firm with manager to firm with no manager
        
        %Firm (a,q) Poaching from firm with no manager
        Egt_n_n=reshape(sum(e_ndist.*gt_en_meet_nm, [1 2]), ats,tpts); %From firm with no manager to firm with no manager
        
        %Firm (a,q) Poaching from firm with team
        Egt_n_t=reshape(sum(e_tdist.*gt_en_meet_t, [1 2 3]), ats, tpts); %From firm with team to firm with no manager
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Expected Gains of poaching for firm with team
        %Firm (a,z,q) Poaching from unemployment
        Egt_t_u=reshape(permute(sum(e_udist.*permute(gt_tf_meet_u, [2 1 3 4]),2),[2 1 3 4]),ats,tpts,tpts); %From unemployment to firm with team
        
        %Firm (a,z,q) Poaching from firm with manager
        Egt_t_m= reshape(sum(e_mdist.*gt_tf_meet_m, [1 2]), ats, tpts, tpts); %From firm with manager to firm with team
        
        %Firm (a,z,q) Poaching from firm with no manager
        Egt_t_n=reshape(sum(e_ndist.*gt_tf_meet_nm, [1 2]), ats,tpts, tpts); %From firm with no manager to firm with team
        
        %Firm (a,z,q) Poaching from firm with team
        Egt_t_t=reshape(sum(e_tdist.*gt_tf_meet_t, [1 2 3]), ats, tpts, tpts); %From firm with team to firm with team
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Expected Gains of being poached
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Expected gains of being poached for unemployed worker 
        %Unemp z Poached from empty firm
        Epgt_u_e=e_edist*permute(gt_ef_meet_u, [2 1]); %From unemployment to empty firm
        
        %Unemp z Poached from firm with manager
        Epgt_u_m=reshape(sum(e_mdist.*permute(gt_em_meet_u, [2 3 1]),[1 2]),1,tpts); %From unemployment to firm with manager
        
        %Unemp z Poached from firm with no manager
        Epgt_u_n=reshape(sum(e_ndist.*permute(gt_en_meet_u, [2 3 1]), [1 2]),1,tpts); %From unemployment to firm with no manager
        
        %Unemp z Poached from firm with team
        Epgt_u_t=reshape(sum(e_tdist.*permute(gt_tf_meet_u, [2 3 4 1]),[1 2 3]),1,tpts); %From unemployment to firm with team
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Expected Gains of being poached for firm with manager
        %Firm (a,z) Poached from empty firm
        Epgt_m_e=sum(e_edist'.*(gt_ef_meet_m),3); %From manager to empty firm 
        
        %Firm (a,z) Poached from firm with manager
        Epgt_m_m=reshape(sum(e_mdist.*permute(gt_em_meet_m, [3 4 1 2]),[1 2]),ats,tpts); %From manager to manager
        
        %Firm (a,z) Poached from firm with no manager
        Epgt_m_n=reshape(sum(e_ndist.*permute(gt_en_meet_m, [3 4 1 2]), [1 2]),ats,tpts); %From manager to no manager
        
        %Firm (a,z) Poached from firm with team
        Epgt_m_t=reshape(sum(e_tdist.*permute(gt_tf_meet_m, [3 4 5 1 2]),[1 2 3]),ats,tpts); %From manager to team
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Expected Gains of being poached for firm with no manager
        %Firm (a,q) Poached from empty firm
        Epgt_n_e=sum(e_edist'.*(gt_ef_meet_nm),3); %From no manager to empty firm
        
        %Firm (a,q) Poached from firm with manager
        Epgt_n_m=reshape(sum(e_mdist.*permute(gt_em_meet_nm, [3 4 1 2]),[1 2]),ats,tpts); %From no manager to manager
        
        %Firm (a,q) Poached from firm with no manager
        Epgt_n_n=reshape(sum(e_ndist.*permute(gt_en_meet_nm, [3 4 1 2]), [1 2]),ats,tpts); %From no manager to no manager
        
        %Firm (a,q) Poached from firm with team
        Epgt_n_t=reshape(sum(e_tdist.*permute(gt_tf_meet_nm, [3 4 5 1 2]),[1 2 3]),ats,tpts); %From no manager to team
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Expected Gains of being poached for firm with team
        %Firm (a,z,q) Poached from empty firm
        Epgt_t_e=sum(e_edist'.*(gt_ef_meet_t),3); %From team to empty firm
        
        %Firm (a,z,q) Poached from firm with manager
        Epgt_t_m=reshape(sum(e_mdist.*permute(gt_em_meet_t, [4 5 1 2 3]),[1 2]),ats,tpts, tpts); %From team to manager   
        
        %Firm (a,z,q) Poached from firm with no manager
        Epgt_t_n=reshape(sum(e_ndist.*permute(gt_en_meet_t, [4 5 1 2 3]), [1 2]),ats,tpts, tpts); %From team to no manager
        
        %Firm (a,z,q) Poached from firm with team
        Epgt_t_t=reshape(sum(e_tdist.*permute(gt_tf_meet_t, [4 5 6 1 2 3]),[1 2 3]),ats,tpts, tpts); %From team to team
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % VF after the S&M
        
        % Empty firm
        Vetl=Veh+ (lamu*bpf/n)*Egt_e_u + (lam*bpf/n)*Egt_e_m + (lam*bpf/n)*Egt_e_n + (lam*bpf/n)*Egt_e_t;
        
        % Firm with manager
        Vmtl=Vmh+ del*(repmat(Veh',1, tpts)+repmat(Uh,ats,1)- Vmh) + (lamu*bpf/n)*Egt_m_u + (lam*bpf/n)*Egt_m_m + (lam*bpf/n)*Egt_m_n + (lam*bpf/n)*Egt_m_t...
            + (lam*bpw/n)*(Epgt_m_e + Epgt_m_m + Epgt_m_n + Epgt_m_t);
        
        % Firm with no manager
        Vntl=Vnh+ del*(repmat(Veh',1, tpts)+repmat(Uh,ats,1)- Vnh) + (lamu*bpf/n)*Egt_n_u + (lam*bpf/n)*Egt_n_m + (lam*bpf/n)*Egt_n_n + (lam*bpf/n)*Egt_n_t...
            + (lam*bpw/n)*(Epgt_n_e + Epgt_n_m + Epgt_n_n + Epgt_n_t);
        
        %Firm with team
        Vttl=Vth+ del*(repmat(Vmh,1,1,tpts)+repmat(Uh,ats,1,tpts)-Vth) + del*(repmat(Vnh,1,1,tpts)+repmat(Uh,ats,1,tpts)-Vth) + (lamu*bpf/n)*Egt_t_u + (lam*bpf/n)*Egt_t_m + (lam*bpf/n)*Egt_t_n + (lam*bpf/n)*Egt_t_t...
            + (lam*bpw/n)*(Epgt_t_e + Epgt_t_m + Epgt_t_n + Epgt_t_t);
        
        %Unemployed
        Utl= (lamu*bpw/n)*(Epgt_u_e + Epgt_u_m + Epgt_u_n + Epgt_u_t);  %Not dead sure this is right
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Vup - Productions and shocks VFs
        
        %Empty firm
        Veup= fe + bt*(a_trans*Vetl')';
        
        %Firm with manager
        Vmup=fman + bt*death*repmat(a_trans*Vetl',1,tpts)+ bt*(1-death)*(a_trans*Vmtl);
        
        % %Firm with no manager (This is q-trans do not depending on a. If it depends will have to include pages)
        % Vnup=fnman+bt*death*repmat(a_trans*Vetl',1,tpts)+ bt*(1-death)*(q_trans*(a_trans*Vntl)')';
        for a=1:ats
            for q=1:tpts
                expectdeath=zeros(1,ats);
                expectlive=zeros(ats,tpts);
                for apirme=1:ats
                    expectdeath(apirme)=a_trans(a,apirme)*Vetl(apirme);
                    for qprime=1:tpts
                        expectlive(apirme,qprime)=a_trans(a,apirme)*q_trans(q,qprime,a)*Vntl(apirme,qprime);
                    end
                end
                Vnup(a,q)=fnman(a,q)+bt*death*sum(expectdeath)+bt*(1-death)*sum(sum(expectlive));
            end
        end

        
        % %Firm with team (old version, q trans is not dependent on a)
        % %Transitions for full firm
        % %Conditional on q, matrix mult to get the expected value for (a,z). Do that for every q 
        % Va_prime=zeros(ats,tpts,tpts);
        % for i=1:tpts 
        %     Va_prime(:,:,i)=a_trans*Vttl(:,:,i);
        % end
        % %Permute to multiply it by q_trans
        % Va_prime=permute(Va_prime, [3 2 1]);
        
        % %Conditional on a, matrix mult to get the expected value for (z,q). Do that for every a
        % Vq_prime=zeros(tpts,tpts,ats);
        % for i=1:ats 
        %     Vq_prime(:,:,i)=q_trans*Va_prime(:,:,i);
        % end
        % Vqa_prime=permute(Vq_prime, [3 2 1]);
        
        % Vtup=fteam+ bt*(death.^2)*repmat(a_trans*Vetl',1,tpts,tpts)+ bt*death*(1-death)*repmat((a_trans*Vmtl),1,1,tpts)...
        %     + bt*death*(1-death)*permute(repmat((q_trans*(a_trans*Vntl)')',1,1,tpts),[1 3 2])...
        %     + bt*((1-death).^2)*Vqa_prime;

        %Firm with team (new version)
        for a=1:ats
            for z=1:tpts
                for q=1:tpts
                    expectdeath_b=zeros(1,ats);
                    expectdeath_n=zeros(1,ats);
                    expectdeath_m=zeros(ats,tpts);
                    expectlive_b=zeros(ats,tpts);
                    for apirme=1:ats
                        expectdeath_b(apirme)=a_trans(a,apirme)*Vetl(apirme);
                        expectdeath_n(apirme)=a_trans(a,apirme)*Vmtl(apirme,z);
                        for qprime=1:tpts
                            expectdeath_m(apirme,qprime)=q_trans(q,qprime,a)*a_trans(a,apirme)*Vntl(apirme,qprime);
                            expectlive_b(apirme,qprime)=q_trans(q,qprime,a)*a_trans(a,apirme)*Vttl(apirme,z,qprime);
                        end
                    end
                    Vtup(a,z,q)=fteam(a,z,q)+bt*(death.^2)*sum(expectdeath_b)+bt*death*(1-death)*sum(expectdeath_n)...
                        + bt*death*(1-death)*sum(sum(expectdeath_m)) + bt*((1-death).^2)*sum(sum(expectlive_b));
                end
            end
        end

    
        %Unemployed worker
        Uplus= Uh + death*(-Uh) + (1-death)*Utl;
        Uup=  b + bt*(u_trans*Uplus')';
    
    
        %diff=max([max(abs(V0 - V0up)),max(max(abs(V1 - V1up))),max(max(abs(V2 - V2up))),max(abs(U-Uup))]);
        diff=max([max(abs(Ve-Veup)), max(abs(Vm-Vmup),[],[1 2]), max(abs(Vn-Vnup),[],[1 2]), max(abs(Vt-Vtup),[],[1 2 3]), max(abs(U-Uup))]);  
    
        %Print every 50 iterations the difference and some text
        % if mod(it,100)==0
        %     fprintf('Iteration %d, error %f \n', it, diff)
        % end
        %In red failed to converge
        if it==itmax
            fprintf(2,'VF Failed to converge\n')
        end
        % if diff<diffmax
        %     cprintf('green','Converged in %d iterations\n',it)
        % end
    end
    end 
    
    Ve=Veup;
    Vm=Vmup;
    Vn=Vnup;
    Vt=Vtup;
    U=Uup;
    
    %Export the hats as well    %Fire/Promotion VFs
    % %Empty firm
    Veh=Ve;    
    %Firm with manager (a,z)
    Vm_alloc(:,:,1)=Vm; %Keep as it is
    Vm_alloc(:,:,2)=Vn-cost_d; %Demote z
    Vm_alloc(:,:,3)=repmat(Ve',1, tpts)+repmat(U,ats,1);%Fire z
    Vmh=max(Vm_alloc,[],3);
    
    
    %Firm with non-manager (a,q) 
    Vn_alloc(:,:,1)=Vn; %Keep as it is
    Vn_alloc(:,:,2)=Vm-cost_p; %Promote q
    Vn_alloc(:,:,3)=repmat(Ve',1, tpts)+repmat(U,ats,1); %Fire q
    Vnh=max(Vn_alloc,[],3);
    
    %Firm with team (a,z,q)
    Vt_alloc(:,:,:,1)=Vt; %Keep as it is
    Vt_alloc(:,:,:,2)=permute(repmat(Vm,1,1,tpts),[1 3 2])+repmat(U,ats,1,tpts)-cost_p; %Fire z and promote q
    Vt_alloc(:,:,:,3)=repmat(Vn,1,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]) - cost_d; %Fire q and demote z
    Vt_alloc(:,:,:,4)=repmat(Vm,1,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]); %Fire q 
    Vt_alloc(:,:,:,5)=permute(repmat(Vn,1,1,tpts),[1 3 2])+repmat(U,ats,1,tpts); %Fire z
    Vt_alloc(:,:,:,6)=permute(Vt, [1 3 2])-cost_d-cost_p; %Swap z and q with promotion and demotion costs
    Vt_alloc(:,:,:,7)=permute(repmat(Ve',1,tpts,tpts),[1 3 2])+repmat(U,ats,1,tpts)+permute(repmat(U,ats,1,tpts), [1 3 2]); %Fire z and q
    Vth=max(Vt_alloc,[],4);

    % And tildas    
    % VF after the S&M
        
    % Empty firm
    Vetl=Veh+ (lamu*bpf/n)*Egt_e_u + (lam*bpf/n)*Egt_e_m + (lam*bpf/n)*Egt_e_n + (lam*bpf/n)*Egt_e_t;
        
     % Firm with manager
    Vmtl=Vmh+ del*(repmat(Veh',1, tpts)+repmat(Uh,ats,1)- Vmh) + (lamu*bpf/n)*Egt_m_u + (lam*bpf/n)*Egt_m_m + (lam*bpf/n)*Egt_m_n + (lam*bpf/n)*Egt_m_t...
        + (lam*bpw/n)*(Epgt_m_e + Epgt_m_m + Epgt_m_n + Epgt_m_t);
        
    % Firm with no manager
    Vntl=Vnh+ del*(repmat(Veh',1, tpts)+repmat(Uh,ats,1)- Vnh) + (lamu*bpf/n)*Egt_n_u + (lam*bpf/n)*Egt_n_m + (lam*bpf/n)*Egt_n_n + (lam*bpf/n)*Egt_n_t...
        + (lam*bpw/n)*(Epgt_n_e + Epgt_n_m + Epgt_n_n + Epgt_n_t);
        
    %Firm with team
    Vttl=Vth+ del*(repmat(Vmh,1,1,tpts)+repmat(Uh,ats,1,tpts)-Vth) + del*(repmat(Vnh,1,1,tpts)+repmat(Uh,ats,1,tpts)-Vth) + (lamu*bpf/n)*Egt_t_u + (lam*bpf/n)*Egt_t_m + (lam*bpf/n)*Egt_t_n + (lam*bpf/n)*Egt_t_t...
        + (lam*bpw/n)*(Epgt_t_e + Epgt_t_m + Epgt_t_n + Epgt_t_t);
        
    %Unemployed
    Utl= (lamu*bpw/n)*(Epgt_u_e + Epgt_u_m + Epgt_u_n + Epgt_u_t);  %Not dead sure this is right

    

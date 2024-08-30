%Function for the VF interation of wages
function [Wm,Wn,Wtm,Wtn,Wmh,Wnh,Wtnh,Wtmh] =wf_iteration(wpts,ats,tpts,Ve,Vm,Vn,Vt,U,Vmh,Vnh,Vth,Veh,Wmini,Wnini,Wtmini,Wtnini,...
    speed,cost_d,cost_p,wgrid,bt,death,del,lam,lamu,bpw,n,a_trans,q_trans,e_udist,e_edist,e_mdist,e_ndist,e_tdist);
    Uh=U;
    %% Policy functions    
    % %Hiring policies, looking at origins not allocations yet 
    [h_e_u, h_e_m, h_e_nm, h_e_t_m, h_e_t_nm, h_m_u, h_m_m, h_m_nm, h_m_t_m, h_m_t_nm, h_nm_u, h_nm_m, h_nm_nm, h_nm_t_m, h_nm_t_nm, h_t_u, h_t_m, h_t_nm, h_t_t_m, h_t_t_nm]...
            =hire_policies(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,true,cost_d,cost_p);
            
    % Allocation Policies after hiring at SM
    [p_e_m, p_e_n, p_m_m_u, p_m_m_d, p_m_n, p_n_n_u, p_n_n_p, p_n_m, p_t_m_u, p_t_m_d, p_t_n_u, p_t_n_p]...
            =alloc_policies(Vmh,Vnh,Uh,Vth,ats,tpts,cost_d,cost_p);
    
    [d_m, r_m, d_n, r_n, d_t_m, d_t_n, r_t_m, r_t_n, d_t_b]...
        =reallocation_policies(Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p);
    %Inital guesses        
    % Wmup=zeros(wpts,ats,tpts); %Value of manaer in a manager only firm
    % Wnup=zeros(wpts,ats,tpts); % Value of non-manager in a non-manager only firm
    % Wtmup=zeros(wpts,ats,tpts,tpts); %Value of manager in a team firm
    % Wtnup=zeros(wpts,ats,tpts,tpts); %Value of non-manager in a team firm 
    Wmup=Wmini;
    Wnup=Wnini;
    Wtmup=Wtmini;
    Wtnup=Wtnini;


    %Iteration parameters
    diff=100;
    diffmax=1e-8;
         
    itw=0; 
    itmax=100000;

    %Loop until convergence
    while diff>diffmax & itw<itmax
        itw=itw+1;
        if itw==1
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %Update stage
           Wm=Wmup;
           Wn=Wnup;
           Wtm=Wtmup;
           Wtn=Wtnup;
        end
       
        if itw>1
             Wm=(1-speed)*Wm + speed*Wmup;
             Wn=(1-speed)*Wn + speed*Wnup;         
             Wtm=(1-speed)*Wtm + speed*Wtmup;
             Wtn=(1-speed)*Wtn + speed*Wtnup;
        end
        %% Reallocation stage
        % %Guesses for W (production values)
        Wmh=zeros(wpts,ats,tpts); %Value of manaer in a manager only firm
        Wnh=zeros(wpts,ats,tpts); % Value of non-manager in a non-manager only firm
        Wtmh=zeros(wpts,ats,tpts,tpts); %Value of manager in a team firm
        Wtnh=zeros(wpts,ats,tpts,tpts); %Value of non-manager in a team firm 
        %value for manager in a manager only firm
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    inner=[U(z), r_m(a,z)*Wn(w,a,z) +(1-r_m(a,z))*min(d_m(a,z)*U(z)+(1-d_m(a,z))*Wm(w,a,z),Vm(a,z)-Ve(a))];
                    Wmh(w,a,z)= max(inner);
                end
            end
        end
        
        %value for non-manager in a non-manager only firm
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    inner=[U(z), r_n(a,z)*Wm(w,a,z) +(1-r_n(a,z))*min(d_n(a,z)*U(z)+(1-d_n(a,z))*Wn(w,a,z),Vn(a,z)-Ve(a))];
                    Wnh(w,a,z)= max(inner);
                end
            end
        end
                    
        %value for manager in a team firm
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    for q=1:tpts
                        inner=[U(z), d_t_n(a,z,q)*(r_t_m(a,z,q)*Wn(w,a,z) +(1-r_t_m(a,z,q))*Wm(w,a,z))+r_t_n(a,z,q)*r_t_m(a,z,q)*Wtn(w,a,q,z) +...
                         (1-d_t_n(a,z,q)-r_t_m(a,z,q)*r_t_n(a,z,q))*min((d_t_m(a,z,q)+d_t_b(a,z,q))*U(z)+(1-d_t_m(a,z,q)-d_t_b(a,z,q))*Wtm(w,a,z,q),Vt(a,z,q)-Vn(a,q))];
                         Wtmh(w,a,z,q)= max(inner);
                    end
                end
            end
        end
        
        
        %value for non-manager in a team firm
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    for q=1:tpts
                        inner=[U(z), d_t_m(a,z,q)*(r_t_n(a,z,q)*Wm(w,a,q) +(1-r_t_n(a,z,q))*Wn(w,a,q))+r_t_n(a,z,q)*r_t_m(a,z,q)*Wtm(w,a,q,z) +...
                         (1-d_t_m(a,z,q)-r_t_n(a,z,q)*r_t_m(a,z,q))*min((d_t_n(a,z,q)+d_t_b(a,z,q))*U(q)+(1-d_t_n(a,z,q)-d_t_b(a,z,q))*Wtn(w,a,z,q),Vt(a,z,q)-Vm(a,z))];
                            Wtnh(w,a,z,q)= max(inner);
                    end
                end
            end
        end
        
        
        
        %% Search and Match Value
        Wmtl=zeros(wpts,ats,tpts); %Value of manaer in a manager only firm
        Wntl=zeros(wpts,ats,tpts); % Value of non-manager in a non-manager only firm
        Wtmtl=zeros(wpts,ats,tpts,tpts); %Value of manager in a team firm
        Wtntl=zeros(wpts,ats,tpts,tpts); %Value of non-manager in a team firm 
        
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    %Value for manager in a manager only firm (a,z)
                    % If meets from unemployed
                    store_u=zeros(1,tpts);
                    for ztilda=1:tpts
                        store_u(ztilda)=e_udist(ztilda)*h_m_u(ztilda,a,z)*(p_m_m_u(a,z,ztilda)*U(z)+p_m_m_d(a,z,ztilda)*Wtnh(w,a,ztilda,z)+p_m_n(a,z,ztilda)*Wtmh(w,a,z,ztilda)-Wmh(w,a,z));
                    end    
                    
                    %If meets from manager
                    store_m=zeros(ats,tpts);
                    for ztilda=1:tpts
                        for atilda=1:ats
                            store_m(atilda,ztilda)=e_mdist(atilda,ztilda)*h_m_m(atilda,ztilda,a,z)*(p_m_m_u(a,z,ztilda)*U(z)+p_m_m_d(a,z,ztilda)*Wtnh(w,a,ztilda,z)+p_m_n(a,z,ztilda)*Wtmh(w,a,z,ztilda)-Wmh(w,a,z));
                        end
                    end
                    
                    %If meets from non-manager
                    store_n=zeros(ats,tpts);
                    for ztilda=1:tpts
                        for atilda=1:ats
                            store_n(atilda,ztilda)=e_ndist(atilda,ztilda)*h_m_nm(atilda,ztilda,a,z)*(p_m_m_u(a,z,ztilda)*U(z)+p_m_m_d(a,z,ztilda)*Wtnh(w,a,ztilda,z)+p_m_n(a,z,ztilda)*Wtmh(w,a,z,ztilda)-Wmh(w,a,z));
                        end
                    end 
                    
                    %If meets from team
                    store_t=zeros(ats,tpts,tpts);
                    for ztilda=1:tpts
                        for atilda=1:ats
                            for qtilda=1:tpts
                                store_t(atilda,ztilda,qtilda)=e_tdist(atilda,ztilda,qtilda)*(h_m_t_m(atilda,ztilda,qtilda,a,z)*(p_m_m_u(a,z,ztilda)*U(z)+p_m_m_d(a,z,ztilda)*Wtnh(w,a,ztilda,z)+p_m_n(a,z,ztilda)*Wtmh(w,a,z,ztilda)-Wmh(w,a,z))...
                                    +h_m_t_nm(atilda,ztilda,qtilda,a,z)*(p_m_m_u(a,z,qtilda)*U(z)+p_m_m_d(a,z,qtilda)*Wtnh(w,a,qtilda,z)+p_m_n(a,z,qtilda)*Wtmh(w,a,z,qtilda)-Wmh(w,a,z)));
                            end 
                        end
                    end
                    
                    %Probablity of z being found and either poached or not (this part is newer)
                    
                    %Poached by empty firm
                    u=Vmh(a,z)-Veh(a); %Outside option for current firm 
                    store_pe=zeros(1,ats);
                    for atilda=1:ats
                        nu=max(Vmh(atilda,z),Vnh(atilda,z))-Veh(atilda); %marginal value for empty firm
                        store_pe(atilda)=e_edist(atilda)*h_e_m(a,z,atilda)*(max(Wmh(w,a,z),min(nu,bpw*nu+(1-bpw)*u))-Wmh(w,a,z));
                    end
                    
                    %Poached by manager firm
                    store_pm=zeros(ats,tpts);
                    for ztilda=1:tpts
                        for atilda=1:ats
                            nu=[Vmh(atilda,z)+Uh(ztilda),Vth(atilda,z,ztilda)-cost_d,Vth(atilda,ztilda,z)]-Vmh(atilda,ztilda);
                            store_pm(atilda,ztilda)=e_mdist(atilda,ztilda)*h_m_m(a,z,atilda,ztilda)*(max(Wmh(w,a,z),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wmh(w,a,z));
                        end
                    end
                    
                    %Poached by non-manager firm
                    store_pn=zeros(ats,tpts);
                    for ztilda=1:tpts
                        for atilda=1:ats
                            nu=[Vnh(atilda,z)+Uh(ztilda),Vth(atilda,z,ztilda),Vth(atilda,ztilda,z)-cost_p]-Vnh(atilda,ztilda);
                            store_pn(atilda,ztilda)=e_ndist(atilda,ztilda)*h_nm_m(a,z,atilda,ztilda)*(max(Wmh(w,a,z),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wmh(w,a,z));
                        end
                    end
                    
                    %Poached by team firm
                    store_pt=zeros(ats,tpts,tpts);
                    for ztilda=1:tpts
                        for atilda=1:ats
                            for qtilda=1:tpts
                                nu=[Vth(atilda,z,ztilda)+Uh(qtilda)-cost_d,Vth(atilda,ztilda,z)+Uh(qtilda),Vth(atilda,qtilda,z)-cost_p+Uh(ztilda), Vth(atilda,z,qtilda)+Uh(ztilda)]-Vth(atilda,ztilda,qtilda);
                                store_pt(atilda,ztilda,qtilda)=e_tdist(atilda,ztilda,qtilda)*h_t_m(a,z,atilda,ztilda,qtilda)*(max(Wmh(w,a,z),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wmh(w,a,z));
                            end
                        end
                    end
                    
                    %Value function
                    Wmtl(w,a,z)= Wmh(w,a,z)+del*(U(z)-Wmh(w,a,z))+...
                        (lamu/n)*sum(store_u)+...
                        (lam/n)*sum(store_m,"all")+...
                        (lam/n)*sum(store_n,"all")+...
                        (lam/n)*sum(store_t,"all")+...
                        (lam/n)*sum(store_pe)+...
                        (lam/n)*sum(store_pm,"all")+...
                        (lam/n)*sum(store_pn,"all")+...
                        (lam/n)*sum(store_pt,"all");
                    
                    %Value for non-manager in a non-manager only firm (a,z)
        
                    % If meets from unemployed
                    store_u=zeros(1,tpts);
                    for ztilda=1:tpts
                        store_u(ztilda)=e_udist(ztilda)*h_nm_u(ztilda,a,z)*(p_n_n_u(a,z,ztilda)*U(z)+p_n_n_p(a,z,ztilda)*Wtmh(w,a,z,ztilda)+p_n_m(a,z,ztilda)*Wtnh(w,a,z,ztilda)-Wnh(w,a,z));
                    end
                    
                    %If meets from manager
                    store_m=zeros(ats,tpts);
                    for ztilda=1:tpts
                        for atilda=1:ats
                            store_m(atilda,ztilda)=e_mdist(atilda,ztilda)*h_nm_m(atilda,ztilda,a,z)*(p_n_n_u(a,z,ztilda)*U(z)+p_n_n_p(a,z,ztilda)*Wtmh(w,a,z,ztilda)+p_n_m(a,z,ztilda)*Wtnh(w,a,z,ztilda)-Wnh(w,a,z));
                        end
                    end
                    
                    %If meets from non-manager
                    store_n=zeros(ats,tpts);
                    for ztilda=1:tpts
                        for atilda=1:ats
                            store_n(atilda,ztilda)=e_ndist(atilda,ztilda)*h_nm_nm(atilda,ztilda,a,z)*(p_n_n_u(a,z,ztilda)*U(z)+p_n_n_p(a,z,ztilda)*Wtmh(w,a,z,ztilda)+p_n_m(a,z,ztilda)*Wtnh(w,a,z,ztilda)-Wnh(w,a,z));
                        end
                    end
                    
                    %If meets from team
                    store_t=zeros(ats,tpts,tpts);
                    for ztilda=1:tpts
                        for atilda=1:ats
                            for qtilda=1:tpts
                                store_t(atilda,ztilda,qtilda)=e_tdist(atilda,ztilda,qtilda)*(h_nm_t_m(atilda,ztilda,qtilda,a,z)*(p_n_n_u(a,z,ztilda)*U(z)+p_n_n_p(a,z,ztilda)*Wtmh(w,a,z,ztilda)+p_n_m(a,z,ztilda)*Wtnh(w,a,z,ztilda)-Wnh(w,a,z))...
                                    +h_nm_t_nm(atilda,z,qtilda,a,z)*(p_n_n_u(a,z,qtilda)*U(z)+p_n_n_p(a,z,qtilda)*Wtmh(w,a,z,qtilda)+p_n_m(a,z,qtilda)*Wtnh(w,a,z,qtilda)-Wnh(w,a,z)));
                            end
                        end
                    end
                    
                    %Probablity of z being found and either poached or not (this part is newer)
                    %Poached by empty firm
                    u=Vnh(a,z)-Veh(a); %Outside option for current firm
                    store_pe=zeros(1,ats);
                    for atilda=1:ats
                        nu=max(Vmh(atilda,z),Vnh(atilda,z))-Veh(atilda); %marginal value for empty firm
                        store_pe(atilda)=e_edist(atilda)*h_e_nm(a,z,atilda)*(max(Wnh(w,a,z),min(nu,bpw*nu+(1-bpw)*u))-Wnh(w,a,z));
                    end
                    
                    %Poached by manager firm
                    store_pm=zeros(ats,tpts);
                    for ztilda=1:tpts
                        for atilda=1:ats
                            nu=[Vmh(atilda,z)+Uh(ztilda),Vth(atilda,z,ztilda)-cost_d,Vth(atilda,ztilda,z)]-Vmh(atilda,ztilda);
                            store_pm(atilda,ztilda)=e_mdist(atilda,ztilda)*h_m_nm(a,z,atilda,ztilda)*(max(Wnh(w,a,z),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wnh(w,a,z));
                        end
                    end
                    
                    %Poached by non-manager firm
                    store_pn=zeros(ats,tpts);   
                    for ztilda=1:tpts
                        for atilda=1:ats
                            nu=[Vnh(atilda,z)+Uh(ztilda),Vth(atilda,z,ztilda),Vth(atilda,ztilda,z)-cost_p]-Vnh(atilda,ztilda);
                            store_pn(atilda,ztilda)=e_ndist(atilda,ztilda)*h_nm_nm(a,z,atilda,ztilda)*(max(Wnh(w,a,z),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wnh(w,a,z));
                        end
                    end
                    
                    %Poached by team firm
                    store_pt=zeros(ats,tpts,tpts);
                    for ztilda=1:tpts
                        for atilda=1:ats
                            for qtilda=1:tpts
                                nu=[Vth(atilda,z,ztilda)+Uh(qtilda)-cost_d,Vth(atilda,ztilda,z)+Uh(qtilda),Vth(atilda,qtilda,z)-cost_p+Uh(ztilda), Vth(atilda,z,qtilda)+Uh(ztilda)]-Vth(atilda,ztilda,qtilda);
                                store_pt(atilda,ztilda,qtilda)=e_tdist(atilda,ztilda,qtilda)*h_t_nm(a,z,atilda,ztilda,qtilda)*(max(Wnh(w,a,z),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wnh(w,a,z));
                            end
                        end
                    end
                    
                    %Value function
                    Wntl(w,a,z)= Wnh(w,a,z)+del*(U(z)-Wnh(w,a,z))+...
                        (lamu/n)*sum(store_u)+...
                        (lam/n)*sum(store_m,"all")+...
                        (lam/n)*sum(store_n,"all")+...
                        (lam/n)*sum(store_t,"all")+...
                        (lam/n)*sum(store_pe)+...
                        (lam/n)*sum(store_pm,"all")+...
                        (lam/n)*sum(store_pn,"all")+...
                        (lam/n)*sum(store_pt,"all");
                    
                    
                    for q=1:tpts
                        % Value for manager in a team firm (a,z,q)
                        % If meets from unemployed
                        store_u=zeros(1,tpts);
                        for ztilda=1:tpts
                            store_u(ztilda)=e_udist(ztilda)*h_t_u(ztilda,a,z,q)*((p_t_m_u(a,z,q,ztilda)+p_t_n_p(a,z,q,ztilda))*U(z)+p_t_m_d(a,z,q,ztilda)*Wtnh(w,a,ztilda,z)+p_t_n_u(a,z,q,ztilda)*Wtmh(w,a,z,ztilda)-Wtmh(w,a,z,q));   
                        end
                        
                        %If meets from manager
                        store_m=zeros(ats,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                store_m(atilda,ztilda)=e_mdist(atilda,ztilda)*h_t_m(atilda,ztilda,a,z,q)*((p_t_m_u(a,z,q,ztilda)+p_t_n_p(a,z,q,ztilda))*U(z)+p_t_m_d(a,z,q,ztilda)*Wtnh(w,a,ztilda,z)+p_t_n_u(a,z,q,ztilda)*Wtmh(w,a,z,ztilda)-Wtmh(w,a,z,q));
                            end
                        end
                        
                        %If meets from non-manager
                        store_n=zeros(ats,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                store_n(atilda,ztilda)=e_ndist(atilda,ztilda)*h_t_nm(atilda,ztilda,a,z,q)*((p_t_m_u(a,z,q,ztilda)+p_t_n_p(a,z,q,ztilda))*U(z)+p_t_m_d(a,z,q,ztilda)*Wtnh(w,a,ztilda,z)+p_t_n_u(a,z,q,ztilda)*Wtmh(w,a,z,ztilda)-Wtmh(w,a,z,q));
                            end
                        end 
                        
                        %If meets from team
                        store_t=zeros(ats,tpts,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                for qtilda=1:tpts
                                    store_t(atilda,ztilda,qtilda)=e_tdist(atilda,ztilda,qtilda)*(h_t_t_m(atilda,ztilda,qtilda,a,z,q)*((p_t_m_u(a,z,q,ztilda)+p_t_n_p(a,z,q,ztilda))*U(z)+p_t_m_d(a,z,q,ztilda)*Wtnh(w,a,ztilda,z)+p_t_n_u(a,z,q,ztilda)*Wtmh(w,a,z,ztilda)-Wtmh(w,a,z,q))...
                                        +h_t_t_nm(atilda,ztilda,qtilda,a,z,q)*((p_t_m_u(a,z,q,qtilda)+p_t_n_p(a,z,q,qtilda))*U(z)+p_t_m_d(a,z,q,qtilda)*Wtnh(w,a,qtilda,z)+p_t_n_u(a,z,q,qtilda)*Wtmh(w,a,z,qtilda)-Wtmh(w,a,z,q)));
                                end 
                            end
                        end
                        
                        %Probablity of z being found and either poached or not (this part is newer)
                        %Poached by empty firm
                        u=Vth(a,z,q)-Vnh(a,q); %Outside option for current firm of losing a manager
                        store_pe=zeros(1,ats);
                        for atilda=1:ats
                            nu=max(Vmh(atilda,z),Vnh(atilda,z))-Veh(atilda); %marginal value for empty firm
                            store_pe(atilda)=e_edist(atilda)*h_e_t_m(a,z,q,atilda)*(max(Wtmh(w,a,z,q),min(nu,bpw*nu+(1-bpw)*u))-Wtmh(w,a,z,q));
                        end
                        
                        %Poached by manager firm
                        store_pm=zeros(ats,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                nu=[Vmh(atilda,z)+Uh(ztilda),Vth(atilda,z,ztilda)-cost_d,Vth(atilda,ztilda,z)]-Vmh(atilda,ztilda);
                                store_pm(atilda,ztilda)=e_mdist(atilda,ztilda)*h_m_t_m(a,z,q,atilda,ztilda)*(max(Wtmh(w,a,z,q),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wtmh(w,a,z,q));
                            end
                        end
                        
                        %Poached by non-manager firm
                        store_pn=zeros(ats,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                nu=[Vnh(atilda,z)+Uh(ztilda),Vth(atilda,z,ztilda),Vth(atilda,ztilda,z)-cost_p]-Vnh(atilda,ztilda);
                                store_pn(atilda,ztilda)=e_ndist(atilda,ztilda)*h_nm_t_m(a,z,q,atilda,ztilda)*(max(Wtmh(w,a,z,q),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wtmh(w,a,z,q));
                            end
                        end
                        
                        %Poached by team firm
                        store_pt=zeros(ats,tpts,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                for qtilda=1:tpts
                                    nu=[Vth(atilda,z,ztilda)+Uh(qtilda)-cost_d,Vth(atilda,ztilda,z)+Uh(qtilda),Vth(atilda,qtilda,z)-cost_p+Uh(ztilda), Vth(atilda,z,qtilda)+Uh(ztilda)]-Vth(atilda,ztilda,qtilda);
                                    store_pt(atilda,ztilda,qtilda)=e_tdist(atilda,ztilda,qtilda)*h_t_t_m(a,z,q,atilda,ztilda,qtilda)*(max(Wtmh(w,a,z,q),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wtmh(w,a,z,q));
                                end
                            end
                        end
                        
                        %Value function
                        Wtmtl(w,a,z,q)= Wtmh(w,a,z,q)+del*(U(z)-Wtmh(w,a,z,q))+...
                            (lamu/n)*sum(store_u)+...
                            (lam/n)*sum(store_m,"all")+...
                            (lam/n)*sum(store_n,"all")+...
                            (lam/n)*sum(store_t,"all")+...
                            (lam/n)*sum(store_pe)+...
                            (lam/n)*sum(store_pm,"all")+...
                            (lam/n)*sum(store_pn,"all")+...
                            (lam/n)*sum(store_pt,"all");
                        
                        
                        
                        
                        %Value for non-manager in a team firm (a,z,q)
                        
                        % If meets from unemployed
                        store_u=zeros(1,tpts);
                        for ztilda=1:tpts
                            store_u(ztilda)=e_udist(ztilda)*h_t_u(ztilda,a,z,q)*((p_t_n_u(a,z,q,ztilda)+p_t_m_d(a,z,q,ztilda))*U(q)+p_t_n_p(a,z,q,ztilda)*Wtmh(w,a,q,ztilda)+p_t_m_u(a,z,q,ztilda)*Wtnh(w,a,ztilda,q)-Wtnh(w,a,z,q));
                        end
                        
                        %If meets from manager
                        store_m=zeros(ats,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                store_m(atilda,ztilda)=e_mdist(atilda,ztilda)*h_t_m(atilda,ztilda,a,z,q)*((p_t_n_u(a,z,q,ztilda)+p_t_m_d(a,z,q,ztilda))*U(q)+p_t_n_p(a,z,q,ztilda)*Wtmh(w,a,q,ztilda)+p_t_m_u(a,z,q,ztilda)*Wtnh(w,a,ztilda,q)-Wtnh(w,a,z,q));
                            end
                        end
                        
                        %If meets from non-manager
                        store_n=zeros(ats,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                store_n(atilda,ztilda)=e_ndist(atilda,ztilda)*h_t_nm(atilda,ztilda,a,z,q)*((p_t_n_u(a,z,q,ztilda)+p_t_m_d(a,z,q,ztilda))*U(q)+p_t_n_p(a,z,q,ztilda)*Wtmh(w,a,q,ztilda)+p_t_m_u(a,z,q,ztilda)*Wtnh(w,a,ztilda,q)-Wtnh(w,a,z,q));
                            end
                        end
                        
                        %If meets from team
                        store_t=zeros(ats,tpts,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                for qtilda=1:tpts
                                    store_t(atilda,ztilda,qtilda)=e_tdist(atilda,ztilda,qtilda)*(h_t_t_m(atilda,ztilda,qtilda,a,z,q)*((p_t_n_u(a,z,q,ztilda)+p_t_m_d(a,z,q,ztilda))*U(q)+p_t_n_p(a,z,q,ztilda)*Wtmh(w,a,q,ztilda)+p_t_m_u(a,z,q,ztilda)*Wtnh(w,a,ztilda,q)-Wtnh(w,a,z,q))...
                                        +h_t_t_nm(atilda,ztilda,qtilda,a,z,q)*((p_t_n_u(a,z,q,qtilda)+p_t_m_d(a,z,q,qtilda))*U(q)+p_t_n_p(a,z,q,qtilda)*Wtmh(w,a,q,qtilda)+p_t_m_u(a,z,q,qtilda)*Wtnh(w,a,qtilda,q)-Wtnh(w,a,z,q)));
                                end 
                            end
                        end
                        
                        %Probablity of q being found and either poached or not (this part is newer)
                        %Poached by empty firm
                        u=Vth(a,z,q)-Vmh(a,z); %Outside option for current firm
                        store_pe=zeros(1,ats);
                        for atilda=1:ats
                            nu=max(Vmh(atilda,q),Vnh(atilda,q))-Veh(atilda); %marginal value for empty firm
                            store_pe(atilda)=e_edist(atilda)*h_e_t_nm(a,z,q,atilda)*(max(Wtnh(w,a,z,q),min(nu,bpw*nu+(1-bpw)*u))-Wtnh(w,a,z,q));
                        end 
                        
                        %Poached by manager firm
                        store_pm=zeros(ats,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                nu=[Vmh(atilda,q)+Uh(ztilda),Vth(atilda,q,ztilda)-cost_d,Vth(atilda,ztilda,q)]-Vmh(atilda,ztilda);
                                store_pm(atilda,ztilda)=e_mdist(atilda,ztilda)*h_m_t_nm(a,z,q,atilda,ztilda)*(max(Wtnh(w,a,z,q),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wtnh(w,a,z,q));
                            end
                        end
                        
                        %Poached by non-manager firm
                        store_pn=zeros(ats,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                nu=[Vnh(atilda,q)+Uh(ztilda),Vth(atilda,q,ztilda),Vth(atilda,ztilda,q)-cost_p]-Vnh(atilda,ztilda);
                                store_pn(atilda,ztilda)=e_ndist(atilda,ztilda)*h_nm_t_nm(a,z,q,atilda,ztilda)*(max(Wtnh(w,a,z,q),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wtnh(w,a,z,q));
                            end
                        end
                        
                        %Poached by team firm
                        store_pt=zeros(ats,tpts,tpts);
                        for ztilda=1:tpts
                            for atilda=1:ats
                                for qtilda=1:tpts
                                    nu=[Vth(atilda,q,ztilda)+Uh(qtilda)-cost_d,Vth(atilda,ztilda,q)+Uh(qtilda),Vth(atilda,qtilda,q)-cost_p+Uh(ztilda), Vth(atilda,q,qtilda)+Uh(ztilda)]-Vth(atilda,ztilda,qtilda);
                                    store_pt(atilda,ztilda,qtilda)=e_tdist(atilda,ztilda,qtilda)*h_t_t_nm(a,z,q,atilda,ztilda,qtilda)*(max(Wtnh(w,a,z,q),min(max(nu),bpw*max(nu)+(1-bpw)*u))-Wtnh(w,a,z,q));
                                end
                            end
                        end
                        
                        %Value function
                        Wtntl(w,a,z,q)= Wtnh(w,a,z,q)+del*(U(q)-Wtnh(w,a,z,q))+...
                            (lamu/n)*sum(store_u)+...
                            (lam/n)*sum(store_m,"all")+...
                            (lam/n)*sum(store_n,"all")+...
                            (lam/n)*sum(store_t,"all")+...
                            (lam/n)*sum(store_pe)+...
                            (lam/n)*sum(store_pm,"all")+...
                            (lam/n)*sum(store_pn,"all")+...
                            (lam/n)*sum(store_pt,"all");
                    end
                end
            end
        end
        
        
        
        
           
        %% Shocks Stage
        %Up values
        Wmup=zeros(wpts,ats,tpts); %Value of manaer in a manager only firm
        Wnup=zeros(wpts,ats,tpts); % Value of non-manager in a non-manager only firm
        Wtmup=zeros(wpts,ats,tpts,tpts); %Value of manager in a team firm
        Wtnup=zeros(wpts,ats,tpts,tpts); %Value of non-manager in a team firm
        
        
        % For firm with manager only
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    expec_m=zeros(1,ats);
                    for aprime=1:ats
                        expec_m(aprime)=a_trans(a,aprime)*Wmtl(w,aprime,z);
                    end
                    Wmup(w,a,z)= wgrid(w)+ bt*(1-death)*sum(expec_m);
                end
            end
        end
        
        % For firm with non-manager only
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    expec_n=zeros(ats,tpts);
                    for aprime=1:ats
                        for qprime=1:tpts
                            expec_n(aprime,qprime)=a_trans(a,aprime)*q_trans(z,qprime,a)*Wntl(w,aprime,qprime);
                        end
                    end
                    Wnup(w,a,z)= wgrid(w)+ bt*(1-death)*sum(expec_n,"all");
                end
            end
        end
        
        
        % For firm with team
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    for q=1:tpts
                        expec_tm=zeros(ats,tpts);
                        expec_tn=zeros(ats,tpts);
                        for aprime=1:ats
                            for qprime=1:tpts
                                expec_tm(aprime,qprime)=a_trans(a,aprime)*q_trans(q,qprime,a)*(death*(1-death)*Wmtl(w,aprime,z)+(1-death)^2*Wtmtl(w,aprime,z,qprime));
                                expec_tn(aprime,qprime)=a_trans(a,aprime)*q_trans(q,qprime,a)*(death*(1-death)*Wntl(w,aprime,qprime)+(1-death)^2*Wtntl(w,aprime,z,qprime));
                            end
                        end
                        Wtmup(w,a,z,q)= wgrid(w)+ bt*sum(expec_tm,"all");
                        Wtnup(w,a,z,q)= wgrid(w)+ bt*sum(expec_tn,"all");
                    end
                end
            end
        end
    
        diff=max([max(abs(Wmup-Wm),[],'all'),max(abs(Wnup-Wn),[],'all'),max(abs(Wtmup-Wtm),[],'all'),max(abs(Wtnup-Wtn),[],'all')]);
        %Print every 50 iterations the difference and some text
        % if mod(itw,100)==0
        %     fprintf('Wage Iteration %d, error %f \n', itw, diff)
        % end
        %In red failed to converge
        if itw==itmax
            fprintf(2,'WF Failed to converge\n')
        end
        if diff<diffmax
            fprintf('Wages Converged in %d iterations\n',itw)
        end
    end
    %Final values
    Wm=Wmup;
    Wn=Wnup;
    Wtm=Wtmup;
    Wtn=Wtnup;
    

    %Hat values
    %% Reallocation stage
        % %Guesses for W (production values)
        Wmh=zeros(wpts,ats,tpts); %Value of manaer in a manager only firm
        Wnh=zeros(wpts,ats,tpts); % Value of non-manager in a non-manager only firm
        Wtmh=zeros(wpts,ats,tpts,tpts); %Value of manager in a team firm
        Wtnh=zeros(wpts,ats,tpts,tpts); %Value of non-manager in a team firm 
        %value for manager in a manager only firm
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    inner=[U(z), r_m(a,z)*Wn(w,a,z) +(1-r_m(a,z))*min(d_m(a,z)*U(z)+(1-d_m(a,z))*Wm(w,a,z),Vm(a,z)-Ve(a))];
                    Wmh(w,a,z)= max(inner);
                end
            end
        end
        
        %value for non-manager in a non-manager only firm
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    inner=[U(z), r_n(a,z)*Wm(w,a,z) +(1-r_n(a,z))*min(d_n(a,z)*U(z)+(1-d_n(a,z))*Wn(w,a,z),Vn(a,z)-Ve(a))];
                    Wnh(w,a,z)= max(inner);
                end
            end
        end
                    
        %value for manager in a team firm
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    for q=1:tpts
                        inner=[U(z), d_t_n(a,z,q)*(r_t_m(a,z,q)*Wn(w,a,z) +(1-r_t_m(a,z,q))*Wm(w,a,z))+r_t_n(a,z,q)*r_t_m(a,z,q)*Wtn(w,a,q,z) +...
                         (1-d_t_n(a,z,q)-r_t_m(a,z,q)*r_t_n(a,z,q))*min((d_t_m(a,z,q)+d_t_b(a,z,q))*U(z)+(1-d_t_m(a,z,q)-d_t_b(a,z,q))*Wtm(w,a,z,q),Vt(a,z,q)-Vn(a,q))];
                         Wtmh(w,a,z,q)= max(inner);
                    end
                end
            end
        end
        
        
        %value for non-manager in a team firm
        for w=1:wpts
            for a=1:ats
                for z=1:tpts
                    for q=1:tpts
                        inner=[U(z), d_t_m(a,z,q)*(r_t_n(a,z,q)*Wm(w,a,q) +(1-r_t_n(a,z,q))*Wn(w,a,q))+r_t_n(a,z,q)*r_t_m(a,z,q)*Wtm(w,a,q,z) +...
                         (1-d_t_m(a,z,q)-r_t_n(a,z,q)*r_t_m(a,z,q))*min((d_t_n(a,z,q)+d_t_b(a,z,q))*U(q)+(1-d_t_n(a,z,q)-d_t_b(a,z,q))*Wtn(w,a,z,q),Vt(a,z,q)-Vm(a,z))];
                            Wtnh(w,a,z,q)= max(inner);
                    end
                end
            end
        end




    
        


        

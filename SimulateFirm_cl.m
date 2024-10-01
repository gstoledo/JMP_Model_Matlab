% %Simulation as a function
function fs=SimulateFirm_cl(p,ps,tg,sp,v,d,w)
    %Opening up the parameters
    fieldNames = fieldnames(p);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = p.' varName ';']);
    end

    fieldNames = fieldnames(ps);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = ps.' varName ';']);
    end


    fieldNames = fieldnames(tg);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = tg.' varName ';']);
    end

    fieldNames = fieldnames(sp);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = sp.' varName ';']);
    end

    fieldNames = fieldnames(v);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = v.' varName ';']);
    end

    fieldNames = fieldnames(d);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = d.' varName ';']);
    end

    fieldNames = fieldnames(w);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = w.' varName ';']);
    end

    %% Derivated from the main parameters
    [a_trans,q_trans,u_trans,fteam,fman,fnman,fe,b,mnew_high,typebirth,wmin,wmax,wgrid] = derivatated_p(p,ps,split_top_bot);

    %Policy functions
    %Get the policy functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Policy functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Hiring policies, looking at origins not allocations yet 
    [h_e_u, h_e_m, h_e_nm, h_e_t_m, h_e_t_nm, h_m_u, h_m_m, h_m_nm, h_m_t_m, h_m_t_nm, h_nm_u, h_nm_m, h_nm_nm, h_nm_t_m, h_nm_t_nm, h_t_u, h_t_m, h_t_nm, h_t_t_m, h_t_t_nm]...
        =hire_policies(Vmh,Vnh,Veh,U,Vth,ats,tpts,0,cost_d,cost_p);
        
    % Allocation Policies after hiring at SM
    [p_e_m, p_e_n, p_m_m_u, p_m_m_d, p_m_n, p_n_n_u, p_n_n_p, p_n_m, p_t_m_u, p_t_m_d, p_t_n_u, p_t_n_p]...
        =alloc_policies(Vmh,Vnh,U,Vth,ats,tpts,cost_d,cost_p);
    
    %Reallocation policies
    [d_m, r_m, d_n, r_n, d_t_m, d_t_n, r_t_m, r_t_n, d_t_b]...
        =reallocation_policies(Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p);
    
    
    
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
    % For a
    %Given a, vector of probabilities to go down, stay or up
    down_a=zeros(1,ats);
    stay_a=zeros(1,ats);
    up_a=zeros(1,ats);
    for a=1:ats
        if a==1
            down_a(a)=0;
        else
            down_a(a)=a_trans(a,a-1);
        end
        stay_a(a)=a_trans(a,a);
        if a==ats
            up_a(a)=0;
        else
            up_a(a)=a_trans(a,a+1);
        end
    end
    
    % For q
    down_q=zeros(ats,tpts);
    stay_q=zeros(ats,tpts);
    up_q=zeros(ats,tpts);
    for a=1:ats
        for q=1:tpts
            if q==1
                down_q(a,q)=0;
            else
                down_q(a,q)=q_trans(q,q-1,a);
            end
            stay_q(a,q)=q_trans(q,q,a);
            if q==tpts
                up_q(a,q)=0;
            else
                up_q(a,q)=q_trans(q,q+1,a);
            end
        end
    end
    
    % For u
    down_u=zeros(1,tpts);
    stay_u=zeros(1,tpts);
    up_u=zeros(1,tpts);
    for u=1:tpts
        if u==1
            down_u(u)=0;
        else
            down_u(u)=u_trans(u,u-1);
        end
        stay_u(u)=u_trans(u,u);
        if u==tpts
            up_u(u)=0;
        else
            up_u(u)=u_trans(u,u+1);
        end
    end
    %Random matrices
    rng(seed,"twister");
    
    r.state=rand(n_months,n_firms); %Inside the meetings, which type of meeting is happening
    r.event=rand(n_months,n_firms); % What kind of meeting is happening
    r.ashock=rand(n_months,n_firms); %Productivity shock
    r.qshock=rand(n_months,n_firms); %Learning shock
    r.type=rand(n_months,n_firms); %Type of match
    
    
    % Storage of firm status
    firm_status=zeros(n_months,n_firms); %(1=empty, 2= Firm with manager, 3= Firm with non manager; 4= Firm with team)
    m_ftp=zeros(n_months,n_firms); %Manager type
    n_ftp=zeros(n_months,n_firms); %Non manager type
    a_ftp=zeros(n_months,n_firms); %Productivity type
    m_fwage=zeros(n_months,n_firms); %Manager wage
    n_fwage=zeros(n_months,n_firms); %Non manager wage
    %Flagging the events
    hire_m=zeros(n_months,n_firms); %Hiring of manager ( 1=hire)
    hire_n=zeros(n_months,n_firms); %Hiring of non manager ( 1=hire)
    prom_sm=zeros(n_months,n_firms); %Promotion at SM (-1=demote, 0=stay, 1=promote)
    prom_realloc=zeros(n_months,n_firms); %Promotion at reallocation (-1=demote, 0=stay, 1=promote)
    fire_m=zeros(n_months,n_firms); %Firing of manager at reallocation ( 1=fire)
    fire_n=zeros(n_months,n_firms); %Firing of non manager at reallocation ( 1=fire)
    
    
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
            m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
            n_fwage(t,i)=0.0;
        elseif (firm_status(t,i)==3) %Firm with non manager
            [a_ftp(t,i),n_ftp(t,i)]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.event(t,i));
            m_ftp(t,i)=0;
            m_fwage(t,i)=0.0;
            target_wage=(1-bpw)*U(n_ftp(t,i))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
            n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
        else
            %Firm with team
            [a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.event(t,i));
            target_wage_m=(1-bpw)*U(m_ftp(t,i))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
            target_wage_n=(1-bpw)*U(n_ftp(t,i))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
            m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage_m,wgrid);
            n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage_n,wgrid);
        end
    end    
    
    %Other months
    for i=1:n_firms
        for t=2:n_months
            % fprintf('Firm %d, Month %d\n',i,t)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Shocks and SM
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % fprintf('SM for status 1')
            if (firm_status(t-1,i)==1) %Empty firm
                %Check productivity shock
                if (r.ashock(t,i)<down_a(a_ftp(t-1,i))) %Firm go down in a
                    a_ftp(t,i)=a_ftp(t-1,i)-1;
                elseif (r.ashock(t,i)< down_a(a_ftp(t-1,i))+stay_a(a_ftp(t-1,i)) && r.ashock(t,i)>=down_a(a_ftp(t-1,i))) %Firm stay in a
                    a_ftp(t,i)=a_ftp(t-1,i);
                else %Firm go up in a
                    a_ftp(t,i)=a_ftp(t-1,i)+1;
                end
               % Does it match someone?
               
               if (r.event(t,i)<prob_matching)
                %Match unemp
                if (r.state(t,i)<=cdf_emp_dist(1))
                    zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
                    % Hire unemp?
                    if (h_e_u(zu,a_ftp(t,i))>0)
                        %Where does it go?
                        if (p_e_m(a_ftp(t,i),zu)>0) %Put as manager
                            firm_status(t,i)=2;
                            hire_m(t,i)=1;
                            m_ftp(t,i)=zu;
                            n_ftp(t,i)=0;
                            target_wage=(1-bpw)*U(m_ftp(t,i))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
                            m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
                            n_fwage(t,i)=0.0;
                        else %Put as non manager
                            firm_status(t,i)=3;
                            hire_n(t,i)=1;
                            n_ftp(t,i)=zu;
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            target_wage=(1-bpw)*U(n_ftp(t,i))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
                            n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                        end
                    else %Firm does not hire
                        firm_status(t,i)=1;
                        m_ftp(t,i)=0;
                        n_ftp(t,i)=0;
                        m_fwage(t,i)=0.0;
                        n_fwage(t,i)=0.0;
                    end
                %Match employed
                elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<=cdf_emp_dist(2))
                    [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
                    % Hire employed manager?
                    if (h_e_m(am,zm,a_ftp(t,i))>0)
                        %Where does it go?
                        if (p_e_m(a_ftp(t,i),zm)>0) %Put as manager
                            firm_status(t,i)=2;
                            hire_m(t,i)=1;
                            m_ftp(t,i)=zm;
                            n_ftp(t,i)=0;
                            target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
                            m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
                            n_fwage(t,i)=0.0;
                        else %Put as non manager
                            firm_status(t,i)=3;
                            hire_n(t,i)=1;
                            n_ftp(t,i)=zm;
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
                            n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                        end
                    else %Firm does not hire
                        firm_status(t,i)=1;
                        m_ftp(t,i)=0;
                        n_ftp(t,i)=0;
                        m_fwage(t,i)=0.0;
                        n_fwage(t,i)=0.0;
                    end
                %Match employed non manager
                elseif (r.state(t,i)>cdf_emp_dist(2) && r.state(t,i)<=cdf_emp_dist(3))
                    [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
                    % Hire employed non manager?
                    if (h_e_nm(an,zn,a_ftp(t,i))>0)
                        %Where does it go?
                        if (p_e_m(a_ftp(i),zn)>0) %Put as manager
                            firm_status(t,i)=2;
                            hire_m(t,i)=1;
                            m_ftp(t,i)=zn;
                            n_ftp(t,i)=0;
                            n_fwage(t,i)=0.0;
                            target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
                            m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
                        else %Put as non manager
                            firm_status(t,i)=3;
                            hire_n(t,i)=1;
                            n_ftp(t,i)=zn;
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
                            n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                        end
                    else %Firm does not hire
                        firm_status(t,i)=1;
                        m_ftp(t,i)=0;
                        n_ftp(t,i)=0;
                        m_fwage(t,i)=0.0;
                        n_fwage(t,i)=0.0;
                    end
                %Match with team
                else 
                    [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
                    % Hire the manager?
                    if (h_e_t_m(at,zt,nt,a_ftp(t,i))>0)
                        %Where does it go?
                        if (p_e_m(a_ftp(t,i),zt)>0) %Put as manager
                            firm_status(t,i)=2;
                            hire_m(t,i)=1;
                            m_ftp(t,i)=zt;
                            n_ftp(t,i)=0;
                            target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
                            m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
                            n_fwage(t,i)=0.0;
                        else %Put as non manager
                            firm_status(t,i)=3;
                            hire_n(t,i)=1;
                            n_ftp(t,i)=zt;
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
                            n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                        end
                    % Hire the non manager?
                    elseif (h_e_t_nm(at,zt,nt,a_ftp(t,i))>0)
                        %Where does it go?
                        if (p_e_m(a_ftp(t,i),nt)>0) %Put as manager
                            firm_status(t,i)=2;
                            hire_m(t,i)=1;
                            m_ftp(t,i)=nt;
                            n_ftp(t,i)=0;
                            target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vmh(a_ftp(t,i),m_ftp(t,i))- Veh(a_ftp(t,i)));
                            m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
                            n_fwage(t,i)=0.0;
                        else %Put as non manager
                            firm_status(t,i)=3;
                            hire_n(t,i)=1;
                            n_ftp(t,i)=nt;
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vnh(a_ftp(t,i),n_ftp(t,i))- Veh(a_ftp(t,i)));
                            n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                        end
                    else %Firm does not hire
                        firm_status(t,i)=1;
                        m_ftp(t,i)=0;
                        n_ftp(t,i)=0;
                        m_fwage(t,i)=0.0;
                        n_fwage(t,i)=0.0;
                    end
                end
                else %Firm does not match
                      firm_status(t,i)=1;
                      m_ftp(t,i)=0;
                      n_ftp(t,i)=0;
                      m_fwage(t,i)=0.0;
                      n_fwage(t,i)=0.0;
                end
            end 
            
            % fprintf('SM for status 2')
            if (firm_status(t-1,i)==2) %Firm with manager
                %Check productivity shock
                if (r.ashock(t,i)<down_a(a_ftp(t-1,i))) %Firm go down in a
                    a_ftp(t,i)=a_ftp(t-1,i)-1;
                elseif (r.ashock(t,i)< down_a(a_ftp(t-1,i))+stay_a(a_ftp(t-1,i)) && r.ashock(t,i)>=down_a(a_ftp(t-1,i))) %Firm stay in a
                    a_ftp(t,i)=a_ftp(t-1,i);
                else %Firm go up in a
                    a_ftp(t,i)=a_ftp(t-1,i)+1;
                end
                %Update the states do we can use below
                m_ftp(t,i)=m_ftp(t-1,i);
                n_ftp(t,i)=n_ftp(t-1,i);
                m_fwage(t,i)=m_fwage(t-1,i);
                n_fwage(t,i)=n_fwage(t-1,i);
                % Check if worker leaves the firm
                if (r.event(t,i)<=del+death)
                    firm_status(t,i)=1;
                    m_ftp(t,i)=0;
                    n_ftp(t,i)=0;
                    m_fwage(t,i)=0.0;
                    n_fwage(t,i)=0.0;
                elseif (r.event(t,i)> del+death && r.event(t,i)<= del+death+prob_matching) %Match with someone
                    %Match unemp
                    if (r.state(t,i)<=cdf_emp_dist(1))
                        zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
                        % Hire unemp?
                        if (h_m_u(zu,a_ftp(t,i),m_ftp(t-1,i))>0)
                            %Where does it go?
                            if (p_m_m_u(a_ftp(t,i),m_ftp(t-1,i),zu)>0) %Put as manager firing current
                                firm_status(t,i)=2;
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*U(zu)+ bpw*(Vmh(a_ftp(t,i),zu)+U(m_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),zu ),target_wage,wgrid);
                                n_fwage(t,i)=0.0;
                                m_ftp(t,i)=zu;
                                n_ftp(t,i)=0;
                            elseif (p_m_m_d(a_ftp(t,i),m_ftp(t-1,i),zu)>0) %Put as manager with demotion
                                firm_status(t,i)=4; %becomes a team
                                hire_m(t,i)=1;
                                prom_sm(t,i)=-1; %Demote the manager
                                target_wage_m=(1-bpw)*U(m_ftp(t,i))+bpw*(Vth(a_ftp(t,i),zu,m_ftp(t,i))-cost_d - Vmh(a_ftp(t,i),m_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zu,m_ftp(t,i)),target_wage_m,wgrid);
                                n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
                                m_ftp(t,i)=zu;
                                n_ftp(t,i)=m_ftp(t-1,i); %Demote the manager
                            else %Put as non manager
                                firm_status(t,i)=4; %becomes a team
                                hire_n(t,i)=1;
                                m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
                                target_wage=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zu)- Vmh(a_ftp(t,i),m_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zu),target_wage,wgrid);
                                m_ftp(t,i)=m_ftp(t-1,i); %Manager stays
                                n_ftp(t,i)=zu;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=2;
                            m_ftp(t,i)=m_ftp(t-1,i);
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=m_fwage(t-1,i);
                            n_fwage(t,i)=0.0;
                        end
                    %Match employed manager
                    elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<=cdf_emp_dist(2))
                        [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
                        % Hire employed manager?
                        if (h_m_m(am,zm,a_ftp(t,i),m_ftp(t,i))>0)
                            %Where does it go?
                            if (p_m_m_u(a_ftp(t,i),m_ftp(t-1,i),zm)>0) %Put as manager firing current
                                firm_status(t,i)=2;
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vmh(a_ftp(t,i),zm)+U(m_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),zm),target_wage,wgrid);
                                n_fwage(t,i)=0.0;
                                m_ftp(t,i)=zm;
                                n_ftp(t,i)=0;
                            elseif (p_m_m_d(a_ftp(t,i),m_ftp(t-1,i),zm)>0) %Put as manager with demotion
                                firm_status(t,i)=4; %becomes a team
                                hire_m(t,i)=1;
                                prom_sm(t,i)=-1; %Demote the manager
                                target_wage_m=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),zm,m_ftp(t,i))-cost_d - Vmh(a_ftp(t,i),m_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zm,m_ftp(t,i)),target_wage_m,wgrid);
                                n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
                                m_ftp(t,i)=zm;
                                n_ftp(t,i)=m_ftp(t-1,i); %Demote the manager
                            else %Put as non manager
                                firm_status(t,i)=4; %becomes a team
                                hire_n(t,i)=1;
                                m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
                                target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zm)- Vmh(a_ftp(t,i),m_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zm),target_wage,wgrid);
                                m_ftp(t,i)=m_ftp(t-1,i); %Manager stays
                                n_ftp(t,i)=zm;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=2;
                            m_ftp(t,i)=m_ftp(t-1,i);
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=m_fwage(t-1,i);
                            n_fwage(t,i)=0.0;
                        end
                    %Match employed non manager
                    elseif (r.state(t,i)>cdf_emp_dist(2) && r.state(t,i)<=cdf_emp_dist(3))
                        [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
                        % Hire employed non manager?
                        if (h_m_nm(an,zn,a_ftp(t,i),m_ftp(t,i))>0)
                            %Where does it go?
                            if (p_m_m_u(a_ftp(t,i),m_ftp(t-1,i),zn)>0) %Put as manager firing current
                                firm_status(t,i)=2;
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vmh(a_ftp(t,i),zn)+U(m_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),zn),target_wage,wgrid);
                                n_fwage(t,i)=0.0;
                                m_ftp(t,i)=zn;
                                n_ftp(t,i)=0;
                            elseif (p_m_m_d(a_ftp(t,i),m_ftp(t-1,i),zn)>0) %Put as manager with demotion
                                firm_status(t,i)=4; %becomes a team
                                hire_m(t,i)=1;
                                prom_sm(t,i)=-1; %Demote the manager
                                target_wage_m=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),zn,m_ftp(t,i))-cost_d - Vmh(a_ftp(t,i),m_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zn,m_ftp(t,i)),target_wage_m,wgrid);
                                n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
                                m_ftp(t,i)=zn;
                                n_ftp(t,i)=m_ftp(t-1,i); %Demote the manager
                            else %Put as non manager
                                firm_status(t,i)=4; %becomes a team
                                hire_n(t,i)=1;
                                m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
                                target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zn)- Vmh(a_ftp(t,i),m_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zn),target_wage,wgrid);
                                m_ftp(t,i)=m_ftp(t-1,i); %Manager stays
                                n_ftp(t,i)=zn;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=2;
                            m_ftp(t,i)=m_ftp(t-1,i);
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=m_fwage(t-1,i);
                            n_fwage(t,i)=0.0;
                        end
                    %Match with team
                    else 
                        [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
                        % Hire the manager?
                        if (h_m_t_m(at,zt,nt,a_ftp(t,i),m_ftp(t,i))>0)
                            %Where does it go?
                            if (p_m_m_u(a_ftp(t,i),m_ftp(t-1,i),zt)>0) %Put as manager firing current
                                firm_status(t,i)=2;
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vmh(a_ftp(t,i),zt)+U(m_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),zt),target_wage,wgrid);
                                n_fwage(t,i)=0.0;
                                m_ftp(t,i)=zt;
                                n_ftp(t,i)=0;
                            elseif (p_m_m_d(a_ftp(t,i),m_ftp(t-1,i),zt)>0) %Put as manager with demotion
                                firm_status(t,i)=4; %becomes a team
                                hire_m(t,i)=1;
                                prom_sm(t,i)=-1; %Demote the manager
                                target_wage_m=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),zt,m_ftp(t,i))-cost_d - Vmh(a_ftp(t,i),m_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,m_ftp(t,i)),target_wage_m,wgrid);
                                n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
                                m_ftp(t,i)=zt;
                                n_ftp(t,i)=m_ftp(t-1,i); %Demote the manager
                            else %Put as non manager
                                firm_status(t,i)=4; %becomes a team
                                hire_n(t,i)=1;
                                m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
                                target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt)) + bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zt)- Vmh(a_ftp(t,i),m_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zt),target_wage,wgrid);
                                n_ftp(t,i)=zt;
                            end
                        % Hire the non manager?
                        elseif (h_m_t_nm(at,zt,nt,a_ftp(t,i),m_ftp(t,i))>0)
                            %Where does it go?
                            if (p_m_m_u(a_ftp(t,i),m_ftp(t-1,i),nt)>0) %Put as manager firing current
                                firm_status(t,i)=2;
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vmh(a_ftp(t,i),nt)+U(m_ftp(t,i))- Vmh(a_ftp(t,i),m_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),nt),target_wage,wgrid);
                                n_fwage(t,i)=0.0;
                                m_ftp(t,i)=nt;
                                n_ftp(t,i)=0;
                            elseif (p_m_m_d(a_ftp(t,i),m_ftp(t-1,i),nt)>0) %Put as manager with demotion
                                firm_status(t,i)=4; %becomes a team
                                hire_m(t,i)=1;
                                prom_sm(t,i)=-1; %Demote the manager
                                target_wage_m=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),nt,m_ftp(t,i))-cost_d - Vmh(a_ftp(t,i),m_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),nt,m_ftp(t,i)),target_wage_m,wgrid);
                                n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
                                m_ftp(t,i)=nt;
                                n_ftp(t,i)=m_ftp(t-1,i); %Demote the manager
                            else %Put as non manager
                                firm_status(t,i)=4; %becomes a team
                                hire_n(t,i)=1;
                                m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
                                target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),nt)- Vmh(a_ftp(t,i),m_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),nt),target_wage,wgrid);
                                m_ftp(t,i)=m_ftp(t-1,i); %Manager stays
                                n_ftp(t,i)=nt;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=2;
                            m_ftp(t,i)=m_ftp(t-1,i);
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=m_fwage(t-1,i);
                            n_fwage(t,i)=0.0;
                        end
                    end
                elseif (r.event(t,i)>del+death+prob_matching && r.event(t,i)<=del+death+prob_matching+ lam) % Firm is met by another
                    if (r.event(t,i)<cdf_firm_dist(1)) %Met By vacant firm
                        av=draw_CDF_1d(cdf_e_edist,r.type(t,i));
                        %Does the poacher hire?
                        if (h_e_m(a_ftp(t,i),m_ftp(t,i),av)>0)
                            firm_status(t,i)=1;
                            m_ftp(t,i)=0;
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=0.0;
                        else %Firm does not hire
                            firm_status(t,i)=2;
                            %Check if wage goes up
                            target_wage=max(Vmh(av,m_ftp(t,i)),Vnh(av,m_ftp(t,i)))-Veh(av); %Marginal of the poaching firm
                            current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wmh(:,a_ftp(t,i),m_ftp(t,i)));
                            if (target_wage>current_wage)
                                m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
                            else
                                m_fwage(t,i)=m_fwage(t-1,i);
                            end
                            n_ftp(t,i)=0;
                            n_fwage(t,i)=0.0;
                        end
                    elseif (r.event(t,i)>=cdf_firm_dist(1) && r.event(t,i)<cdf_firm_dist(2)) %Met by firm with manager
                        [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
                        %Does the poacher hire?
                        if (h_m_m(a_ftp(t,i),m_ftp(t,i),am,zm)>0)
                            firm_status(t,i)=1;
                            m_ftp(t,i)=0;
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=0.0;
                        else %Firm does not hire
                            firm_status(t,i)=2;
                            %Check if wage goes up
                            nu=[Vmh(am,m_ftp(t,i))+U(zm),Vth(am,zm,m_ftp(t,i)), Vth(am,m_ftp(t,i),zm)-cost_d];
                            target_wage_m=max(nu)-Vmh(am,zm); %Marginal of the poaching firm
                            current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wmh(:,a_ftp(t,i),m_ftp(t,i)));
                            if (target_wage_m>current_wage)
                                m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage_m,wgrid);
                            else
                                m_fwage(t,i)=m_fwage(t-1,i);
                            end
                            n_ftp(t,i)=0;
                            n_fwage(t,i)=0.0;
                        end
                    elseif (r.event(t,i)>=cdf_firm_dist(2) && r.event(t,i)<cdf_firm_dist(3)) %Met by firm with non manager
                        [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
                        %Does the poacher hire?
                        if (h_nm_m(a_ftp(t,i),m_ftp(t,i),an,zn)>0)
                            firm_status(t,i)=1;
                            m_ftp(t,i)=0;
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=0.0;
                        else %Firm does not hire
                            firm_status(t,i)=2;
                            %Check if wage goes up
                            nu=[Vnh(an,m_ftp(t,i))+U(zn),Vth(an,zn,m_ftp(t,i))-cost_p, Vth(an,m_ftp(t,i),zn)];
                            target_wage_m=max(nu)-Vnh(an,zn); %Marginal of the poaching firm
                            current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wmh(:,a_ftp(t,i),m_ftp(t,i)));
                            if (target_wage_m>current_wage)
                                m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage_m,wgrid);
                            else
                                m_fwage(t,i)=m_fwage(t-1,i);
                            end
                            n_ftp(t,i)=0;
                            n_fwage(t,i)=0.0;
                        end
                    else %Met by team
                        [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
                        %Does the poacher hire?
                        if (h_t_m(a_ftp(t,i),m_ftp(t,i),at,zt,nt)>0)
                            firm_status(t,i)=1;
                            m_ftp(t,i)=0;
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=0.0;
                        else %Firm does not hire
                            firm_status(t,i)=2;
                            %Check if wage goes up
                            nu=[Vth(at,zt,m_ftp(t,i))+U(nt),Vth(at,m_ftp(t,i),zt)-cost_d+U(nt), Vth(at,nt,m_ftp(t,i))+U(zt)-cost_p, Vth(at,m_ftp(t,i),nt)+U(zt)];
                            target_wage_m=max(nu)-Vth(at,zt,nt) %Marginal of the poaching firm
                            current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wmh(:,a_ftp(t,i),m_ftp(t,i)));
                            if (target_wage_m>current_wage)
                                m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage_m,wgrid);
                            else
                                m_fwage(t,i)=m_fwage(t-1,i);
                            end
                            n_ftp(t,i)=0;
                            n_fwage(t,i)=0.0;
                        end
                    end
                else %Nothing happens to the firm at SM
                      firm_status(t,i)=2;
                      m_ftp(t,i)=m_ftp(t-1,i);
                      n_ftp(t,i)=0;
                      m_fwage(t,i)=m_fwage(t-1,i);
                      n_fwage(t,i)=0.0;
                end
                % %Reallocation post SM
                % if (d_m(a_ftp(t,i),m_ftp(t,i))>0) %Fire the manager
                %     firm_status(t,i)=1;
                %     m_ftp(t,i)=0;
                %     n_ftp(t,i)=0;
                %     m_fwage(t,i)=0.0;
                %     n_fwage(t,i)=0.0;
                % elseif (r_m(a_ftp(t,i),m_ftp(t,i))>0) %reallocate the manager
                %     firm_status(t,i)=3;
                %     n_ftp(t,i)=m_ftp(t,i);
                %     n_fwage(t,i)=m_fwage(t,i); %Wage does not change upon reallocation
                %     m_ftp(t,i)=0;
                %     m_fwage(t,i)=0.0;
                % else %Manager stays
                %     firm_status(t,i)=2;
                %     %Does the wage goes down?
                %     current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wmh(:,a_ftp(t,i),m_ftp(t,i)));
                %     target_wage=Vmh(a_ftp(t,i),m_ftp(t,i))-Veh(a_ftp(t,i));
                %     if (target_wage<current_wage)
                %         m_fwage(t,i)= InterpolateWage(Wmh(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
                %     else
                %         m_fwage(t,i)=m_fwage(t-1,i);
                %     end
                %     n_ftp(t,i)=0;
                %     n_fwage(t,i)=0.0;
                % end
            end
            
            % fprintf('SM for status 3')
            if (firm_status(t-1,i)==3) %Firm with non manager
                %Check productivity shock
                if (r.ashock(t,i)<down_a(a_ftp(t-1,i))) %Firm go down in a
                    a_ftp(t,i)=a_ftp(t-1,i)-1;
                elseif (r.ashock(t,i)< down_a(a_ftp(t-1,i))+stay_a(a_ftp(t-1,i)) && r.ashock(t,i)>=down_a(a_ftp(t-1,i))) %Firm stay in a
                    a_ftp(t,i)=a_ftp(t-1,i);
                else %Firm go up in a
                    a_ftp(t,i)=a_ftp(t-1,i)+1;
                end
                %Check Learning shock on the non manager
                if (r.qshock(t,i)<down_q(a_ftp(t-1,i), n_ftp(t-1,i))) %Firm go down in q
                    n_ftp(t,i)=n_ftp(t-1,i)-1;    
                elseif (r.qshock(t,i)< down_q(a_ftp(t-1,i), n_ftp(t-1,i))+stay_q(a_ftp(t-1,i), n_ftp(t-1,i)) && r.qshock(t,i)>=down_q(a_ftp(t-1,i), n_ftp(t-1,i))) %Firm stay in q
                    n_ftp(t,i)=n_ftp(t-1,i);
                else %Firm go up in q
                    n_ftp(t,i)=n_ftp(t-1,i)+1;
                end
                %Update the states do we can use below
                m_ftp(t,i)=m_ftp(t-1,i);
                m_fwage(t,i)=m_fwage(t-1,i);
                n_fwage(t,i)=n_fwage(t-1,i);
                % Check if worker leaves the firm
                if (r.event(t,i)<=del+death)
                    firm_status(t,i)=1;
                    m_ftp(t,i)=0;
                    n_ftp(t,i)=0;
                    m_fwage(t,i)=0.0;
                    n_fwage(t,i)=0.0;
                elseif (r.event(t,i)> del+death && r.event(t,i)<= del+death+prob_matching) %Match with someone
                    %Match unemp
                    if (r.state(t,i)<=cdf_emp_dist(1)) zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
                        % Hire unemp?
                        if (h_nm_u(zu,a_ftp(t,i),n_ftp(t,i))>0)
                            %Where does it go?
                            if (p_n_n_u(a_ftp(t,i),n_ftp(t,i),zu)>0) %Put as non manager firing current
                                firm_status(t,i)=3;
                                hire_n(t,i)=1;
                                target_wage=(1-bpw)*U(zu)+bpw*(Vnh(a_ftp(t,i),zu)+U(n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),zu),target_wage,wgrid);
                                m_fwage(t,i)=0.0;
                                m_ftp(t,i)=0;
                                n_ftp(t,i)=zu;
                            elseif (p_n_n_p(a_ftp(t,i),n_ftp(t,i),zu)>0) %Put as non manager with promotion
                                firm_status(t,i)=4; %becomes a team
                                hire_n(t,i)=1;
                                prom_sm(t,i)=1; %Promote the non manager 
                                target_wage_n=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zu)-cost_p - Vnh(a_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zu),target_wage_n,wgrid);
                                m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
                                m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
                                n_ftp(t,i)=zu;
                            else %Put as manager
                                firm_status(t,i)=4; %becomes a team
                                hire_m(t,i)=1;
                                n_fwage(t,i)=n_fwage(t-1,i); %Wage of the non manager do not change now
                                target_wage=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),zu,n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zu,n_ftp(t,i)),target_wage,wgrid);
                                m_ftp(t,i)=zu;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=3;
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=n_fwage(t-1,i);
                        end
                    %Match employed manager
                    elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<=cdf_emp_dist(2))
                        %Match employed manager
                        [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
                        % Hire employed manager?
                        if (h_nm_m(am,zm,a_ftp(t,i),n_ftp(t,i))>0)
                            %Where does it go?
                            if (p_n_n_u(a_ftp(t,i),n_ftp(t,i),zm)>0) %Put as non manager firing current
                                firm_status(t,i)=3;
                                hire_n(t,i)=1;
                                target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vnh(a_ftp(t,i),zm)+U(n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),zm),target_wage,wgrid);
                                m_fwage(t,i)=0.0;
                                m_ftp(t,i)=0;
                                n_ftp(t,i)=zm;
                            elseif (p_n_n_p(a_ftp(t,i),n_ftp(t,i),zm)>0) %Put as non manager with promotion
                                firm_status(t,i)=4; %becomes a team
                                hire_n(t,i)=1;
                                prom_sm(t,i)=1; %Promote the non manager
                                target_wage_n=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zm)-cost_p - Vnh(a_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zm),target_wage_n,wgrid);
                                m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
                                m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
                                n_ftp(t,i)=zm;
                            else %Put as manager
                                firm_status(t,i)=4; %becomes a team
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),zm,n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zm,n_ftp(t,i)),target_wage,wgrid);
                                n_fwage(t,i)=n_fwage(t-1,i); %Wage of the non manager do not change now
                                m_ftp(t,i)=zm;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=3;
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=n_fwage(t-1,i);
                        end
                    %Match employed non manager
                    elseif (r.state(t,i)>cdf_emp_dist(2) && r.state(t,i)<=cdf_emp_dist(3))
                        [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
                        % Hire employed non manager?
                        if (h_nm_nm(an,zn,a_ftp(t,i),n_ftp(t,i))>0)
                            %Where does it go?
                            if (p_n_n_u(a_ftp(t,i),n_ftp(t,i),zn)>0) %Put as non manager firing current
                                firm_status(t,i)=3;
                                hire_n(t,i)=1;
                                target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vnh(a_ftp(t,i),zn)+U(n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),zn),target_wage,wgrid);
                                m_fwage(t,i)=0.0;
                                m_ftp(t,i)=0;
                                n_ftp(t,i)=zn;
                            elseif (p_n_n_p(a_ftp(t,i),n_ftp(t,i),zn)>0) %Put as non manager with promotion
                                firm_status(t,i)=4; %becomes a team
                                hire_n(t,i)=1;
                                prom_sm(t,i)=1; %Promote the non manager
                                target_wage_n=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zn)-cost_p - Vnh(a_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zn),target_wage_n,wgrid);
                                m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
                                m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
                                n_ftp(t,i)=zn;
                            else %Put as manager
                                firm_status(t,i)=4; %becomes a team
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),zn,n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zn,n_ftp(t,i)),target_wage,wgrid);
                                n_fwage(t,i)=n_fwage(t-1,i); %Wage of the non
                                m_ftp(t,i)=zn;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=3;
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=n_fwage(t-1,i);
                        end
                    %Match with team
                    else 
                        [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
                        % Hire the manager?
                        if (h_nm_t_m(at,zt,nt,a_ftp(t,i),n_ftp(t,i))>0)
                            %Where does it go?
                            if (p_n_n_u(a_ftp(t,i),n_ftp(t,i),zt)>0) %Put as non manager firing current
                                firm_status(t,i)=3;
                                hire_n(t,i)=1;
                                target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vnh(a_ftp(t,i),zt)+U(n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),zt),target_wage,wgrid);
                                m_fwage(t,i)=0.0;
                                m_ftp(t,i)=0;
                                n_ftp(t,i)=zt;
                            elseif (p_n_n_p(a_ftp(t,i),n_ftp(t,i),zt)>0) %Put as non manager with promotion
                                firm_status(t,i)=4; %becomes a team
                                hire_n(t,i)=1;
                                prom_sm(t,i)=1; %Promote the non manager
                                target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zt)-cost_p - Vnh(a_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zt),target_wage_n,wgrid);
                                m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
                                m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
                                n_ftp(t,i)=zt;
                            else %Put as manager
                                firm_status(t,i)=4; %becomes a team
                                hire_m(t,i)=1;
                                n_fwage(t,i)=n_fwage(t-1,i); %Wage of the non manager do not change now
                                target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),zt,n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,n_ftp(t,i)),target_wage,wgrid);
                                m_ftp(t,i)=zt;
                            end
                        % Hire the non manager?
                        elseif (h_nm_t_nm(at,zt,nt,a_ftp(t,i),n_ftp(t,i))>0)
                            %Where does it go?
                            if (p_n_n_u(a_ftp(t,i),n_ftp(t,i),nt)>0) %Put as non manager firing current
                                firm_status(t,i)=3;
                                hire_n(t,i)=1;
                                target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vnh(a_ftp(t,i),nt)+U(n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),nt),target_wage,wgrid);
                                m_fwage(t,i)=0.0;
                                m_ftp(t,i)=0;
                                n_ftp(t,i)=nt;
                            elseif (p_n_n_p(a_ftp(t,i),n_ftp(t,i),nt)>0) %Put as non manager with promotion
                                firm_status(t,i)=4; %becomes a team
                                hire_n(t,i)=1;
                                prom_sm(t,i)=1; %Promote the non manager
                                target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),nt)-cost_p - Vnh(a_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),nt),target_wage_n,wgrid);
                                m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
                                m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
                                n_ftp(t,i)=nt;
                            else %Put as manager
                                firm_status(t,i)=4; %becomes a team
                                hire_m(t,i)=1;
                                n_fwage(t,i)=n_fwage(t-1,i); %Wage of the non manager do not change now
                                target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),nt,n_ftp(t,i))- Vnh(a_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),nt,n_ftp(t,i)),target_wage,wgrid);
                                m_ftp(t,i)=nt;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=3;
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=n_fwage(t-1,i);
                        end
                    end
                elseif (r.event(t,i)>del+death+prob_matching && r.event(t,i)<=del+death+prob_matching+ lam) % Firm is met by another
                    if (r.event(t,i)<cdf_firm_dist(1)) %Met By vacant firm
                        av=draw_CDF_1d(cdf_e_edist,r.type(t,i));
                        %Does the poacher hire?
                        if (h_e_nm(a_ftp(t,i),n_ftp(t,i),av)>0)
                            firm_status(t,i)=1;
                            m_ftp(t,i)=0;
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=0.0;
                        else %Firm does not hire
                            firm_status(t,i)=3;
                            %Check if wage goes up
                            target_wage=max(Vnh(av,n_ftp(t,i)),Vmh(av,n_ftp(t,i)))-Veh(av); %Marginal of the poaching firm
                            current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wnh(:,a_ftp(t,i),n_ftp(t,i)));
                            if (target_wage>current_wage)
                                n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                            else
                                n_fwage(t,i)=n_fwage(t-1,i);
                            end
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                        end
                    elseif (r.event(t,i)>=cdf_firm_dist(1) && r.event(t,i)<cdf_firm_dist(2)) %Met by firm with manager
                        [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
                        %Does the poacher hire?
                        if (h_m_nm(a_ftp(t,i),n_ftp(t,i),am,zm)>0)
                            firm_status(t,i)=1;
                            m_ftp(t,i)=0;
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=0.0;
                        else %Firm does not hire
                            firm_status(t,i)=3;
                            %Check if wage goes up
                            nu=[Vmh(am,n_ftp(t,i))+U(zm),Vth(am,zm,n_ftp(t,i)), Vth(am,n_ftp(t,i),zm)-cost_d];
                            target_wage_m=max(nu)-Vmh(am,zm); %Marginal of the poaching firm
                            current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wnh(:,a_ftp(t,i),n_ftp(t,i)));
                            if (target_wage_m>current_wage)
                                n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage_m,wgrid);
                            else
                                n_fwage(t,i)=n_fwage(t-1,i);
                            end
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                        end
                    elseif (r.event(t,i)>=cdf_firm_dist(2) && r.event(t,i)<cdf_firm_dist(3)) %Met by firm with non manager
                        [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
                        %Does the poacher hire?
                        if (h_nm_nm(a_ftp(t,i),n_ftp(t,i),an,zn)>0)
                            firm_status(t,i)=1;
                            m_ftp(t,i)=0;
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=0.0;
                        else %Firm does not hire
                            firm_status(t,i)=3;
                            %Check if wage goes up
                            nu=[Vnh(an,n_ftp(t,i))+U(zn),Vth(an,zn,n_ftp(t,i))-cost_p, Vth(an,n_ftp(t,i),zn)];
                            target_wage_m=max(nu)-Vnh(an,zn); %Marginal of the poaching firm
                            current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wnh(:,a_ftp(t,i),n_ftp(t,i)));
                            if (target_wage_m>current_wage)
                                n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage_m,wgrid);
                            else
                                n_fwage(t,i)=n_fwage(t-1,i);
                            end
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                        end
                    else %Met by team
                        [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
                        %Does the poacher hire?
                        if (h_t_nm(a_ftp(t,i),n_ftp(t,i),at,zt,nt)>0)
                            firm_status(t,i)=1;
                            m_ftp(t,i)=0;
                            n_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=0.0;
                        else %Firm does not hire
                            firm_status(t,i)=3;
                            %Check if wage goes up
                            nu=[Vth(at,zt,n_ftp(t,i))+U(nt),Vth(at,n_ftp(t,i),zt)-cost_d+U(nt), Vth(at,nt,n_ftp(t,i))+U(zt)-cost_p, Vth(at,n_ftp(t,i),nt)+U(zt)];
                            target_wage_m=max(nu)-Vth(at,zt,nt) %Marginal of the poaching firm
                            current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wnh(:,a_ftp(t,i),n_ftp(t,i)));
                            if (target_wage_m>current_wage)
                                n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage_m,wgrid);
                            else
                                n_fwage(t,i)=n_fwage(t-1,i);
                            end
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                        end
                    end
                else %Nothing happens to the firm at SM
                      firm_status(t,i)=3;
                      m_ftp(t,i)=0;
                      n_ftp(t,i)=n_ftp(t-1,i);
                      m_fwage(t,i)=0.0;
                      n_fwage(t,i)=n_fwage(t-1,i);
                end
                % %Reallocation post SM
                % if (d_n(a_ftp(t,i),n_ftp(t,i))>0) %Fire the non manager
                %     firm_status(t,i)=1;
                %     m_ftp(t,i)=0;
                %     n_ftp(t,i)=0;
                %     m_fwage(t,i)=0.0;
                %     n_fwage(t,i)=0.0;
                % elseif (r_n(a_ftp(t,i),n_ftp(t,i))>0) %reallocate the non manager
                %     firm_status(t,i)=2; %Firm with manager
                %     m_ftp(t,i)=n_ftp(t,i);
                %     m_fwage(t,i)=n_fwage(t,i); %Wage does not change upon reallocation
                %     n_ftp(t,i)=0;
                %     n_fwage(t,i)=0.0;
                % else %Non manager stays
                %     firm_status(t,i)=3;
                %     %Does the wage goes down?
                %     current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wnh(:,a_ftp(t,i),n_ftp(t,i)));
                %     target_wage=Vnh(a_ftp(t,i),n_ftp(t,i))-Veh(a_ftp(t,i));
                %     if (target_wage<current_wage)
                %         n_fwage(t,i)= InterpolateWage(Wnh(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                %     else
                %         n_fwage(t,i)=n_fwage(t-1,i);
                %     end
                %     m_ftp(t,i)=0;
                %     m_fwage(t,i)=0.0;
                % end
            end
                            
            % fprintf('SM for status 4')
            if (firm_status(t-1,i)==4) %Firm with team
                %Check productivity shock
                if (r.ashock(t,i)<down_a(a_ftp(t-1,i))) %Firm go down in a
                    a_ftp(t,i)=a_ftp(t-1,i)-1;
                elseif (r.ashock(t,i)< down_a(a_ftp(t-1,i))+stay_a(a_ftp(t-1,i)) && r.ashock(t,i)>=down_a(a_ftp(t-1,i))) %Firm stay in a
                    a_ftp(t,i)=a_ftp(t-1,i);
                else %Firm go up in a
                    a_ftp(t,i)=a_ftp(t-1,i)+1;
                end
                %Check Learning shock on the non manager
                if (r.qshock(t,i)<down_q(a_ftp(t-1,i), n_ftp(t-1,i))) %Firm go down in q
                    n_ftp(t,i)=n_ftp(t-1,i)-1;
                elseif (r.qshock(t,i)< down_q(a_ftp(t-1,i), n_ftp(t-1,i)) + stay_q(a_ftp(t-1,i),n_ftp(t-1,i)) && r.qshock(t,i)>=down_q(a_ftp(t-1,i), n_ftp(t-1,i))) %Firm stay in q
                    n_ftp(t,i)=n_ftp(t-1,i);
                else %Firm go up in q
                    n_ftp(t,i)=n_ftp(t-1,i)+1;
                end
                %Update the states do we can use below
                m_ftp(t,i)=m_ftp(t-1,i);
                m_fwage(t,i)=m_fwage(t-1,i);
                n_fwage(t,i)=n_fwage(t-1,i);
                % Check if manager leaves the firm
                if (r.event(t,i)<=del+death)
                    firm_status(t,i)=3;
                    m_ftp(t,i)=0;
                    n_ftp(t,i)=n_ftp(t-1,i);
                    m_fwage(t,i)=0.0;
                    n_fwage(t,i)=n_fwage(t-1,i);
                %check if non manager leaves the firm
                elseif (r.event(t,i)> del+death && r.event(t,i)<= del+death+del+death) 
                    firm_status(t,i)=2;
                    m_ftp(t,i)=m_ftp(t-1,i);
                    n_ftp(t,i)=0;
                    m_fwage(t,i)=m_fwage(t-1,i);
                    n_fwage(t,i)=0.0;
                elseif (r.event(t,i)> del+death+del+death && r.event(t,i)<= del+death+del+death+prob_matching) %Match with someone
                    %Match unemp
                    if (r.state(t,i)<=cdf_emp_dist(1)) 
                        zu=draw_CDF_1d(cdf_e_udist,r.type(t,i));
                        % Hire unemp?
                        if (h_t_u(zu,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0)
                            %Where does it go?
                            if (p_t_m_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zu)>0) %Put as manager firing current
                                firm_status(t,i)=4;
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),zu,n_ftp(t,i))+U(m_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zu,n_ftp(t,i)),target_wage,wgrid);
                                n_fwage(t,i)=n_fwage(t-1,i);
                                m_ftp(t,i)=zu;
                            elseif (p_t_m_d(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zu)>0) %Put as manager with demotion
                                firm_status(t,i)=4; %still a team
                                hire_m(t,i)=1;
                                prom_sm(t,i)=-1; %Demote the manager
                                target_wage_m=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),zu,m_ftp(t,i))+U(n_ftp(t,i))-cost_d- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zu,m_ftp(t,i)),target_wage_m,wgrid);
                                n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
                                n_ftp(t,i)=m_ftp(t,i); %Demote the manager
                                m_ftp(t,i)=zu;
                            elseif (p_t_n_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zu)>0) %Put as non manager firing current
                                firm_status(t,i)=4; %still a team
                                hire_n(t,i)=1;
                                target_wage_n=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zu)+U(n_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zu),target_wage_n,wgrid);
                                m_fwage(t,i)=m_fwage(t-1,i); %Wage of the non manager do not change now
                                n_ftp(t,i)=zu;
                            else %Put as non manager with promotion
                                firm_status(t,i)=4; %still a team
                                hire_n(t,i)=1;
                                prom_sm(t,i)=1; %Promote the non manager
                                target_wage_n=(1-bpw)*U(zu)+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zu)-cost_p+ U(m_ftp(t,i)) -Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zu),target_wage_n,wgrid);
                                m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
                                m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
                                n_ftp(t,i)=zu;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=4;
                            m_ftp(t,i)=m_ftp(t-1,i);
                            m_fwage(t,i)=m_fwage(t-1,i);
                            n_fwage(t,i)=n_fwage(t-1,i);
                        end
                    %Match employed manager
                    elseif (r.state(t,i)>=cdf_emp_dist(1) && r.state(t,i)<=cdf_emp_dist(2))
                        %Match employed manager
                        [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
                        % Hire employed manager?
                        if (h_t_m(am,zm,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0)
                            %Where does it go?
                            if (p_t_m_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zm)>0) %Put as manager firing current
                                firm_status(t,i)=4;
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),zm,n_ftp(t,i))+U(m_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zm,n_ftp(t,i)),target_wage,wgrid);
                                n_fwage(t,i)=n_fwage(t-1,i);
                                m_ftp(t,i)=zm;
                            elseif (p_t_m_d(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zm)>0) %Put as manager with demotion
                                firm_status(t,i)=4; %still a team
                                hire_m(t,i)=1;
                                prom_sm(t,i)=-1; %Demote the manager
                                target_wage_m=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),zm,m_ftp(t,i))+U(n_ftp(t,i))-cost_d- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zm,m_ftp(t,i)),target_wage_m,wgrid);
                                n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
                                n_ftp(t,i)=m_ftp(t,i); %Demote the manager
                                m_ftp(t,i)=zm;
                            elseif (p_t_n_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zm)>0) %Put as non manager firing current
                                firm_status(t,i)=4; %still a team
                                hire_n(t,i)=1;
                                target_wage_n=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zm)+U(n_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zm),target_wage_n,wgrid);
                                m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
                                n_ftp(t,i)=zm;
                            else %Put as non manager with promotion
                                firm_status(t,i)=4; %still a team
                                hire_n(t,i)=1;
                                prom_sm(t,i)=1; %Promote the non manager
                                target_wage_n=(1-bpw)*(Vmh(am,zm)-Veh(am))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zm)-cost_p + U(m_ftp(t,i))-Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zm),target_wage_n,wgrid);
                                m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
                                m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
                                n_ftp(t,i)=zm;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=4;
                            m_ftp(t,i)=m_ftp(t-1,i);
                            m_fwage(t,i)=m_fwage(t-1,i);
                            n_fwage(t,i)=n_fwage(t-1,i);
                        end
                    %Match employed non manager
                    elseif (r.state(t,i)>cdf_emp_dist(2) && r.state(t,i)<=cdf_emp_dist(3))
                        [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
                        % Hire employed non manager?
                        if (h_t_nm(an,zn,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0)
                            %Where does it go?
                            if (p_t_m_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zn)>0) %Put as manager firing current
                                firm_status(t,i)=4;
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),zn,n_ftp(t,i))+U(m_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zn,n_ftp(t,i)),target_wage,wgrid);
                                n_fwage(t,i)=n_fwage(t-1,i);
                                m_ftp(t,i)=zn;
                            elseif (p_t_m_d(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zn)>0) %Put as manager with demotion
                                firm_status(t,i)=4; %still a team
                                hire_m(t,i)=1;
                                prom_sm(t,i)=-1; %Demote the manager
                                target_wage_m=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),zn,m_ftp(t,i))+U(n_ftp(t,i))-cost_d- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zn,m_ftp(t,i)),target_wage_m,wgrid);
                                n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
                                n_ftp(t,i)=m_ftp(t,i); %Demote the manager
                                m_ftp(t,i)=zn;
                            elseif (p_t_n_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zn)>0) %Put as non manager firing current
                                firm_status(t,i)=4; %still a team
                                hire_n(t,i)=1;
                                target_wage_n=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zn)+U(n_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zn),target_wage_n,wgrid);
                                m_fwage(t,i)=m_fwage(t-1,i); %Wage of the non manager do not change now
                                n_ftp(t,i)=zn;
                            else %Put as non manager with promotion
                                firm_status(t,i)=4; %still a team
                                hire_n(t,i)=1;
                                prom_sm(t,i)=1; %Promote the non manager
                                target_wage_n=(1-bpw)*(Vnh(an,zn)-Veh(an))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zn)-cost_p + U(m_ftp(t,i))-Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zn),target_wage_n,wgrid);
                                m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
                                m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
                                n_ftp(t,i)=zn;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=4;
                            m_ftp(t,i)=m_ftp(t-1,i);
                            m_fwage(t,i)=m_fwage(t-1,i);
                            n_fwage(t,i)=n_fwage(t-1,i);
                        end
                    %Match with team
                    else 
                        [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
                        % Hire the manager?
                        if (h_t_t_m(at,zt,nt,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0)
                            %Where does it go?
                            if (p_t_m_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zt)>0) %Put as manager firing current
                                firm_status(t,i)=4;
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),zt,n_ftp(t,i))+U(m_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,n_ftp(t,i)),target_wage,wgrid);
                                n_fwage(t,i)=n_fwage(t-1,i);
                                m_ftp(t,i)=zt;
                            elseif (p_t_m_d(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zt)>0) %Put as manager with demotion
                                firm_status(t,i)=4; %still a team
                                hire_m(t,i)=1;
                                prom_sm(t,i)=-1; %Demote the manager
                                target_wage_m=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),zt,m_ftp(t,i))+U(n_ftp(t,i))-cost_d- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,m_ftp(t,i)),target_wage_m,wgrid);
                                n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
                                n_ftp(t,i)=m_ftp(t,i); %Demote the manager
                                m_ftp(t,i)=zt;
                            elseif (p_t_n_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),zt)>0) %Put as non manager firing current
                                firm_status(t,i)=4; %still a team
                                hire_n(t,i)=1;
                                target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zt)+U(n_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zt),target_wage_n,wgrid);
                                m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
                                n_ftp(t,i)=zt;
                            else %Put as non manager with promotion
                                firm_status(t,i)=4; %still a team
                                hire_n(t,i)=1;
                                prom_sm(t,i)=1; %Promote the non manager
                                target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vnh(at,nt))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zt)-cost_p + U(m_ftp(t,i))-Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zt),target_wage_n,wgrid);
                                m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
                                m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
                                n_ftp(t,i)=zt;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=4;
                            m_ftp(t,i)=m_ftp(t-1,i);
                            m_fwage(t,i)=m_fwage(t-1,i);
                            n_fwage(t,i)=n_fwage(t-1,i);
                        end
            
                        % Hire the non manager?
                        if (h_t_t_nm(at,zt,nt,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0)
                            %Where does it go?
                            if (p_t_m_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),nt)>0) %Put as manager firing current
                                firm_status(t,i)=4;
                                hire_m(t,i)=1;
                                target_wage=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),zt,n_ftp(t,i))+U(m_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,n_ftp(t,i)),target_wage,wgrid);
                                n_fwage(t,i)=n_fwage(t-1,i);
                                m_ftp(t,i)=zt;
                            elseif (p_t_m_d(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),nt)>0) %Put as manager with demotion
                                firm_status(t,i)=4; %still a team
                                hire_m(t,i)=1;
                                prom_sm(t,i)=-1; %Demote the manager
                                target_wage_m=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),zt,m_ftp(t,i))+U(n_ftp(t,i))-cost_d- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),zt,m_ftp(t,i)),target_wage_m,wgrid);
                                n_fwage(t,i)=m_fwage(t-1,i); %Wage of the demoted manager do not change now
                                n_ftp(t,i)=m_ftp(t,i); %Demote the manager
                                m_ftp(t,i)=zt;
                            elseif (p_t_n_u(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),nt)>0) %Put as non manager firing current
                                firm_status(t,i)=4; %still a team
                                hire_n(t,i)=1;
                                target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),m_ftp(t,i),zt)+U(n_ftp(t,i))- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),zt),target_wage_n,wgrid);
                                m_fwage(t,i)=m_fwage(t-1,i); %Wage of the manager do not change now
                                n_ftp(t,i)=zt;
                            else %Put as non manager with promotion
                                firm_status(t,i)=4; %still a team
                                hire_n(t,i)=1;
                                prom_sm(t,i)=1; %Promote the non manager
                                target_wage_n=(1-bpw)*(Vth(at,zt,nt)-Vmh(at,zt))+bpw*(Vth(a_ftp(t,i),n_ftp(t,i),zt)-cost_p + U(m_ftp(t,i))-Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),n_ftp(t,i),zt),target_wage_n,wgrid);
                                m_fwage(t,i)=n_fwage(t-1,i); %Wage of the promoted non manager do not change now
                                m_ftp(t,i)=n_ftp(t,i); %Promote the non manager
                                n_ftp(t,i)=zt;
                            end
                        else %Firm does not hire
                            firm_status(t,i)=4;
                            m_ftp(t,i)=m_ftp(t-1,i);
                            m_fwage(t,i)=m_fwage(t-1,i);
                            n_fwage(t,i)=n_fwage(t-1,i);
                        end
                    end
                elseif (r.event(t,i)>del+death+prob_matching && r.event(t,i)<=del+death+prob_matching+ lam) % Firm is met by another
                    if (r.event(t,i)<cdf_firm_dist(1)) %Met By vacant firm
                        av=draw_CDF_1d(cdf_e_edist,r.type(t,i));
                        %Does the poacher hire the manager?
                        if (h_e_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),av)>0)
                            firm_status(t,i)=3; %becomes a firm with non manager
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=n_fwage(t-1,i);
                        %Does the poacher hire the non manager?    
                        elseif (h_e_t_nm(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),av)>0)
                            firm_status(t,i)=2; %becomes a firm with manager
                            n_ftp(t,i)=0;
                            n_fwage(t,i)=0.0;
                            m_fwage(t,i)=m_fwage(t-1,i);
                        else %Firm does not hire
                            firm_status(t,i)=4;
                            %Check if wage goes up
                            %Who had the biggest gain from trade? This will determine which, if any, wage goes up
                            gain_m=max(Vmh(av,m_ftp(t,i)),Vnh(av,m_ftp(t,i)))-Veh(av) - Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vnh(a_ftp(t,i),n_ftp(t,i));
                            gain_n=max(Vmh(av,n_ftp(t,i)),Vnh(av,n_ftp(t,i)))-Veh(av) - Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vmh(a_ftp(t,i),m_ftp(t,i));
                            if (gain_m>=gain_n) %Manager might get a raise
                                target_wage=max(Vmh(av,m_ftp(t,i)),Vnh(av,m_ftp(t,i)))-Veh(av); %Marginal of the poaching firm
                                current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                if (target_wage>current_wage)
                                    m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                                else
                                    m_fwage(t,i)=m_fwage(t-1,i);
                                end
                                n_fwage(t,i)=n_fwage(t-1,i);
                            else %Non manager might get a raise
                                target_wage=max(Vmh(av,n_ftp(t,i)),Vnh(av,n_ftp(t,i)))-Veh(av); %Marginal of the poaching firm
                                current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                if (target_wage>current_wage)
                                    n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                                else
                                    n_fwage(t,i)=n_fwage(t-1,i);
                                end
                                m_fwage(t,i)=m_fwage(t-1,i);
                            end
                        end
                    elseif (r.event(t,i)>=cdf_firm_dist(1) && r.event(t,i)<cdf_firm_dist(2)) %Met by firm with manager
                        [am,zm]=draw_CDF_2d(cdf_e_mdist,ats,tpts,r.type(t,i));
                        %Does the poacher hire the manager?
                        if (h_m_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),am,zm)>0)
                            firm_status(t,i)=3; %becomes a firm with non manager
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=n_fwage(t-1,i);
                        %Does the poacher hire the non manager?    
                        elseif (h_m_t_nm(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),am,zm)>0)
                            firm_status(t,i)=2; %becomes a firm with manager
                            n_ftp(t,i)=0;
                            n_fwage(t,i)=0.0;
                            m_fwage(t,i)=m_fwage(t-1,i);
                        else %Firm does not hire
                            firm_status(t,i)=4;
                            %Check if wage goes up
                            %Who had the biggest gain from trade? This will determine which, if any, wage goes up
                            nu_m=[Vmh(am,m_ftp(t,i))+U(zm),Vth(am,zm,m_ftp(t,i)), Vth(am,m_ftp(t,i),zm)-cost_d];
                            nu_n=[Vmh(am,n_ftp(t,i))+U(zm),Vth(am,zm,n_ftp(t,i)), Vth(am,n_ftp(t,i),zm)];
                            gain_m=max(nu_m)-Vmh(am,zm)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vnh(a_ftp(t,i),n_ftp(t,i));
                            gain_n=max(nu_n)-Vnh(am,zm)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vmh(a_ftp(t,i),m_ftp(t,i));
                            if (gain_m>=gain_n) %Manager might get a raise
                                target_wage=max(nu_m)-Vmh(am,zm); %Marginal of the poaching firm
                                current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                if (target_wage>current_wage)
                                    m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                                else
                                    m_fwage(t,i)=m_fwage(t-1,i);
                                end
                                n_fwage(t,i)=n_fwage(t-1,i);
                            else %Non manager might get a raise
                                target_wage=max(nu_n)-Vnh(am,zm); %Marginal of the poaching firm
                                current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                if (target_wage>current_wage)
                                    n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                                else
                                    n_fwage(t,i)=n_fwage(t-1,i);
                                end
                                m_fwage(t,i)=m_fwage(t-1,i);
                            end
                        end
                    elseif (r.event(t,i)>=cdf_firm_dist(2) && r.event(t,i)<cdf_firm_dist(3)) %Met by firm with non manager
                        [an,zn]=draw_CDF_2d(cdf_e_ndist,ats,tpts,r.type(t,i));
                        %Does the poacher hire the manager?
                        if (h_nm_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),an,zn)>0)
                            firm_status(t,i)=3; %becomes a firm with non manager
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=n_fwage(t-1,i);
                        %Does the poacher hire the non manager?    
                        elseif (h_nm_t_nm(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),an,zn)>0)
                            firm_status(t,i)=2; %becomes a firm with manager
                            n_ftp(t,i)=0;
                            n_fwage(t,i)=0.0;
                            m_fwage(t,i)=m_fwage(t-1,i);
                        else %Firm does not hire
                            firm_status(t,i)=4;
                            %Check if wage goes up
                            %Who had the biggest gain from trade? This will determine which, if any, wage goes up
                            nu_m=[Vnh(an,m_ftp(t,i))+U(zn),Vth(an,zn,m_ftp(t,i)), Vth(an,m_ftp(t,i),zn)-cost_d];
                            nu_n=[Vnh(an,n_ftp(t,i))+U(zn),Vth(an,zn,n_ftp(t,i))-cost_p, Vth(an,n_ftp(t,i),zn)];
                            gain_m=max(nu_m)-Vnh(an,zn)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vnh(a_ftp(t,i),n_ftp(t,i));
                            gain_n=max(nu_n)-Vnh(an,zn)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vmh(a_ftp(t,i),m_ftp(t,i));
                            if (gain_m>=gain_n) %Manager might get a raise
                                target_wage=max(nu_m)-Vnh(an,zn); %Marginal of the poaching firm
                                current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                if (target_wage>current_wage)
                                    m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                                else
                                    m_fwage(t,i)=m_fwage(t-1,i);
                                end
                                n_fwage(t,i)=n_fwage(t-1,i);
                            else %Non manager might get a raise
                                target_wage=max(nu_n)-Vnh(an,zn); %Marginal of the poaching firm
                                current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                if (target_wage>current_wage)
                                    n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                                else
                                    n_fwage(t,i)=n_fwage(t-1,i);
                                end
                                m_fwage(t,i)=m_fwage(t-1,i);
                            end
                        end
                    else %Met by firm with team
                        [at,zt,nt]=draw_CDF_3d(cdf_e_tdist,ats,tpts,r.type(t,i));
                        %Does the poacher hire the manager?
                        if (h_t_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),at,zt,nt)>0)
                            firm_status(t,i)=3; %becomes a firm with non manager
                            m_ftp(t,i)=0;
                            m_fwage(t,i)=0.0;
                            n_fwage(t,i)=n_fwage(t-1,i);
                        %Does the poacher hire the non manager?
                        elseif (h_t_t_nm(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i),at,zt,nt)>0)
                            firm_status(t,i)=2; %becomes a firm with manager
                            n_ftp(t,i)=0;
                            n_fwage(t,i)=0.0;
                            m_fwage(t,i)=m_fwage(t-1,i);
                        else %Firm does not hire
                            firm_status(t,i)=4;
                            %Check if wage goes up
                            %Who had the biggest gain from trade? This will determine which, if any, wage goes up
                            nu_m=[Vth(at,zt,m_ftp(t,i))+U(nt),Vth(at,m_ftp(t,i),zt)-cost_d+U(nt), Vth(at,m_ftp(t,i),nt)+U(zt), Vth(at,nt,m_ftp(t,i))+U(zt)-cost_p];
                            nu_n=[Vth(at,zt,n_ftp(t,i))+U(nt),Vth(at,n_ftp(t,i),zt)-cost_d+U(nt), Vth(at,n_ftp(t,i),nt)+U(zt), Vth(at,nt,n_ftp(t,i))+U(zt)-cost_p];
                            gain_m=max(nu_m)-Vth(at,zt,nt)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vnh(a_ftp(t,i),n_ftp(t,i));
                            gain_n=max(nu_n)-Vth(at,zt,nt)- Vth(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))+Vmh(a_ftp(t,i),m_ftp(t,i));
                            if (gain_m>=gain_n) %Manager might get a raise
                                target_wage=max(nu_m)-Vth(at,zt,nt); %Marginal of the poaching firm
                                current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                if (target_wage>current_wage)
                                    m_fwage(t,i)= InterpolateWage(Wtmh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                                else
                                    m_fwage(t,i)=m_fwage(t-1,i);
                                end
                                n_fwage(t,i)=n_fwage(t-1,i);
                            else %Non manager might get a raise
                                target_wage=max(nu_n)-Vth(at,zt,nt); %Marginal of the poaching firm
                                current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                                if (target_wage>current_wage)
                                    n_fwage(t,i)= InterpolateWage(Wtnh(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                                else
                                    n_fwage(t,i)=n_fwage(t-1,i);
                                end
                                m_fwage(t,i)=m_fwage(t-1,i);
                            end
                        end
                    end
                else %Nothing happens at SM
                    firm_status(t,i)=4;
                    m_ftp(t,i)=m_ftp(t-1,i);
                    m_fwage(t,i)=m_fwage(t-1,i);
                    n_fwage(t,i)=n_fwage(t-1,i);
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Reallocation post SM
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % fprintf('R for status 2')
            if (firm_status(t,i)==2) %Gets to the end of SM with a manager
                if (d_m(a_ftp(t,i),m_ftp(t,i))>0)  %Fire the manager
                    firm_status(t,i)=1;
                    fire_m(t,i)=1;
                    m_ftp(t,i)=0;
                    n_ftp(t,i)=0;
                    m_fwage(t,i)=0.0;
                    n_fwage(t,i)=0.0;
                elseif (r_m(a_ftp(t,i),m_ftp(t,i))>0) %reallocate the manager
                    firm_status(t,i)=3;
                    prom_sm(t,i)=-1; %Demote the manager
                    n_ftp(t,i)=m_ftp(t,i);
                    n_fwage(t,i)=m_fwage(t,i); %Wage does not change upon reallocation
                    m_ftp(t,i)=0;
                    m_fwage(t,i)=0.0;
                else %Manager stays
                    firm_status(t,i)=2;
                    %Does the wage goes down?
                    current_wage=InterpolateWage(wgrid, m_fwage(t,i),Wm(:,a_ftp(t,i),m_ftp(t,i)));
                    target_wage=Vm(a_ftp(t,i),m_ftp(t,i))-Ve(a_ftp(t,i));
                    if (target_wage<current_wage)
                        m_fwage(t,i)= InterpolateWage(Wm(:,a_ftp(t,i),m_ftp(t,i)),target_wage,wgrid);
                    end
                    n_ftp(t,i)=0;
                    n_fwage(t,i)=0.0;
                end
            end
    
            % fprintf('R for status 3')
            if (firm_status(t,i)==3) %Gets to the end of SM with a non manager
                %Reallocation post SM
                if (d_n(a_ftp(t,i),n_ftp(t,i))>0) %Fire the non manager
                    firm_status(t,i)=1;
                    fire_n(t,i)=1;
                    m_ftp(t,i)=0;
                    n_ftp(t,i)=0;
                    m_fwage(t,i)=0.0;
                    n_fwage(t,i)=0.0;
                elseif (r_n(a_ftp(t,i),n_ftp(t,i))>0) %reallocate the non manager
                    firm_status(t,i)=2; %Firm with manager
                    prom_sm(t,i)=1; %Promote the non manager
                    m_ftp(t,i)=n_ftp(t,i);
                    m_fwage(t,i)=n_fwage(t,i); %Wage does not change upon reallocation
                    n_ftp(t,i)=0;
                    n_fwage(t,i)=0.0;
                else %Non manager stays
                    firm_status(t,i)=3;
                    %Does the wage goes down?
                    current_wage=InterpolateWage(wgrid, n_fwage(t,i),Wn(:,a_ftp(t,i),n_ftp(t,i)));
                    target_wage=Vn(a_ftp(t,i),n_ftp(t,i))-Ve(a_ftp(t,i));
                    if (target_wage<current_wage)
                        n_fwage(t,i)= InterpolateWage(Wn(:,a_ftp(t,i),n_ftp(t,i)),target_wage,wgrid);
                    end
                    m_ftp(t,i)=0;
                    m_fwage(t,i)=0.0;
                end
            end
            
            % fprintf('R for status 4')
            if (firm_status(t,i)==4) % Gets to the end of Sm with a team 
                if (d_t_b(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0) %Fire the team
                    firm_status(t,i)=1;
                    fire_m(t,i)=1;
                    fire_n(t,i)=1;
                    m_ftp(t,i)=0;
                    n_ftp(t,i)=0;
                    m_fwage(t,i)=0.0;
                    n_fwage(t,i)=0.0;
                elseif (r_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))*r_t_n(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0) %reallocate the team
                    %Lets be careful here to ensure the swap
                    firm_status(t,i)=4;
                    prom_realloc(t,i)=1; %Promote the non manager
                    temp_m=m_ftp(t,i);
                    temp_n=n_ftp(t,i);
                    temp_mw=m_fwage(t,i);
                    temp_nw=n_fwage(t,i);
                    m_ftp(t,i)=temp_n;
                    n_ftp(t,i)=temp_m;
                    m_fwage(t,i)=temp_nw;
                    n_fwage(t,i)=temp_mw;
                    clear temp_m temp_n temp_mw temp_nw;
                elseif (d_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))*(1-r_t_n(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)))>0) %Fire the manager do not reallocate nm
                    firm_status(t,i)=3;
                    fire_m(t,i)=1;
                    m_ftp(t,i)=0;
                    m_fwage(t,i)=0.0;
                elseif (d_t_n(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))*(1-r_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)))>0) %Fire the non manager do not reallocate m
                    firm_status(t,i)=2;
                    fire_n(t,i)=1;
                    n_ftp(t,i)=0;
                    n_fwage(t,i)=0.0;
                elseif (d_t_n(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))*r_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0) %Fire the non manager reallocate m
                    firm_status(t,i)=3;
                    fire_n(t,i)=1;
                    prom_realloc(t,i)=-1; %Demote the manager
                    n_ftp(t,i)=m_ftp(t,i);
                    n_fwage(t,i)=m_fwage(t,i); %Wage does not change upon reallocation
                    m_ftp(t,i)=0;
                    m_fwage(t,i)=0.0;
                elseif (d_t_m(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))*r_t_n(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))>0) %Fire the manager reallocate nm
                    firm_status(t,i)=2;
                    fire_m(t,i)=1;
                    prom_realloc(t,i)=1; %Promote the non manager
                    m_ftp(t,i)=n_ftp(t,i);
                    m_fwage(t,i)=n_fwage(t,i); %Wage does not change upon reallocation
                    n_ftp(t,i)=0;
                    n_fwage(t,i)=0.0;
                else %Team stays
                    %Does the wage goes down?
                    %For the manager
                    current_wage_m=InterpolateWage(wgrid, m_fwage(t,i),Wtm(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                    target_wage_m=Vt(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))-Vn(a_ftp(t,i),n_ftp(t,i));
                    if (target_wage_m<current_wage_m)
                        m_fwage(t,i)= InterpolateWage(Wtm(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage_m,wgrid);
                    end
                    %For the non manager
                    current_wage_n=InterpolateWage(wgrid, n_fwage(t,i),Wtn(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)));
                    target_wage_n=Vt(a_ftp(t,i),m_ftp(t,i),n_ftp(t,i))-Vm(a_ftp(t,i),m_ftp(t,i));
                    if (target_wage_n<current_wage_n)
                        n_fwage(t,i)= InterpolateWage(Wtn(:,a_ftp(t,i),m_ftp(t,i),n_ftp(t,i)),target_wage_n,wgrid);
                    end
                end
            end
        end
    end

    fs=struct('firm_status',firm_status,'a_ftp',a_ftp,'m_ftp',m_ftp,'n_ftp',n_ftp,'m_fwage',m_fwage,'n_fwage',n_fwage,'hire_m',hire_m,'hire_n',hire_n,'fire_m',fire_m,'fire_n',fire_n,'prom_sm',prom_sm,'prom_realloc',prom_realloc);

    

    
    


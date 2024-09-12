function model_moments= model_moments(ps,sp,fs,ws,moments_selection)
    %Unpack sp, fs and ws, ps 
    fieldNames = fieldnames(ps);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with
        % the same
        eval([varName ' = ps.' varName ';']);
    end
    
    
    fieldNames = fieldnames(sp);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with
        % the same
        eval([varName ' = sp.' varName ';']);
    end
    
    fieldNames = fieldnames(fs);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = fs.' varName ';']);
    end
    
    fieldNames = fieldnames(ws);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = ws.' varName ';']);
    end
    
    
    %% Data Wide Format

    % Data for firms 
    
    prod=fs.a_ftp;
    status=fs.firm_status;
    man_type=fs.m_ftp;
    nman_type=fs.n_ftp;
    man_wage=fs.m_fwage;
    nman_wage=fs.n_fwage;
    man_hire=fs.hire_m;
    nman_hire=fs.hire_n;
    man_fire=fs.fire_m;
    nman_fire=fs.fire_n;
    prom_sm=fs.prom_sm;
    prom_realloc=fs.prom_realloc;
    
    %12 variables for the firm path
    prod(1:t_burn,:)=[]; %Burn the first t_burn months
    status(1:t_burn,:)=[]; 
    man_type(1:t_burn,:)=[];
    nman_type(1:t_burn,:)=[];
    man_wage(1:t_burn,:)=[];
    nman_wage(1:t_burn,:)=[];
    man_hire(1:t_burn,:)=[];
    nman_hire(1:t_burn,:)=[];
    man_fire(1:t_burn,:)=[];
    nman_fire(1:t_burn,:)=[];
    prom_sm(1:t_burn,:)=[];
    prom_realloc(1:t_burn,:)=[];
    
    %Year variable size of n_months-t_burn and change every 12 months
    year=zeros(n_months-t_burn,1);
    for i=1:n_months-t_burn
        year(i)=ceil(i/12);
    end
    
    %"Monthly" data for now
    data=zeros(n_firms*(n_months-t_burn),15);
    data_row=0;
    id=0;
    
    for j=1:n_firms
        period=0;
        id=id+1;
        for t=1:(n_months-t_burn)
            period=period+1;
            data_row=data_row+1;
            data(data_row,:)=[id,period,year(t),prod(t,j),status(t,j),man_type(t,j),nman_type(t,j),man_wage(t,j),nman_wage(t,j),man_hire(t,j),nman_hire(t,j), man_fire(t,j),nman_fire(t,j),prom_sm(t,j),prom_realloc(t,j)];
        end
    end
    
    I_remove=find(sum(data==0,2)==28);
    data(I_remove,:)=[];
            
    %Upack data
    id = data(:,1);
    period = data(:,2);
    year = data(:,3);
    prod = data(:,4);
    status = data(:,5);
    man_type = data(:,6);
    nman_type = data(:,7);
    man_wage = data(:,8);
    nman_wage =data(:,9);
    man_hire = data(:,10);
    nman_hire = data(:,11);
    man_fire = data(:,12);
    nman_fire = data(:,13);
    prom_sm = data(:,14);
    prom_realloc = data(:,15);
    
    clearvars data

    %data for workers
    w_status=ws.worker_status;
    wage=ws.self_wage;
    type=ws.self_z;
    tenure=ws.self_tenure;
    cw_wage=ws.co_wage;
    cw_tenure=ws.co_tenure;
    cw_type=ws.co_z;
    firm_type=ws.a_wtp;
    W_hire_manager=ws.hire_m;
    W_hire_nonmanager=ws.hire_n;
    W_fire_manager=ws.fire_m;
    W_fire_nonmanager=ws.fire_n;
    W_promote_sm=ws.prom_sm;
    W_promote_realloc=ws.prom_realloc;

    %14 variables for the worker path
    w_status(1:t_burn,:)=[]; %Burn the first t_burn months
    wage(1:t_burn,:)=[];
    type(1:t_burn,:)=[];
    tenure(1:t_burn,:)=[];
    cw_wage(1:t_burn,:)=[];
    cw_tenure(1:t_burn,:)=[];
    cw_type(1:t_burn,:)=[];
    firm_type(1:t_burn,:)=[];
    W_hire_manager(1:t_burn,:)=[];
    W_hire_nonmanager(1:t_burn,:)=[];
    W_fire_manager(1:t_burn,:)=[];
    W_fire_nonmanager(1:t_burn,:)=[];
    W_promote_sm(1:t_burn,:)=[];
    W_promote_realloc(1:t_burn,:)=[];


    %Year variable size of n_months-t_burn and change every 12 months
    year=zeros(n_months-t_burn,1);
    for i=1:n_months-t_burn
        year(i)=ceil(i/12);
    end

    %"Monthly" data for now
    dataw=zeros(n_workers*(n_months-t_burn),17);
    data_row=0;
    idw=0;

    for j=1:n_workers
        periodw=0;
        idw=idw+1;
        for t=1:(n_months-t_burn)
            periodw=periodw+1;
            data_row=data_row+1;
            dataw(data_row,:)=[idw,periodw,year(t),w_status(t,j),wage(t,j),type(t,j),tenure(t,j),cw_wage(t,j),cw_tenure(t,j),cw_type(t,j),firm_type(t,j),W_hire_manager(t,j),W_hire_nonmanager(t,j), W_fire_manager(t,j),W_fire_nonmanager(t,j),W_promote_sm(t,j),W_promote_realloc(t,j)];
        end
    end


    I_remove=find(sum(dataw==0,2)==28);
    dataw(I_remove,:)=[];

    %Unpack data
    idw = dataw(:,1);
    periodw = dataw(:,2);
    year = dataw(:,3);
    w_status = dataw(:,4);
    wage = dataw(:,5);
    type = dataw(:,6);
    tenure = dataw(:,7);
    cw_wage = dataw(:,8);
    cw_tenure = dataw(:,9);
    cw_type = dataw(:,10);
    firm_type = dataw(:,11);
    W_hire_manager = dataw(:,12);
    W_hire_nonmanager = dataw(:,13);
    W_fire_manager = dataw(:,14);
    W_fire_nonmanager = dataw(:,15);
    W_promote_sm = dataw(:,16);
    W_promote_realloc = dataw(:,17);

    
    clearvars dataw

    %Dummy Manager 
    dummy_manager=zeros(n_workers*(n_months-t_burn),1);
    for i=1:n_workers*(n_months-t_burn)
        if w_status(i)==2 || w_status(i)==4
            dummy_manager(i)=1;
        end
    end

    %Dummy for firm type (excluding the lowest level of productivity)
    dummy_firm2=zeros(n_workers*(n_months-t_burn),1);
    for i=1:n_workers*(n_months-t_burn)
        if firm_type(i)==2
            dummy_firm2(i)=1;
        end
    end

    dummy_firm3=zeros(n_workers*(n_months-t_burn),1);
    for i=1:n_workers*(n_months-t_burn)
        if firm_type(i)==3
            dummy_firm3(i)=1;
        end
    end

    dummy_firm4=zeros(n_workers*(n_months-t_burn),1);
    for i=1:n_workers*(n_months-t_burn)
        if firm_type(i)==4
            dummy_firm4(i)=1;
        end
    end

    dummy_firm5=zeros(n_workers*(n_months-t_burn),1);
    for i=1:n_workers*(n_months-t_burn)
        if firm_type(i)==5
            dummy_firm5(i)=1;
        end
    end

            

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Model moments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %From firm data

    [Iq1]=find(prod==1);
    [Iq2]=find(prod==2);
    I_man=find(status==2 | status==4); %manager alone or with team
    I_nman=find(status==3 | status==4); %non manager
    %Average wage of managers in firms at q1
    avg_q1=sum(man_wage(Iq1))/length(Iq1);
    avg_q2=sum(man_wage(Iq2))/length(Iq2);
    
    %Promotion
    I_prom=find(prom_sm==1 | prom_realloc==1);
    rate_prom=length(I_prom)/length(I_man);
    
    
    %Model moments
    %Moments
    %Promotion
    I_prom=find(prom_sm==1 | prom_realloc==1);
    mm.NiMi = length(I_prom)/length(I_man);
    
    %Demotion
    I_dem=find(prom_sm==-1 | prom_realloc==-1);
    mm.MiNi = length(I_dem)/length(I_nman);
    
    % Split the sample by firms with each level of productivity
    
    % Iq5=find(prod==5);
    nm_avg_q=zeros(1,ats);
    for a=1:ats
        Iq{a}=find(prod==a);
        nm_avg_q(a)=sum(nman_wage(Iq{a}))/length(Iq{a});
    end


    % Wage ratios by quintile of firm productivity
    if ats== 5
        mm.Q5Q1nm_wage = nm_avg_q(5)/nm_avg_q(1);                    
        mm.Q5Q2nm_wage = nm_avg_q(5)/nm_avg_q(2);                   
        mm.Q5Q3nm_wage = nm_avg_q(5)/nm_avg_q(3);                     
        mm.Q5Q4nm_wage = nm_avg_q(5)/nm_avg_q(4);                     
        mm.Q5Q5nm_wage = nm_avg_q(5)/nm_avg_q(5);
    else 
        mm.Q5Q1nm_wage =  nm_avg_q(2)/nm_avg_q(1);
        mm.Q5Q5nm_wage= nm_avg_q(2)/nm_avg_q(2);
    end


    % Manager to worker wage ratio
    I_team=find(status==4); %team
    ratio= man_wage(I_team)./nman_wage(I_team);
    mm.ManWorkerRatio = mean(ratio);

    %From worker data


    %Crosse section reg inside the model 
    %Sample of employed workers
    I_samp=find(w_status~=1 & wage>0);
    length(I_samp ); 
   
    %Winsorize variables
    wins_lb_cut=.5;
    wins_ub_cut=99.5;

    %Winsorize wages
    % lw=winsorize(log(wage(I_samp)),wins_lb_cut,wins_ub_cut);
    %Manually winsorize
    lw=log(wage(I_samp));
    lw(lw<prctile(lw,wins_lb_cut))=prctile(lw,wins_lb_cut);
    lw(lw>prctile(lw,wins_ub_cut))=prctile(lw,wins_ub_cut);
    

    if ats==5
        X=[ones(length(I_samp),1),dummy_manager(I_samp),dummy_firm2(I_samp),dummy_firm3(I_samp),dummy_firm4(I_samp),dummy_firm5(I_samp),tenure(I_samp)];
    else
        X=[ones(length(I_samp),1),dummy_manager(I_samp),dummy_firm2(I_samp),tenure(I_samp)];
    end
    [~,nX]=size(X);
    Y=lw;
    b_cross=(X'*X)\X'*Y;   
    eps=Y-X*b_cross;
    sig_sq=sum(eps.^2)/length(I_samp);
    std_error_cross=(sig_sq*eye(nX)/(X'*X))^.5;

    mm.b_cross=b_cross(2);

    % Loop over the fields specified in moments_selection
    model_moments = zeros(1, length(moments_selection));
    for i = 1:length(moments_selection)
        field_name = moments_selection{i};  % Get the field name from the selection
        model_moments(i) = mm.(field_name);  % Access the corresponding field in the struct
    end
    




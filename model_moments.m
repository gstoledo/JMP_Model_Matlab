%From test script long_data

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
        
    
    % See if understand what is happening
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

    % Loop over the fields specified in moments_selection
    model_moments = zeros(1, length(moments_selection));
    for i = 1:length(moments_selection)
        field_name = moments_selection{i};  % Get the field name from the selection
        model_moments(i) = mm.(field_name);  % Access the corresponding field in the struct
    end
    


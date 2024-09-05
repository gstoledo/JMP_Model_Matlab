function [v,e,w]=joint_loop(p,ps,tg,spec_name)
    %% Opening up the parameters
    fieldNames = fieldnames(p);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = p.' varName ';']);
    end
    %Opening up the pre set parameters
    fieldNames = fieldnames(ps);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = ps. ' varName ';']);
    end
    %Opening up the toggles
    fieldNames = fieldnames(tg);
    % Loop over each field and assign it to a variable in the workspace
    for i = 1:length(fieldNames)% Dynamically create the variable name
        varName = fieldNames{i};
        % Use eval to assign the value to the variable with the same name
        eval([varName ' = tg.' varName ';']);
    end

    %% Derivated from the main parameters
    [a_trans,q_trans,u_trans,fteam,fman,fnman,fe,b,mnew_high,typebirth,wmin,wmax,wgrid] = derivatated_p(p,ps);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model_file = 'solved_model.mat';
    if exist(model_file, 'file')==2 && use_guess(1)=='y'
        load(model_file, 'v', 'e', 'w');
        %Making sure the initial guesses are the same size as the model and they are not empty
        if length(v.Ve)==ats && length(v.U)==tpts && ~isempty(v.Ve) && ~isempty(v.U) && ~isempty(v.Vm) && ~isempty(v.Vn) && ~isempty(v.Vt)
            Veini=v.Ve;
            Vmini=v.Vm;
            Vnini=v.Vn;
            Vtini=v.Vt;
            Uini=v.U;
            eplus_udist=e.eplus_udist;
            eplus_edist=e.eplus_edist;
            eplus_mdist=e.eplus_mdist;
            eplus_ndist=e.eplus_ndist;
            eplus_tdist=e.eplus_tdist;
        else
            clear v e w    
            %%Initial guesses for the LOM
            eplus_udist=(1/tpts)*ones(1,tpts);
            eplus_edist=(n/ats)*ones(1,ats);    
            eplus_mdist=zeros(ats,tpts); %Distribution of firms with manager e_m(a,z)
            eplus_ndist=zeros(ats,tpts); %Distribution of firms with no manager e_n(a,q)
            eplus_tdist=zeros(ats,tpts,tpts); %Distribution of firms with team e_t(a,z,q)

            %Inital guesses for value functionsz
            Veini=zeros(1,ats);
            Vmini=fman/(1-bt);
            Vnini=fnman/(1-bt);
            Vtini=fteam/(1-bt);
            Uini=b/(1-bt);
        end 
    else        
        %%Initial guesses for the LOM
        eplus_udist=(1/tpts)*ones(1,tpts);
        eplus_edist=(n/ats)*ones(1,ats);    
        eplus_mdist=zeros(ats,tpts); %Distribution of firms with manager e_m(a,z)
        eplus_ndist=zeros(ats,tpts); %Distribution of firms with no manager e_n(a,q)
        eplus_tdist=zeros(ats,tpts,tpts); %Distribution of firms with team e_t(a,z,q)

        %Inital guesses for value functionsz
        Veini=zeros(1,ats);
        Vmini=fman/(1-bt);
        Vnini=fnman/(1-bt);
        Vtini=fteam/(1-bt);
        Uini=b/(1-bt);
    end

    %Initialize the loop with a guess for value functions
    [Ve, Vm, Vn, Vt, U, Veh, Vmh, Vnh, Vth, Vetl, Vmtl, Vntl, Vttl, Utl]= vf_iterationV2(eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist,eplus_udist,Veini,Vmini,Vnini,Vtini,Uini,ats,tpts,cost_d,cost_p,tr,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);


    diff_joint_lag1=0;
    diff_joint_lag2=0;

    %Iterate on value functions and type distribution
    diff_joint    =100;
    diff_joint_max=1e-8; %Value func/distribution max tolerance
    
    it_joint      =0;
    it_joint_min  =250;
    it_joint_max  =2000;
    %% Joint loop
    %Start timer
    % tic;
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
        [Veplus, Vmplus, Vnplus, Vtplus, Uplus, Vehplus, Vmhplus, Vnhplus, Vthplus, Vetlplus, Vmtlplus, Vntlplus, Vttlplus, Utlplus] = vf_iterationV2(e_edist,e_mdist,e_ndist,e_tdist,e_udist,Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p,tr,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed);
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

        % %Print every 100 iterations
        % if mod(it_joint,10)==0
        %     fprintf(string(spec_name)+' Joint Iteration %d, error %f \n', it_joint, diff_joint);
        % end
        if diff_joint<diff_joint_max && diff_joint_lag1<diff_joint_max  && diff_joint_lag2<diff_joint_max && it_joint>=it_joint_min
            fprintf(string(spec_name)+' Joint Converged in %d iterations\n',it_joint);
        end
        if it_joint==it_joint_max
            fprintf(2,'Joint Failed to converge\n')
        end
    end
    %Check sum of the eplus final distributions
    nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
    popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");

    % toc;

    % %% Wage iteration
    % %Inital guesses for wages
    % model_file = 'solution_wages.mat';
    % if exist(model_file, 'file')==2 && use_guess(1)=='y'
    %     load(model_file, 'Wm', 'Wn', 'Wtm', 'Wtn');
    %     if size(Wm,1)==wpts && size(Wm,2)==ats && size(Wm,3)==tpts 
    %         Wmini=Wm;
    %         Wnini=Wn;
    %         Wtmini=Wtm;
    %         Wtnini=Wtn;
    %     else
    %         Wmini=ones(wpts,ats,tpts)*wmin;
    %         Wnini=ones(wpts,ats,tpts)*wmin;
    %         Wtmini=ones(wpts,ats,tpts,tpts)*wmin;
    %         Wtnini=ones(wpts,ats,tpts,tpts)*wmin;
    %     end
    % else
    %     Wmini=ones(wpts,ats,tpts)*wmin;
    %     Wnini=ones(wpts,ats,tpts)*wmin;
    %     Wtmini=ones(wpts,ats,tpts,tpts)*wmin;
    %     Wtnini=ones(wpts,ats,tpts,tpts)*wmin;
    % end

    if exist(model_file, 'file')==2 && use_guess(1)=='y'
        load(model_file, 'w');
        if size(w.Wm,1)==wpts && size(w.Wm,2)==ats && size(w.Wm,3)==tpts && ~isempty(w.Wm) && ~isempty(w.Wn) && ~isempty(w.Wtm) && ~isempty(w.Wtn)
            Wmini=w.Wm;
            Wnini=w.Wn;
            Wtmini=w.Wtm;
            Wtnini=w.Wtn;
        else
            clear w
            Wmini=ones(wpts,ats,tpts)*wmin;
            Wnini=ones(wpts,ats,tpts)*wmin;
            Wtmini=ones(wpts,ats,tpts,tpts)*wmin;
            Wtnini=ones(wpts,ats,tpts,tpts)*wmin;
        end
    else
        Wmini=ones(wpts,ats,tpts)*wmin;
        Wnini=ones(wpts,ats,tpts)*wmin;
        Wtmini=ones(wpts,ats,tpts,tpts)*wmin;
        Wtnini=ones(wpts,ats,tpts,tpts)*wmin;
    end



    % tic
    [Wm,Wn,Wtm,Wtn,Wmh,Wnh,Wtnh,Wtmh]=wf_iteration(wpts,ats,tpts,Ve,Vm,Vn,Vt,U,Vmh,Vnh,Vth,Veh,Wmini,Wnini,Wtmini,Wtnini,...
    speed,cost_d,cost_p,wgrid,bt,death,del,lam,lamu,bpw,n,a_trans,q_trans,e_udist,e_edist,e_mdist,e_ndist,e_tdist);
    % toc
    %Output the results
    v=struct('Ve',Ve,'Vm',Vm,'Vn',Vn,'Vt',Vt,'U',U,'Veh',Veh,'Vmh',Vmh,'Vnh',Vnh,'Vth',Vth);
    e=struct('eplus_udist',eplus_udist,'eplus_edist',eplus_edist,'eplus_mdist',eplus_mdist,'eplus_ndist',eplus_ndist,'eplus_tdist',eplus_tdist,'nplus',nplus,'popplus',popplus);
    w=struct('Wm',Wm,'Wn',Wn,'Wtm',Wtm,'Wtn',Wtn,'Wmh',Wmh,'Wnh',Wnh,'Wtmh',Wtmh,'Wtnh',Wtnh);




    
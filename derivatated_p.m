%Refere to this function to change prod fucntion and trantisiton matrices
function [a_trans,q_trans,u_trans,fteam,fman,fnman,fe,b,mnew_high,typebirth,wmin,wmax,wgrid] = derivatated_p(p,ps)
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


    %% Derivated from the main parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Process for a
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    astay=1-aup-adown;                             %Probability of staying
    a_trans=create_trans(adown,astay,aup,ats); %Transition matrix for productivity

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Process for q (Later we need one value for each a)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qstay=1-qup-qdown;                        %Probability of staying
    q_trans=zeros(tpts,tpts,ats);         %Transition matrix for q

    for a=1:ats
        q_trans(:,:,a)=create_trans(qdown(a),qstay(a),qup(a),tpts);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Process for u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ugain=0.00;                                           %Probability of moving up
    ustay=1-ulose-ugain;                              %Probability of staying
    u_trans=create_trans(ulose,ustay,ugain,tpts); %Transition matrix for unemployment 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Production
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    typemin=1;                                                   %Lowest type
    typemax=5 ;                                                  %Highest type
    amin=1;                                                      %Lowest productivity
    amax=10;                                                     %Highest productivity
    type=[typemin:(typemax-typemin)/(tpts-1):typemax];%Types
    a_type=[amin:(amax-amin)/(ats-1):amax];           %Productivity
    fteam=zeros(ats,tpts,tpts);                         %Production of team
    fman=zeros(ats,tpts);                                  %Production of manager
    fnman=zeros(ats,tpts);                                 %Production of non manager
    fe=zeros(1,ats);                                          %Production of unemployed
    for a=1:ats
        for z=1:tpts
            for q=1:tpts
                fteam(a,z,q)=A*((alpha_m*a_type(a)*type(z)^fcomp + (1-alpha_m)*type(q)^fcomp))^(1/fcomp);
            end
            fman(a,z)=(1/2)*A*a_type(a)*(alpha_m*type(z)^fcomp+ (1-alpha_m)*type(1)^fcomp)^(1/fcomp);
            fnman(a,z)=0;
        end
        fe(a)=0;
    end
    b=homeprod*fman(1,:) ;                                    %Type home production vector

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rates and Trabnsitions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mnew_high=tpts;                                            %Number of types
    typebirth=zeros(1,tpts);                                   %Distributions of birth types
    for i=1:mnew_high
        typebirth(i)= mnew*exp(-type(i)*mnew)./sum(mnew*exp(-type(1:mnew_high)*mnew));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Wages 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wmin  =  -6*min(fteam(:));                                 %Minimum wage in the grid
    wmax = max(fteam(:));                                      %Maximum wage in the grid
    wgrid = [wmin:(wmax - wmin)/(wpts - 1):wmax]; %Wage grid


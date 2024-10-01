

load('xmin_results_combined.mat')
x=xmin; %Use some values of HLPM to start
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Toggles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tg.use_guess='no';                %Use guess for value functions
    tg.zero_tol=1e-10;                  %Tolerance for zero
    tg.nm_penal=0;                            %NManager penalty toggle
    tg.speed=1;                        %Convergence speed
    tg.speed_dist=1;                         %Convergenge speed distribution
    tg.update_speed_v=1;                %Update speed for value functions
    tg.update_speed=1;                  %Update speed for wages
    tg.display_iter=0;                      %Display iterations in joint loop
    tg.display_iter_v=0;                  %Display value function iterations
    tg.display_iter_dist=0;              %Display distribution iterations
    tg.split_top_bot=1;                 %Split top and bottom of distribution for the qup
    tg.annual=0;                        %Annual calibration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Simulation parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sp.seed=1;                          %Seed for random number generator
    sp.n_firms=20000;                      %Number of firms
    sp.n_workers=sp.n_firms;             %Number of workers
    sp.n_months=240;                     %Number of months
    sp.burn_years=1;                  %Number of years to burn
    sp.t_burn=12*sp.burn_years;        %Number of months to burn
    sp.n_years=sp.n_months/12;        %Number of years total
    sp.stable_years = sp.n_years - sp.burn_years; %Number of years to keep the simulation


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Data Moments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if tg.annual==1
        dm.NiMi                             =0.047+0.009;  %Rate internal promotions to managers (as a fraction of total managers)
        dm.MiNi                             =0.021+0.001;  %Rate internal demotions to non-managers (as a fraction of total non-managers)
    else
        dm.NiMi                             = (0.047+0.009)/12;  %Rate internal promotions to managers (as a fraction of total managers)
        dm.MiNi                             = (0.021+0.001)/12;  %Rate internal demotions to non-managers (as a fraction of total non-managers)
    end
    dm.Q5Q1nm_wage                      =6.84;          %Non Manager Wage Ratio quintile 5 to quintile 1 of total avg wage
    dm.Q5Q2nm_wage                      =3.35;          %Non Manager Wage Ratio quintile 5 to quintile 2 of total avg wage
    dm.Q5Q3nm_wage                      =2.30;          %Non Manager Wage Ratio quintile 5 to quintile 3 of total avg wage
    dm.Q5Q4nm_wage                      =1.68;          %Non Manager Wage Ratio quintile 5 to quintile 4 of total avg wage
    dm.Q5Q5nm_wage                      =1;           %Non Manager Wage Ratio quintile 5 to quintile 4 of total avg wage
    dm.ManWorkerRatio                   =1.53;           % Ratio between managers and workers wage/employment ratio
    dm.bcross                           =0.137;         %Cross Section Wage Premium

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Pre-set parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta.tpts = 3;                        %Type points space
    theta.ats = 2;                         %Productivity space
    if location == "hpc"
        theta.ats=5;
        theta.tpts=7;
    end
    theta.wpts=75;                       %Wage type
    if tg.annual==1
        theta.bt=1/(1.13^(1/12))^12;              %beta, discount factor (annual)
        theta.death=1/(40);                 %Probability agent dies (annual)
    else
        theta.bt=1/(1.13^(1/12));                %beta, discount factor (monthly)
        theta.death=1/(40*12);                 %Probability agent dies (monthly)
    end
    theta.bpw=1-xmin(11);                  %bpw is bargaining power of worker
    theta.bpf=1-theta.bpw;                    %Firm bargaining power (1-gamma)
    theta.n=1;                             %Measure of firms

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    theta.aup=0.05;                          %Probability of moving up
    theta.adown=0.05;                        %Probability of moving down
    theta.ulose=x(6);                       %Probability of moving donw in u
    theta.A=x(14);                          %Aggegate productivity parameter, still unsure hor to add this here
    theta.fcomp=x(3);                       %Complementarity in F
    theta.alpha_m=0.67;                      %Manager share
    theta.homeprod=0.7;                   %Home production
    theta.mnew=x(9);                        %Paramter of exponetial dist for newborn workers
    if tg.annual==1
        theta.lamu=1-(1-x(1))^12;                        %Prob of U find a firm
        theta.lam=1-(1-x(2))^12;                        %Prob of firm find another firm
        theta.del=1-(1-x(5))^12;                         %Separation rate
    else
        theta.lamu=1-(1-x(1));                        %Prob of U find a firm
        theta.lam=1-(1-x(2));                        %Prob of firm find another firm
        theta.del=1-(1-x(5));                         %Separation rate
    end
    if tg.split_top_bot==0
        theta.qup=0.05*ones(theta.ats,1);       %Probability of moving up in q. The qup for the last level of productivity is irrelevant
        theta.qdown=0*ones(theta.ats,1);        %Probability of moving down in q
    else
        theta.quptop=0.02;                     %Probability of moving up in q. The qup for the last level of productivity is irrelevant
        theta.qupbot=0.02;                     %Probability of moving up in q. The qup for the last level of productivity is irrelevant
        theta.qdowntop=0;                      %Probability of moving down in q
        theta.qdownbot=0;                      %Probability of moving down in q
    end
    theta.cost_d=2.25;                       %cost of demoting a manager to non manager  
    theta.cost_p=1.09;                       %cost of promoting a non manager to manager
    theta.phim=0.8;                         %Prob of being free manager (COUNTERFACTUAL)
    theta.phin=0.8;                         %Prob of being free non manager (COUNTERFACTUAL)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Grids
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cgrid.cost_d=linspace(0.5,1.4,15);
    cgrid.cost_p=linspace(0.02,0.2,15);
    if tg.split_top_bot==0
        cgrid.qup = repmat({linspace(0.08, 0.4, 5)},length(theta.qup),1);
    else
        cgrid.quptop=linspace(0.02,0.06,5);
        cgrid.qupbot=linspace(0.02,0.06,5);
    end
    cgrid.lamu=linspace(0.2,0.4,10);
    cgrid.lam=linspace(0.2,0.4,10);
    cgrid.alpha_m=linspace(0.5,0.9,10);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calibration weights  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cw.NiMi=2; 
    cw.MiNi=1; 
    cw.ManWorkerRatio=0.1; 
    cw.bcross=0.1; %Weights for moments

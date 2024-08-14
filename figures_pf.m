clear all
close all
location="local";
% location="hpc";
cd('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit')
% load('model_solutions.mat')
load('model_solutions_sp3local.mat')
load('xmin_results_combined.mat')




%Check sum of the eplus final distributions
nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");


x=xmin;

tpts=size(Vm,2);
ats=size(Vm,1);
cost_p=0.5                   ; %cost of promoting a non manager to manager
cost_d=0.5                   ; %cost of demoting a manager to non manager






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Predefine qualitative color palettes
load('color_style.mat'); % Ran in color_style.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%% Hiring Policies only
% Fix a ztilde and see if gets hires by all types of firm, coming from some states
%It will be a matriz for each type of origin of \ztilde

ztilde=1;
a=1;
atilde=1;
% Who hires ztilde from unemployed?
H_u=zeros(tpts+1,tpts+1);
H_u(1,1)=h_e_u(ztilde,1);
for z=2:tpts+1
    for q=2:tpts+1
        H_u(z,q)=h_t_u(ztilde,a,z-1,q-1);
    end
    H_u(z,1)=h_m_u(ztilde,a,z-1);
    H_u(1,z)=h_nm_u(ztilde,a,z-1);
end
plot_hire_mat(H_u, 'Hire z='+string(ztilde)+'  from u', 'Worker Type', 'Manager Type', 20, 'o', 'k', 'Figures_PFs/Hire_z1_u.png');
% H_u = [1 0 1;
%        0 1 0;
%        1 1 0];
  

% Who hires ztilde from firm with manager?
H_m=zeros(tpts+1,tpts+1);
H_m(1,1)=h_e_m(atilde,ztilde,a);
for z=2:tpts+1
    for q=2:tpts+1
        H_m(z,q)=h_t_m(atilde,ztilde,a,z-1,q-1);
    end
    H_m(z,1)=h_m_m(atilde,ztilde,a,z-1);
    H_m(1,z)=h_nm_m(atilde,ztilde,a,z-1);
end

plot_hire_mat(H_m, 'Hire z='+string(ztilde)+' from m', 'Worker Type', 'Manager Type', 20, 'o', 'k', 'Figures_PFs/Hire_z1_m.png');

% Who hires ztilde from firm with no manager?
H_nm=zeros(tpts+1,tpts+1);
H_nm(1,1)=h_e_nm(atilde,ztilde,a);
for z=2:tpts+1
    for q=2:tpts+1
        H_nm(z,q)=h_t_nm(atilde,ztilde,a,z-1,q-1);
    end
    H_nm(z,1)=h_m_nm(atilde,ztilde,a,z-1);
    H_nm(1,z)=h_nm_nm(atilde,ztilde,a,z-1);
end

plot_hire_mat(H_nm, 'Hire z='+string(ztilde)+' from nm', 'Worker Type', 'Manager Type', 20, 'o', 'k', 'Figures_PFs/Hire_z1_nm.png');

% Who hires ztilde who is a  manager from team (atilde,ztilde,qtilde)?

H_t_m=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    H_t_m(1,1,qtilde)=h_e_t_m(atilde,ztilde,qtilde,a);
    for z=2:tpts+1
        for q=2:tpts+1
            H_t_m(z,q,qtilde)=h_t_t_m(atilde,ztilde,qtilde,a,z-1,q-1);
        end
        H_t_m(z,1,qtilde)=h_m_t_m(atilde,ztilde,qtilde,a,z-1);
        H_t_m(1,z,qtilde)=h_nm_t_m(atilde,ztilde,qtilde,a,z-1);
    end
end
%Plot for lowerr and upper bounds of qtilde
plot_hire_mat(H_t_m(:,:,1), 'Hire z='+string(ztilde)+' from manager in team q=1', 'Worker Type', 'Manager Type', 20, 'o', 'k', 'Figures_PFs/Hire_z1_t_m_q1.png');
plot_hire_mat(H_t_m(:,:,2), 'Hire z='+string(ztilde)+' from manager in team q=2', 'Worker Type', 'Manager Type', 20, 'o', 'k', 'Figures_PFs/Hire_z1_t_m_q2.png');
plot_hire_mat(H_t_m(:,:,tpts), 'Hire z='+string(ztilde)+' from manager in team q='+string(tpts), 'Worker Type', 'Manager Type', 20, 'o', 'k', 'Figures_PFs/Hire_z1_t_m_qpng'); 
 



% Who hires ztilde who is a non manager from team (atilde,qtilde,ztilde)?
H_t_nm=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    H_t_nm(1,1,qtilde)=h_e_t_nm(atilde,qtilde,ztilde,a);
    for z=2:tpts+1
        for q=2:tpts+1
            H_t_nm(z,q,qtilde)=h_t_t_nm(atilde,qtilde,ztilde,a,z-1,q-1);
        end
        H_t_nm(z,1,qtilde)=h_m_t_nm(atilde,qtilde,ztilde,a,z-1);
        H_t_nm(1,z,qtilde)=h_nm_t_nm(atilde,qtilde,ztilde,a,z-1);
    end
end

%Plot for lowerr and upper bounds of qtilde
plot_hire_mat(H_t_nm(:,:,1), 'Hire z='+string(ztilde)+' from non manager in team q=1', 'Worker Type', 'Manager Type', 20, 'o', 'k', 'Figures_PFs/Hire_z1_t_nm_q1.png');
plot_hire_mat(H_t_nm(:,:,2), 'Hire z='+string(ztilde)+' from non manager in team q=2', 'Worker Type', 'Manager Type', 20, 'o', 'k', 'Figures_PFs/Hire_z1_t_nm_q2.png');
plot_hire_mat(H_t_nm(:,:,tpts), 'Hire z='+string(ztilde)+' from non manager in team q='+string(tpts), 'Worker Type', 'Manager Type', 20, 'o', 'k', 'Figures_PFs/Hire_z1_t_nm_q5.png');



%Look at all ztilde
for ztilde=1:tpts
    % Who hires ztilde from unemployed?
    H_u=zeros(tpts+1,tpts+1);
    H_u(1,1)=h_e_u(ztilde,1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_u(z,q)=h_t_u(ztilde,a,z-1,q-1);
        end
        H_u(z,1)=h_m_u(ztilde,a,z-1);
        H_u(1,z)=h_nm_u(ztilde,a,z-1);
    end

    % Who hires ztilde from firm with manager?
    H_m=zeros(tpts+1,tpts+1);
    H_m(1,1)=h_e_m(atilde,ztilde,a);
    for z=2:tpts+1
        for q=2:tpts+1
            H_m(z,q)=h_t_m(atilde,ztilde,a,z-1,q-1);
        end
        H_m(z,1)=h_m_m(atilde,ztilde,a,z-1);
        H_m(1,z)=h_nm_m(atilde,ztilde,a,z-1);
    end

    % Who hires ztilde from firm with no manager?
    H_nm=zeros(tpts+1,tpts+1);
    H_nm(1,1)=h_e_nm(atilde,ztilde,a);
    for z=2:tpts+1
        for q=2:tpts+1
            H_nm(z,q)=h_t_nm(atilde,ztilde,a,z-1,q-1);
        end
        H_nm(z,1)=h_m_nm(atilde,ztilde,a,z-1);
        H_nm(1,z)=h_nm_nm(atilde,ztilde,a,z-1);
    end
    % Who hires ztilde who is a  manager from team (atilde,ztilde,qtilde)?

    H_t_m=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        H_t_m(1,1,qtilde)=h_e_t_m(atilde,ztilde,qtilde,a);
        for z=2:tpts+1
            for q=2:tpts+1
                H_t_m(z,q,qtilde)=h_t_t_m(atilde,ztilde,qtilde,a,z-1,q-1);
            end
            H_t_m(z,1,qtilde)=h_m_t_m(atilde,ztilde,qtilde,a,z-1);
            H_t_m(1,z,qtilde)=h_nm_t_m(atilde,ztilde,qtilde,a,z-1);
        end
    end

    % Who hires ztilde who is a non manager from team (atilde,qtilde,ztilde)?
    H_t_nm=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        H_t_nm(1,1,qtilde)=h_e_t_nm(atilde,qtilde,ztilde,a);
        for z=2:tpts+1
            for q=2:tpts+1
                H_t_nm(z,q,qtilde)=h_t_t_nm(atilde,qtilde,ztilde,a,z-1,q-1);
            end
            H_t_nm(z,1,qtilde)=h_m_t_nm(atilde,qtilde,ztilde,a,z-1);
            H_t_nm(1,z,qtilde)=h_nm_t_nm(atilde,qtilde,ztilde,a,z-1);
        end
    end
    % Create a new figure with specific size
    figure('Units', 'normalized', 'Position', [0, 0, 1, 1]); % Fullscreen figure
    
    % Use tiledlayout for better subplot spacing
    tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % Generate each plot in the appropriate position
    plot_hire_matcombine(H_u, ['Hire z=' + string(ztilde) + ' from u'], 'Worker Type', 'Manager Type', 20, 'o', 'k', [], [3, 3, 1]);
    plot_hire_matcombine(H_m, ['Hire z=' + string(ztilde) + ' from m'], 'Worker Type', 'Manager Type', 20, 'o', 'k', [], [3, 3, 2]);
    plot_hire_matcombine(H_nm, ['Hire z=' + string(ztilde) + ' from nm'], 'Worker Type', 'Manager Type', 20, 'o', 'k', [], [3, 3, 3]);
    plot_hire_matcombine(H_t_m(:,:,1), ['Hire z=' + string(ztilde) + ' from manager in team q=1'], 'Worker Type', 'Manager Type', 20, 'o', 'k', [], [3, 3, 4]);
    plot_hire_matcombine(H_t_m(:,:,2), ['Hire z=' + string(ztilde) + ' from manager in team q=2'], 'Worker Type', 'Manager Type', 20, 'o', 'k', [], [3, 3, 5]);
    plot_hire_matcombine(H_t_m(:,:,tpts), ['Hire z=' + string(ztilde) + ' from manager in team q=' + string(tpts)], 'Worker Type', 'Manager Type', 20, 'o', 'k', [], [3, 3, 6]);
    plot_hire_matcombine(H_t_nm(:,:,1), ['Hire z=' + string(ztilde) + ' from non-manager in team q=1'], 'Worker Type', 'Manager Type', 20, 'o', 'k', [], [3, 3, 7]);
    plot_hire_matcombine(H_t_nm(:,:,2), ['Hire z=' + string(ztilde) + ' from non-manager in team q=2'], 'Worker Type', 'Manager Type', 20, 'o', 'k', [], [3, 3, 8]);
    plot_hire_matcombine(H_t_nm(:,:,tpts), ['Hire z=' + string(ztilde) + ' from non-manager in team q=' + string(tpts)], 'Worker Type', 'Manager Type', 20, 'o', 'k', [], [3, 3, 9]);
    
    % Save the figure as an image
    saveas(gcf, 'Figures_PFs/Hire_z' + string(ztilde) + '_all.png');
end


ztilde=1;
q_tilda=2;
a=1;
a_tilda=1;
q=3;
z=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the allocation policies (that are always conditional on hiring)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Who hires ztilde from unemployed?
% Putting as a manager (simbol black star)
H_u_pm=zeros(tpts+1,tpts+1);
H_u_pm(1,1)=h_e_u(ztilde,a)*p_e_m(a,ztilde);
for z=2:tpts+1
    H_u_pm(1,z)=h_nm_u(ztilde,a,z-1)*p_n_m(a,z-1,ztilde);
end

%Putting as a non manager (simbol black square)
H_u_pn=zeros(tpts+1,tpts+1);
H_u_pn(1,1)=h_e_u(ztilde,a)*p_e_n(a,ztilde);
for z=2:tpts+1
    H_u_pn(z,1)=h_m_u(ztilde,a,z-1)*p_m_n(a,z-1,ztilde);
end

%Putting as a manager firing the current manager (simbol magenta star)
H_u_pmu=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_u_pmu(z,q)=h_t_u(ztilde,a,z-1,q-1)*p_t_m_u(a,z-1,q-1,ztilde);
    end
    H_u_pmu(z,1)=h_m_u(ztilde,a,z-1)*p_m_m_u(a,z-1,ztilde);
end

%Putting as a manager demoting the current manager (simbol orange star)
H_u_pmd=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_u_pmd(z,q)=h_t_u(ztilde,a,z-1,q-1)*p_t_m_d(a,z-1,q-1,ztilde);
    end
    H_u_pmd(z,1)=h_m_u(ztilde,a,z-1)*p_m_m_d(a,z-1,ztilde);
end

%Putting as a non manager firing the current non manager (simbol magenta square)
H_u_pnu=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_u_pnu(z,q)=h_t_u(ztilde,a,z-1,q-1)*p_t_n_u(a,z-1,q-1,ztilde);
    end
    H_u_pnu(1,z)=h_nm_u(ztilde,a,z-1)*p_n_n_u(a,z-1,ztilde);
end

%Putting as a non manager promoting the current non manager (simbol blue square)
H_u_pnp=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_u_pnp(z,q)=h_t_u(ztilde,a,z-1,q-1)*p_t_n_p(a,z-1,q-1,ztilde);
    end
    H_u_pnp(1,z)=h_nm_u(ztilde,a,z-1)*p_n_n_p(a,z-1,ztilde);
end

%Plot the overlay of all these matrices
plot_hire_alloc({H_u_pm,H_u_pn,H_u_pmu,H_u_pmd,H_u_pnu,H_u_pnp}, 'Alloc of z='+string(ztilde)+' from u',...
    'Worker Type', 'Manager Type', 20, {'d', 's', 'd', 'd', 's', 's'}, {'k', 'k','m','#FFA500','m','b'},...
        {'Hire as M', 'Hire as NM', 'Hire as M w/ fire', 'Hire as M w/ demotion','Hire as NM w/ fire', 'Hire as NM w/ promotion'},...
    'Figures_PFs/Alloc_z1_u.png');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Who hires ztilde from firm with manager?
% Putting as a manager (simbol black star)
H_m_pm=zeros(tpts+1,tpts+1);
H_m_pm(1,1)=h_e_m(a_tilda,ztilde,a)*p_e_m(a,ztilde);
for z=2:tpts+1
    H_m_pm(1,z)=h_nm_m(a_tilda,ztilde,a,z-1)*p_n_m(a,z-1,ztilde);
end

%Putting as a non manager (simbol black square)
H_m_pn=zeros(tpts+1,tpts+1);
H_m_pn(1,1)=h_e_m(a_tilda,ztilde,a)*p_e_n(a,ztilde);
for z=2:tpts+1
    H_m_pn(z,1)=h_m_m(a_tilda,ztilde,a,z-1)*p_m_n(a,z-1,ztilde);
end

%Putting as a manager firing the current manager (simbol magenta star)
H_m_pmu=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_m_pmu(z,q)=h_t_m(a_tilda,ztilde,a,z-1,q-1)*p_t_m_u(a,z-1,q-1,ztilde);
    end
    H_m_pmu(z,1)=h_m_m(a_tilda,ztilde,a,z-1)*p_m_m_u(a,z-1,ztilde);
end

%Putting as a manager demoting the current manager (simbol orange star) 
H_m_pmd=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_m_pmd(z,q)=h_t_m(a_tilda,ztilde,a,z-1,q-1)*p_t_m_d(a,z-1,q-1,ztilde);
    end
    H_m_pmd(z,1)=h_m_m(a_tilda,ztilde,a,z-1)*p_m_m_d(a,z-1,ztilde);
end

%Putting as a non manager firing the current non manager (simbol magenta square)
H_m_pnu=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_m_pnu(z,q)=h_t_m(a_tilda,ztilde,a,z-1,q-1)*p_t_n_u(a,z-1,q-1,ztilde);
    end
    H_m_pnu(1,z)=h_nm_m(a_tilda,ztilde,a,z-1)*p_n_n_u(a,z-1,ztilde);
end

%Putting as a non manager promoting the current non manager (simbol blue square)
H_m_pnp=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_m_pnp(z,q)=h_t_m(a_tilda,ztilde,a,z-1,q-1)*p_t_n_p(a,z-1,q-1,ztilde);
    end
    H_m_pnp(1,z)=h_nm_m(a_tilda,ztilde,a,z-1)*p_n_n_p(a,z-1,ztilde);
end

%Plot the overlay of all these matrices
plot_hire_alloc({H_m_pm,H_m_pn,H_m_pmu,H_m_pmd,H_m_pnu,H_m_pnp}, 'Alloc of z='+string(ztilde)+' from m',...
    'Worker Type', 'Manager Type', 20, {'d', 's', 'd', 'd', 's', 's'}, {'k', 'k','m','#FFA500','m','b'},...
    {'Hire as M', 'Hire as NM', 'Hire as M w/ fire', 'Hire as M w/ demotion','Hire as NM w/ fire', 'Hire as NM w/ promotion'},...
    'Figures_PFs/Alloc_z1_m.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Who hires ztilde from firm with no manager?
% Putting as a manager (simbol black star)
H_nm_pm=zeros(tpts+1,tpts+1);
H_nm_pm(1,1)=h_e_nm(a_tilda,ztilde,a)*p_e_m(a,ztilde);
for z=2:tpts+1
    H_nm_pm(1,z)=h_nm_nm(a_tilda,ztilde,a,z-1)*p_n_m(a,z-1,ztilde);
end

%Putting as a non manager (simbol black square)
H_nm_pn=zeros(tpts+1,tpts+1);
H_nm_pn(1,1)=h_e_nm(a_tilda,ztilde,a)*p_e_n(a,ztilde);
for z=2:tpts+1
    H_nm_pn(z,1)=h_m_nm(a_tilda,ztilde,a,z-1)*p_m_n(a,z-1,ztilde);
end

%Putting as a manager firing the current manager (simbol magenta star)
H_nm_pmu=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_nm_pmu(z,q)=h_t_nm(a_tilda,ztilde,a,z-1,q-1)*p_t_m_u(a,z-1,q-1,ztilde);
    end
    H_nm_pmu(z,1)=h_m_nm(a_tilda,ztilde,a,z-1)*p_m_m_u(a,z-1,ztilde);
end

%Putting as a manager demoting the current manager (simbol orange star)
H_nm_pmd=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_nm_pmd(z,q)=h_t_nm(a_tilda,ztilde,a,z-1,q-1)*p_t_m_d(a,z-1,q-1,ztilde);
    end
    H_nm_pmd(z,1)=h_m_nm(a_tilda,ztilde,a,z-1)*p_m_m_d(a,z-1,ztilde);
end

%Putting as a non manager firing the current non manager (simbol magenta square)
H_nm_pnu=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_nm_pnu(z,q)=h_t_nm(a_tilda,ztilde,a,z-1,q-1)*p_t_n_u(a,z-1,q-1,ztilde);
    end
    H_nm_pnu(1,z)=h_nm_nm(a_tilda,ztilde,a,z-1)*p_n_n_u(a,z-1,ztilde);
end

%Putting as a non manager promoting the current non manager (simbol blue square)
H_nm_pnp=zeros(tpts+1,tpts+1);
for z=2:tpts+1
    for q=2:tpts+1
        H_nm_pnp(z,q)=h_t_nm(a_tilda,ztilde,a,z-1,q-1)*p_t_n_p(a,z-1,q-1,ztilde);
    end
    H_nm_pnp(1,z)=h_nm_nm(a_tilda,ztilde,a,z-1)*p_n_n_p(a,z-1,ztilde);
end

%Plot the overlay of all these matrices
plot_hire_alloc({H_nm_pm,H_nm_pn,H_nm_pmu,H_nm_pmd,H_nm_pnu,H_nm_pnp}, 'Alloc of z='+string(ztilde)+' from nm',...
    'Worker Type', 'Manager Type', 20, {'d', 's', 'd', 'd', 's', 's'}, {'k', 'k','m','#FFA500','m','b'},...
    {'Hire as M', 'Hire as NM', 'Hire as M w/ fire', 'Hire as M w/ demotion','Hire as NM w/ fire', 'Hire as NM w/ promotion'},...
    'Figures_PFs/Alloc_z1_nm.png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Who hires ztilde who is a manager from team (a_tilda,ztilde,q_tilda)?
% Putting as a manager (simbol black star)
H_t_m_pm=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    H_t_m_pm(1,1,qtilde)=h_e_t_m(a_tilda,ztilde,qtilde,a)*p_e_m(a,ztilde);
    for z=2:tpts+1
        H_t_m_pm(1,z,qtilde)=h_nm_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_n_m(a,z-1,ztilde);
    end
end

%Putting as a non manager (simbol black square)
H_t_m_pn=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    H_t_m_pn(1,1,qtilde)=h_e_t_m(a_tilda,ztilde,qtilde,a)*p_e_n(a,ztilde);
    for z=2:tpts+1
        H_t_m_pn(z,1,qtilde)=h_m_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_m_n(a,z-1,ztilde);
    end
end

%Putting as a manager firing the current manager (simbol magenta star)
H_t_m_pmu=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    for z=2:tpts+1
        for q=2:tpts+1
            H_t_m_pmu(z,q,qtilde)=h_t_t_m(a_tilda,ztilde,qtilde,a,z-1,q-1)*p_t_m_u(a,z-1,q-1,ztilde);
        end
        H_t_m_pmu(z,1,qtilde)=h_m_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_m_m_u(a,z-1,ztilde);
    end
end

%Putting as a manager demoting the current manager (simbol orange star)
H_t_m_pmd=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    for z=2:tpts+1
        for q=2:tpts+1
            H_t_m_pmd(z,q,qtilde)=h_t_t_m(a_tilda,ztilde,qtilde,a,z-1,q-1)*p_t_m_d(a,z-1,q-1,ztilde);
        end
        H_t_m_pmd(z,1,qtilde)=h_m_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_m_m_d(a,z-1,ztilde);
    end
end

%Putting as a non manager firing the current non manager (simbol magenta square)
H_t_m_pnu=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    for z=2:tpts+1
        for q=2:tpts+1
            H_t_m_pnu(z,q,qtilde)=h_t_t_m(a_tilda,ztilde,qtilde,a,z-1,q-1)*p_t_n_u(a,z-1,q-1,ztilde);
        end
        H_t_m_pnu(1,z,qtilde)=h_nm_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_n_n_u(a,z-1,ztilde);
    end
end

%Putting as a non manager promoting the current non manager (simbol blue square)
H_t_m_pnp=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    for z=2:tpts+1
        for q=2:tpts+1
            H_t_m_pnp(z,q,qtilde)=h_t_t_m(a_tilda,ztilde,qtilde,a,z-1,q-1)*p_t_n_p(a,z-1,q-1,ztilde);
        end
        H_t_m_pnp(1,z,qtilde)=h_nm_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_n_n_p(a,z-1,ztilde);
    end
end

%Plot the overlay of all these matrices
plot_hire_alloc({H_t_m_pm(:,:,1),H_t_m_pn(:,:,1),H_t_m_pmu(:,:,1),H_t_m_pmd(:,:,1),H_t_m_pnu(:,:,1),H_t_m_pnp(:,:,1)}, 'Alloc of z='+string(ztilde)+' from manager in team q=1',...
    'Worker Type', 'Manager Type', 20, {'d', 's', 'd', 'd', 's', 's'}, {'k', 'k','m','#FFA500','m','b'},...
    {'Hire as M', 'Hire as NM', 'Hire as M w/ fire', 'Hire as M w/ demotion','Hire as NM w/ fire', 'Hire as NM w/ promotion'},...
    'Figures_PFs/Alloc_z1_t_m_q1.png');

plot_hire_alloc({H_t_m_pm(:,:,2),H_t_m_pn(:,:,2),H_t_m_pmu(:,:,2),H_t_m_pmd(:,:,2),H_t_m_pnu(:,:,2),H_t_m_pnp(:,:,2)}, 'Alloc of z='+string(ztilde)+' from manager in team q=2',...
    'Worker Type', 'Manager Type', 20, {'d', 's', 'd', 'd', 's', 's'}, {'k', 'k','m','#FFA500','m','b'},...
    {'Hire as M', 'Hire as NM', 'Hire as M w/ fire', 'Hire as M w/ demotion','Hire as NM w/ fire', 'Hire as NM w/ promotion'},...
    'Figures_PFs/Alloc_z1_t_m_q2.png');

plot_hire_alloc({H_t_m_pm(:,:,tpts),H_t_m_pn(:,:,tpts),H_t_m_pmu(:,:,tpts),H_t_m_pmd(:,:,tpts),H_t_m_pnu(:,:,tpts),H_t_m_pnp(:,:,tpts)}, 'Alloc of z='+string(ztilde)+' from manager in team q='+string(tpts),...
    'Worker Type', 'Manager Type', 20, {'d', 's', 'd', 'd', 's', 's'}, {'k', 'k','m','#FFA500','m','b'},...
    {'Hire as M', 'Hire as NM', 'Hire as M w/ fire', 'Hire as M w/ demotion','Hire as NM w/ fire', 'Hire as NM w/ promotion'},...
    'Figures_PFs/Alloc_z1_t_m_q3.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Who hires ztilde who is a non manager from team (a_tilda,q_tilda,ztilde)?
% Putting as a manager (simbol black star)
H_t_nm_pm=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    H_t_nm_pm(1,1,qtilde)=h_e_t_nm(a_tilda,qtilde,ztilde,a)*p_e_m(a,ztilde);
    for z=2:tpts+1
        H_t_nm_pm(1,z,qtilde)=h_nm_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_n_m(a,z-1,ztilde);
    end
end

%Putting as a non manager (simbol black square)
H_t_nm_pn=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    H_t_nm_pn(1,1,qtilde)=h_e_t_nm(a_tilda,qtilde,ztilde,a)*p_e_n(a,ztilde);
    for z=2:tpts+1
        H_t_nm_pn(z,1,qtilde)=h_m_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_m_n(a,z-1,ztilde);
    end
end

%Putting as a manager firing the current manager (simbol magenta star)
H_t_nm_pmu=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    for z=2:tpts+1
        for q=2:tpts+1
            H_t_nm_pmu(z,q,qtilde)=h_t_t_nm(a_tilda,qtilde,ztilde,a,z-1,q-1)*p_t_m_u(a,z-1,q-1,ztilde);
        end
        H_t_nm_pmu(z,1,qtilde)=h_m_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_m_m_u(a,z-1,ztilde);
    end
end

%Putting as a manager demoting the current manager (simbol orange star)
H_t_nm_pmd=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    for z=2:tpts+1
        for q=2:tpts+1
            H_t_nm_pmd(z,q,qtilde)=h_t_t_nm(a_tilda,qtilde,ztilde,a,z-1,q-1)*p_t_m_d(a,z-1,q-1,ztilde);
        end
        H_t_nm_pmd(z,1,qtilde)=h_m_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_m_m_d(a,z-1,ztilde);
    end
end

%Putting as a non manager firing the current non manager (simbol magenta square)
H_t_nm_pnu=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    for z=2:tpts+1
        for q=2:tpts+1
            H_t_nm_pnu(z,q,qtilde)=h_t_t_nm(a_tilda,qtilde,ztilde,a,z-1,q-1)*p_t_n_u(a,z-1,q-1,ztilde);
        end
        H_t_nm_pnu(1,z,qtilde)=h_nm_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_n_n_u(a,z-1,ztilde);
    end
end

%Putting as a non manager promoting the current non manager (simbol blue square)
H_t_nm_pnp=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
for qtilde=1:tpts
    for z=2:tpts+1
        for q=2:tpts+1
            H_t_nm_pnp(z,q,qtilde)=h_t_t_nm(a_tilda,qtilde,ztilde,a,z-1,q-1)*p_t_n_p(a,z-1,q-1,ztilde);
        end
        H_t_nm_pnp(1,z,qtilde)=h_nm_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_n_n_p(a,z-1,ztilde);
    end
end

%Plot the overlay of all these matrices
plot_hire_alloc({H_t_nm_pm(:,:,1),H_t_nm_pn(:,:,1),H_t_nm_pmu(:,:,1),H_t_nm_pmd(:,:,1),H_t_nm_pnu(:,:,1),H_t_nm_pnp(:,:,1)}, 'Alloc of z='+string(ztilde)+' from non manager in team q=1',...
    'Worker Type', 'Manager Type', 20, {'d', 's', 'd', 'd', 's', 's'}, {'k', 'k','m','#FFA500','m','b'},...
    {'Hire as M', 'Hire as NM', 'Hire as M w/ fire', 'Hire as M w/ demotion','Hire as NM w/ fire', 'Hire as NM w/ promotion'},...
    'Figures_PFs/Alloc_z1_t_nm_q1.png');

plot_hire_alloc({H_t_nm_pm(:,:,2),H_t_nm_pn(:,:,2),H_t_nm_pmu(:,:,2),H_t_nm_pmd(:,:,2),H_t_nm_pnu(:,:,2),H_t_nm_pnp(:,:,2)}, 'Alloc of z='+string(ztilde)+' from non manager in team q=2',...
    'Worker Type', 'Manager Type', 20, {'d', 's', 'd', 'd', 's', 's'}, {'k', 'k','m','#FFA500','m','b'},...
    {'Hire as M', 'Hire as NM', 'Hire as M w/ fire', 'Hire as M w/ demotion','Hire as NM w/ fire', 'Hire as NM w/ promotion'},...
    'Figures_PFs/Alloc_z1_t_nm_q2.png');

plot_hire_alloc({H_t_nm_pm(:,:,tpts),H_t_nm_pn(:,:,tpts),H_t_nm_pmu(:,:,tpts),H_t_nm_pmd(:,:,tpts),H_t_nm_pnu(:,:,tpts),H_t_nm_pnp(:,:,tpts)}, 'Alloc of z='+string(ztilde)+' from non manager in team q='+string(tpts),...
    'Worker Type', 'Manager Type', 20, {'d', 's', 'd', 'd', 's', 's'}, {'k', 'k','m','#FFA500','m','b'},...
    {'Hire as M', 'Hire as NM', 'Hire as M w/ fire', 'Hire as M w/ demotion','Hire as NM w/ fire', 'Hire as NM w/ promotion'},...
    'Figures_PFs/Alloc_z1_t_nm_q3.png');
    

a=1;
a_tilda=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% For every ztilde, plot the allocation policies
for ztilde=1:tpts
    %Who hires ztilde from unemployed?
    % Putting as a manager (simbol black star)
    H_u_pm=zeros(tpts+1,tpts+1);
    H_u_pm(1,1)=h_e_u(ztilde,a)*p_e_m(a,ztilde);
    for z=2:tpts+1
        H_u_pm(1,z)=h_nm_u(ztilde,a,z-1)*p_n_m(a,z-1,ztilde);
    end
    
    %Putting as a non manager (simbol black square)
    H_u_pn=zeros(tpts+1,tpts+1);
    H_u_pn(1,1)=h_e_u(ztilde,a)*p_e_n(a,ztilde);
    for z=2:tpts+1
        H_u_pn(z,1)=h_m_u(ztilde,a,z-1)*p_m_n(a,z-1,ztilde);
    end
    
    %Putting as a manager firing the current manager (simbol magenta star)
    H_u_pmu=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_u_pmu(z,q)=h_t_u(ztilde,a,z-1,q-1)*p_t_m_u(a,z-1,q-1,ztilde);
        end
        H_u_pmu(z,1)=h_m_u(ztilde,a,z-1)*p_m_m_u(a,z-1,ztilde);
    end
    
    %Putting as a manager demoting the current manager (simbol orange star)
    H_u_pmd=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_u_pmd(z,q)=h_t_u(ztilde,a,z-1,q-1)*p_t_m_d(a,z-1,q-1,ztilde);
        end
        H_u_pmd(z,1)=h_m_u(ztilde,a,z-1)*p_m_m_d(a,z-1,ztilde);
    end
    
    %Putting as a non manager firing the current non manager (simbol magenta square)
    H_u_pnu=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_u_pnu(z,q)=h_t_u(ztilde,a,z-1,q-1)*p_t_n_u(a,z-1,q-1,ztilde);
        end
        H_u_pnu(1,z)=h_nm_u(ztilde,a,z-1)*p_n_n_u(a,z-1,ztilde);
    end
    
    %Putting as a non manager promoting the current non manager (simbol blue square)
    H_u_pnp=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_u_pnp(z,q)=h_t_u(ztilde,a,z-1,q-1)*p_t_n_p(a,z-1,q-1,ztilde);
        end
        H_u_pnp(1,z)=h_nm_u(ztilde,a,z-1)*p_n_n_p(a,z-1,ztilde);
    end
    %Who hires ztilde from firm with manager?
    % Putting as a manager (simbol black star)
    H_m_pm=zeros(tpts+1,tpts+1);
    H_m_pm(1,1)=h_e_m(a_tilda,ztilde,a)*p_e_m(a,ztilde);
    for z=2:tpts+1
        H_m_pm(1,z)=h_nm_m(a_tilda,ztilde,a,z-1)*p_n_m(a,z-1,ztilde);
    end
    
    %Putting as a non manager (simbol black square)
    H_m_pn=zeros(tpts+1,tpts+1);
    H_m_pn(1,1)=h_e_m(a_tilda,ztilde,a)*p_e_n(a,ztilde);
    for z=2:tpts+1
        H_m_pn(z,1)=h_m_m(a_tilda,ztilde,a,z-1)*p_m_n(a,z-1,ztilde);
    end
    
    %Putting as a manager firing the current manager (simbol magenta star)
    H_m_pmu=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_m_pmu(z,q)=h_t_m(a_tilda,ztilde,a,z-1,q-1)*p_t_m_u(a,z-1,q-1,ztilde);
        end
        H_m_pmu(z,1)=h_m_m(a_tilda,ztilde,a,z-1)*p_m_m_u(a,z-1,ztilde);
    end
    
    %Putting as a manager demoting the current manager (simbol orange star) 
    H_m_pmd=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_m_pmd(z,q)=h_t_m(a_tilda,ztilde,a,z-1,q-1)*p_t_m_d(a,z-1,q-1,ztilde);
        end
        H_m_pmd(z,1)=h_m_m(a_tilda,ztilde,a,z-1)*p_m_m_d(a,z-1,ztilde);
    end
    
    %Putting as a non manager firing the current non manager (simbol magenta square)
    H_m_pnu=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_m_pnu(z,q)=h_t_m(a_tilda,ztilde,a,z-1,q-1)*p_t_n_u(a,z-1,q-1,ztilde);
        end
        H_m_pnu(1,z)=h_nm_m(a_tilda,ztilde,a,z-1)*p_n_n_u(a,z-1,ztilde);
    end
    
    %Putting as a non manager promoting the current non manager (simbol blue square)
    H_m_pnp=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_m_pnp(z,q)=h_t_m(a_tilda,ztilde,a,z-1,q-1)*p_t_n_p(a,z-1,q-1,ztilde);
        end
        H_m_pnp(1,z)=h_nm_m(a_tilda,ztilde,a,z-1)*p_n_n_p(a,z-1,ztilde);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Who hires ztilde from firm with no manager?
    % Putting as a manager (simbol black star)
    H_nm_pm=zeros(tpts+1,tpts+1);
    H_nm_pm(1,1)=h_e_nm(a_tilda,ztilde,a)*p_e_m(a,ztilde);
    for z=2:tpts+1
        H_nm_pm(1,z)=h_nm_nm(a_tilda,ztilde,a,z-1)*p_n_m(a,z-1,ztilde);
    end
    
    %Putting as a non manager (simbol black square)
    H_nm_pn=zeros(tpts+1,tpts+1);
    H_nm_pn(1,1)=h_e_nm(a_tilda,ztilde,a)*p_e_n(a,ztilde);
    for z=2:tpts+1
        H_nm_pn(z,1)=h_m_nm(a_tilda,ztilde,a,z-1)*p_m_n(a,z-1,ztilde);
    end
    
    %Putting as a manager firing the current manager (simbol magenta star)
    H_nm_pmu=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_nm_pmu(z,q)=h_t_nm(a_tilda,ztilde,a,z-1,q-1)*p_t_m_u(a,z-1,q-1,ztilde);
        end
        H_nm_pmu(z,1)=h_m_nm(a_tilda,ztilde,a,z-1)*p_m_m_u(a,z-1,ztilde);
    end
    
    %Putting as a manager demoting the current manager (simbol orange star)
    H_nm_pmd=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_nm_pmd(z,q)=h_t_nm(a_tilda,ztilde,a,z-1,q-1)*p_t_m_d(a,z-1,q-1,ztilde);
        end
        H_nm_pmd(z,1)=h_m_nm(a_tilda,ztilde,a,z-1)*p_m_m_d(a,z-1,ztilde);
    end
    
    %Putting as a non manager firing the current non manager (simbol magenta square)
    H_nm_pnu=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_nm_pnu(z,q)=h_t_nm(a_tilda,ztilde,a,z-1,q-1)*p_t_n_u(a,z-1,q-1,ztilde);
        end
        H_nm_pnu(1,z)=h_nm_nm(a_tilda,ztilde,a,z-1)*p_n_n_u(a,z-1,ztilde);
    end
    
    %Putting as a non manager promoting the current non manager (simbol blue square)
    H_nm_pnp=zeros(tpts+1,tpts+1);
    for z=2:tpts+1
        for q=2:tpts+1
            H_nm_pnp(z,q)=h_t_nm(a_tilda,ztilde,a,z-1,q-1)*p_t_n_p(a,z-1,q-1,ztilde);
        end
        H_nm_pnp(1,z)=h_nm_nm(a_tilda,ztilde,a,z-1)*p_n_n_p(a,z-1,ztilde);
    end

    %Who hires ztilde who is a manager from team (a_tilda,ztilde,q_tilda)?
    % Putting as a manager (simbol black star)
    H_t_m_pm=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        H_t_m_pm(1,1,qtilde)=h_e_t_m(a_tilda,ztilde,qtilde,a)*p_e_m(a,ztilde);
        for z=2:tpts+1
            H_t_m_pm(1,z,qtilde)=h_nm_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_n_m(a,z-1,ztilde);
        end
    end
    
    %Putting as a non manager (simbol black square)
    H_t_m_pn=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        H_t_m_pn(1,1,qtilde)=h_e_t_m(a_tilda,ztilde,qtilde,a)*p_e_n(a,ztilde);
        for z=2:tpts+1
            H_t_m_pn(z,1,qtilde)=h_m_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_m_n(a,z-1,ztilde);
        end
    end
    
    %Putting as a manager firing the current manager (simbol magenta star)
    H_t_m_pmu=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        for z=2:tpts+1
            for q=2:tpts+1
                H_t_m_pmu(z,q,qtilde)=h_t_t_m(a_tilda,ztilde,qtilde,a,z-1,q-1)*p_t_m_u(a,z-1,q-1,ztilde);
            end
            H_t_m_pmu(z,1,qtilde)=h_m_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_m_m_u(a,z-1,ztilde);
        end
    end
    
    %Putting as a manager demoting the current manager (simbol orange star)
    H_t_m_pmd=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        for z=2:tpts+1
            for q=2:tpts+1
                H_t_m_pmd(z,q,qtilde)=h_t_t_m(a_tilda,ztilde,qtilde,a,z-1,q-1)*p_t_m_d(a,z-1,q-1,ztilde);
            end
            H_t_m_pmd(z,1,qtilde)=h_m_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_m_m_d(a,z-1,ztilde);
        end
    end
    
    %Putting as a non manager firing the current non manager (simbol magenta square)
    H_t_m_pnu=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        for z=2:tpts+1
            for q=2:tpts+1
                H_t_m_pnu(z,q,qtilde)=h_t_t_m(a_tilda,ztilde,qtilde,a,z-1,q-1)*p_t_n_u(a,z-1,q-1,ztilde);
            end
            H_t_m_pnu(1,z,qtilde)=h_nm_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_n_n_u(a,z-1,ztilde);
        end
    end
    
    %Putting as a non manager promoting the current non manager (simbol blue square)
    H_t_m_pnp=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        for z=2:tpts+1
            for q=2:tpts+1
                H_t_m_pnp(z,q,qtilde)=h_t_t_m(a_tilda,ztilde,qtilde,a,z-1,q-1)*p_t_n_p(a,z-1,q-1,ztilde);
            end
            H_t_m_pnp(1,z,qtilde)=h_nm_t_m(a_tilda,ztilde,qtilde,a,z-1)*p_n_n_p(a,z-1,ztilde);
        end
    end%Who hires ztilde who is a non manager from team (a_tilda,q_tilda,ztilde)?
    % Putting as a manager (simbol black star)
    H_t_nm_pm=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        H_t_nm_pm(1,1,qtilde)=h_e_t_nm(a_tilda,qtilde,ztilde,a)*p_e_m(a,ztilde);
        for z=2:tpts+1
            H_t_nm_pm(1,z,qtilde)=h_nm_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_n_m(a,z-1,ztilde);
        end
    end
    
    %Putting as a non manager (simbol black square)
    H_t_nm_pn=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        H_t_nm_pn(1,1,qtilde)=h_e_t_nm(a_tilda,qtilde,ztilde,a)*p_e_n(a,ztilde);
        for z=2:tpts+1
            H_t_nm_pn(z,1,qtilde)=h_m_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_m_n(a,z-1,ztilde);
        end
    end
    
    %Putting as a manager firing the current manager (simbol magenta star)
    H_t_nm_pmu=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        for z=2:tpts+1
            for q=2:tpts+1
                H_t_nm_pmu(z,q,qtilde)=h_t_t_nm(a_tilda,qtilde,ztilde,a,z-1,q-1)*p_t_m_u(a,z-1,q-1,ztilde);
            end
            H_t_nm_pmu(z,1,qtilde)=h_m_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_m_m_u(a,z-1,ztilde);
        end
    end
    
    %Putting as a manager demoting the current manager (simbol orange star)
    H_t_nm_pmd=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        for z=2:tpts+1
            for q=2:tpts+1
                H_t_nm_pmd(z,q,qtilde)=h_t_t_nm(a_tilda,qtilde,ztilde,a,z-1,q-1)*p_t_m_d(a,z-1,q-1,ztilde);
            end
            H_t_nm_pmd(z,1,qtilde)=h_m_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_m_m_d(a,z-1,ztilde);
        end
    end
    
    %Putting as a non manager firing the current non manager (simbol magenta square)
    H_t_nm_pnu=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        for z=2:tpts+1
            for q=2:tpts+1
                H_t_nm_pnu(z,q,qtilde)=h_t_t_nm(a_tilda,qtilde,ztilde,a,z-1,q-1)*p_t_n_u(a,z-1,q-1,ztilde);
            end
            H_t_nm_pnu(1,z,qtilde)=h_nm_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_n_n_u(a,z-1,ztilde);
        end
    end
    
    %Putting as a non manager promoting the current non manager (simbol blue square)
    H_t_nm_pnp=zeros(tpts+1,tpts+1,tpts); %do one for all possible qtilde
    for qtilde=1:tpts
        for z=2:tpts+1
            for q=2:tpts+1
                H_t_nm_pnp(z,q,qtilde)=h_t_t_nm(a_tilda,qtilde,ztilde,a,z-1,q-1)*p_t_n_p(a,z-1,q-1,ztilde);
            end
            H_t_nm_pnp(1,z,qtilde)=h_nm_t_nm(a_tilda,qtilde,ztilde,a,z-1)*p_n_n_p(a,z-1,ztilde);
        end
    end

    %All the matrices and titles for the plots
    H_matrices={{H_u_pm,H_u_pn,H_u_pmu,H_u_pmd,H_u_pnu,H_u_pnp},...
        {H_m_pm,H_m_pn,H_m_pmu,H_m_pmd,H_m_pnu,H_m_pnp},...
        {H_nm_pm,H_nm_pn,H_nm_pmu,H_nm_pmd,H_nm_pnu,H_nm_pnp},...
        {H_t_m_pm(:,:,1),H_t_m_pn(:,:,1),H_t_m_pmu(:,:,1),H_t_m_pmd(:,:,1),H_t_m_pnu(:,:,1),H_t_m_pnp(:,:,1)},...
        {H_t_m_pm(:,:,4),H_t_m_pn(:,:,4),H_t_m_pmu(:,:,4),H_t_m_pmd(:,:,4),H_t_m_pnu(:,:,4),H_t_m_pnp(:,:,4)},...
        {H_t_m_pm(:,:,tpts),H_t_m_pn(:,:,tpts),H_t_m_pmu(:,:,tpts),H_t_m_pmd(:,:,tpts),H_t_m_pnu(:,:,tpts),H_t_m_pnp(:,:,tpts)},...
        {H_t_nm_pm(:,:,1),H_t_nm_pn(:,:,1),H_t_nm_pmu(:,:,1),H_t_nm_pmd(:,:,1),H_t_nm_pnu(:,:,1),H_t_nm_pnp(:,:,1)},...
        {H_t_nm_pm(:,:,4),H_t_nm_pn(:,:,4),H_t_nm_pmu(:,:,4),H_t_nm_pmd(:,:,4),H_t_nm_pnu(:,:,4),H_t_nm_pnp(:,:,4)},...
        {H_t_nm_pm(:,:,tpts),H_t_nm_pn(:,:,tpts),H_t_nm_pmu(:,:,tpts),H_t_nm_pmd(:,:,tpts),H_t_nm_pnu(:,:,tpts),H_t_nm_pnp(:,:,tpts)}
        };    
    
    titles={'from Unemp.',...
        'From firm with Manager',...
        'From frim with Non-Manager',...
        'From team (z='+string(ztilde)+',q=1)',...
        'From team (z='+string(ztilde)+',q=4)',...
        'From team (z='+string(ztilde)+',q='+string(tpts)+')',...
        'From team (q=1,z='+string(ztilde)+')',...
        'From team (q=4,z='+string(ztilde)+')',...
        'From team (q='+string(tpts)+',z='+string(ztilde)+')'};
    
    plot_alloc_matcombine(H_matrices, titles,'Allocation of hire z='+string(ztilde) ,'Worker Type', 'Manager Type', 16, {'d', 's', 'd', 'd', 's', 's'}, {blue5, orange5,greenColor,'#FFA500',pinkColor,'b'},...
        {'Hire as M', 'Hire as NM', 'Hire as M w/ fire', 'Hire as M w/ demotion','Hire as NM w/ fire', 'Hire as NM w/ promotion'},...
        'Figures_PFs/Alloc_z'+string(ztilde)+'_all.png');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reallocation policies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d_m, r_m, d_n, r_n, d_t_m, d_t_n, r_t_m, r_t_n, d_t_b]...
    =reallocation_policies(Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p);

for a=1:ats
    %Conditional on a
    %Who fires a manager? symbol blue right arrow
    R_d_m=zeros(tpts+1,tpts+1,ats);
    for a=1:ats
        for z=2:tpts+1
            for q=2:tpts+1
                R_d_m(z,q,a)=d_t_m(a,z-1,q-1)*(1-r_t_n(a,z-1,q-1));
            end
            R_d_m(z,1,a)=d_m(a,z-1);
        end
    end
    
    %Who fires a non manager? symbol orange right arrow
    R_d_n=zeros(tpts+1,tpts+1,ats);
    for a=1:ats
        for z=2:tpts+1
            for q=2:tpts+1
                R_d_n(z,q,a)=d_t_n(a,z-1,q-1)*(1-r_t_m(a,z-1,q-1));
            end
            R_d_n(z,1,a)=d_n(a,z-1);
        end
    end
    
    % Who demotes a manager? symbol blue down arrow
    R_r_m=zeros(tpts+1,tpts+1,ats);
    for a=1:ats
        for z=2:tpts+1
            R_r_m(z,1,a)=r_m(a,z-1);
        end
    end
    
    % Who Promotes a non manager? symbol orange up arrow
    R_r_n=zeros(tpts+1,tpts+1,ats);
    for a=1:ats
        for z=2:tpts+1
            R_r_n(1,z,a)=r_n(a,z-1);
        end
    end
    
    %Who fires M with promotion symbol green up arrow
    R_dp_m=zeros(tpts+1,tpts+1,ats);
    for a=1:ats
        for z=2:tpts
            for q=2:tpts
                R_dp_m(z,q,a)=d_t_m(a,z-1,q-1)*r_t_n(a,z-1,q-1);
            end
        end
    end
    
    %Who fires NM with demotion 
    R_dp_n=zeros(tpts+1,tpts+1,ats);
    for a=1:ats
        for z=2:tpts
            for q=2:tpts
                R_dp_n(z,q,a)=d_t_n(a,z-1,q-1)*r_t_m(a,z-1,q-1);
            end
        end
    end
    
    %Who reallocastes both M and NM purple circle
    R_rr_b=zeros(tpts+1,tpts+1,ats);
    for a=1:ats
        for z=2:tpts
            for q=2:tpts
                R_rr_b(z,q,a)=r_t_m(a,z-1,q-1)*r_t_n(a,z-1,q-1);
            end
        end
    end
    
    %Who fires both M and NM red circle
    R_dd_b=zeros(tpts+1,tpts+1,ats);
    for a=1:ats
        for z=2:tpts
            for q=2:tpts
                R_dd_b(z,q,a)=d_t_b(a,z-1,q-1);
            end
        end
    end
    %Plot the reallocation policies
    plot_hire_alloc({R_d_m(:,:,a),R_d_n(:,:,a),R_r_m(:,:,a),R_r_n(:,:,a),R_dp_m(:,:,a),R_dp_n(:,:,a),R_rr_b(:,:,a),R_dd_b(:,:,a)}, 'Reallocation policies conditional on a='+string(a),...
        'Worker Type', 'Manager Type', 16, {'d', 's', 'd', 'd', 's', 's', 'o', 'o'}, {'b', '#FFA500','b','#FFA500','g','g','m','r'},...
        {'Fire M', 'Fire NM', 'Demote M', 'Promote NM','Fire M w/ promotion', 'Fire NM w/ demotion','Reallocate both','Fire both'},...
        'Figures_PFs/Reallocation_a'+string(a)+'.png');
end

%Plot them combined
R_matrices={{R_d_m(:,:,1),R_d_n(:,:,1),R_r_m(:,:,1),R_r_n(:,:,1),R_dp_m(:,:,1),R_dp_n(:,:,1),R_rr_b(:,:,1),R_dd_b(:,:,1)},...
    {R_d_m(:,:,3),R_d_n(:,:,3),R_r_m(:,:,3),R_r_n(:,:,3),R_dp_m(:,:,3),R_dp_n(:,:,3),R_rr_b(:,:,3),R_dd_b(:,:,3)},...
    {R_d_m(:,:,ats),R_d_n(:,:,ats),R_r_m(:,:,ats),R_r_n(:,:,ats),R_dp_m(:,:,ats),R_dp_n(:,:,ats),R_rr_b(:,:,ats),R_dd_b(:,:,ats)}};

titles={'a=1',...
    'a=3',...
    'a='+string(ats)};
    
    plot_alloc_matcombine(R_matrices, titles,'Reallocation policies','Worker Type', 'Manager Type', 16, {'d', 's', 'd', 'd', 's', 's', 'o', 'o'}, {blue4,redDark ,blue4,orange4,'g','g',purple7,'r'},...
    {'Fire M', 'Fire NM', 'Demote M', 'Promote NM','Fire M w/ promotion', 'Fire NM w/ demotion','Reallocate both','Fire both'},...
    'Figures_PFs/Reallocation_all.png');
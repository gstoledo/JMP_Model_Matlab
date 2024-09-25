% Compute all gains from trade from firms with non-managers
function [gt_meet_u, gt_meet_m, gt_meet_nm, gt_meet_t_m, gt_meet_t_nm, gt_meet_t] = nf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,nm_penal,cost_p)
    tol=1e-8;
    %Gains from trade from meeting unemployed
    gt_meet_u = zeros(tpts,ats,tpts);
    % Compute gt_meet_u
    %gt_meet_u(z_tilda,a,q)
    for a = 1:ats
        for q = 1:tpts
            for z_tilda = 1:tpts
                nu= [Vnh(a,z_tilda)+Uh(q), Vth(a,q,z_tilda)-cost_p ,Vth(a,z_tilda,q)];
                gt_meet_u(z_tilda,a,q) = max(max(nu) - Vnh(a,q)   -Uh(z_tilda),tol);
            end
        end
    end
    
    %Gains from trade from meeting firm with manager
    gt_meet_m = zeros(ats,tpts,ats,tpts);
    % Compute gt_meet_m
    %gt_meet_m(a_tilda,z_tilda,a,q)
    for a = 1:ats
        for q = 1:tpts
            for z_tilda = 1:tpts
                for a_tilda = 1:ats
                    nu=[Vnh(a,z_tilda)+Uh(q), Vth(a,q,z_tilda)-cost_p ,Vth(a,z_tilda,q)];
                    gt_meet_m(a_tilda,z_tilda,a,q) = max(max(nu) - Vnh(a,q) - Vmh(a_tilda,z_tilda) + Veh(a_tilda), tol);
                end
            end
        end
    end
    
    %Gains from trade from meeting firm with non-manager
    gt_meet_nm = zeros(ats,tpts,ats,tpts);
    % Compute gt_meet_nm
    %gt_meet_nm(a_tilda,z_tilda,a,q)    
    %Vmhp=Non_man_penalty_loop(Vmh,tpts,ats,nm_penal);
    Vnhp=Non_man_penalty_loop(Vnh,tpts,ats,nm_penal);
    Vthpm=Non_man_penalty_full_m_loop(Vth,tpts,ats,nm_penal);
    Vthpn=Non_man_penalty_full_n_loop(Vth,tpts,ats,nm_penal);
    for a = 1:ats
        for q = 1:tpts
            for z_tilda = 1:tpts
                for a_tilda = 1:ats
                    nu=[Vnhp(a,z_tilda)+Uh(q),Vthpm(a,z_tilda,q),Vthpn(a,q,z_tilda)-cost_p];
                    gt_meet_nm(a_tilda,z_tilda,a,q) = max(max(nu) - Vnh(a,q) - Vnh(a_tilda,z_tilda) + Veh(a_tilda), tol);
                end
            end
        end
    end
    
    %Gains from trade from meeting firm with team
    gt_meet_t = zeros(ats,tpts,tpts,ats,tpts);
    % Grabbing the manager
    %gt_meet_t_m(a_tilda,z_tilda,q_tilda,a,q)   
    gt_meet_t_m= zeros(ats,tpts,tpts,ats,tpts);
    for a = 1:ats
        for q = 1:tpts
            for q_tilda = 1:tpts
                for z_tilda = 1:tpts
                    for a_tilda = 1:ats
                        nu=[Vnh(a,z_tilda)+Uh(q),Vth(a,q,z_tilda)-cost_p,Vth(a,z_tilda,q)];
                        gt_meet_t_m(a_tilda,z_tilda,q_tilda,a,q) = max(max(nu) - Vnh(a,q) - Vth(a_tilda,z_tilda,q_tilda) + Vnh(a_tilda,q_tilda), tol);
                    end
                end
            end
        end
    end
    
    %Getting the non-manager
    %gt_meet_t_nm(a_tilda,z_tilda,q_tilda,a,q)
    gt_meet_t_nm= zeros(ats,tpts,tpts,ats,tpts);
    for a = 1:ats
        for q = 1:tpts
            for q_tilda = 1:tpts
                for z_tilda = 1:tpts
                    for a_tilda = 1:ats
                        nu=[Vnhp(a,z_tilda)+Uh(q),Vthpn(a,q,z_tilda)-cost_p,Vthpm(a,z_tilda,q)];
                        gt_meet_t_nm(a_tilda,z_tilda,q_tilda,a,q) = max(max(nu) - Vnh(a,q) - Vth(a_tilda,z_tilda,q_tilda) + Vmh(a_tilda,z_tilda), tol);
                    end
                end
            end
        end
    end
    
    % Final Gains from trade from meeting team

    gt_meet_t=max(gt_meet_t_m,gt_meet_t_nm);
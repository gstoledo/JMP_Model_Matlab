% Compute all gains from trade from empty firms
function [gt_meet_u, gt_meet_m, gt_meet_nm, gt_meet_t_m, gt_meet_t_nm, gt_meet_t] = ef_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,nm_penal)
    tol=1e-8;
    % Initialize gt_meet_u
    gt_meet_u = zeros(tpts,ats);
    % Compute gt_meet_u
    %gt_meet_u(z_tilda,a)
    for a = 1:ats
        for z_tilda = 1:tpts
            gt_meet_u(z_tilda,a) = max(max(Vmh(a,z_tilda),Vnh(a,z_tilda)) - Veh(a) -Uh(z_tilda), tol);
        end
    end
    % Initialize gt_meet_m
    gt_meet_m = zeros(ats,tpts,ats);
    % Compute gt_meet_m
    %gt_meet_m(a_tilda,z_tilda,a)
    for a = 1:ats
        for z_tilda = 1:tpts
            for a_tilda = 1:ats
                gt_meet_m(a_tilda,z_tilda,a) = max(max(Vmh(a,z_tilda),Vnh(a,z_tilda)) - Veh(a) - Vmh(a_tilda,z_tilda) + Veh(a_tilda), tol);
            end
        end
    end
    % Initialize gt_meet_nm
    gt_meet_nm = zeros(ats,tpts,ats);
    % Compute gt_meet_nm
    %gt_meet_nm(a_tilda,z_tilda,a)
    Vmhp=Non_man_penalty_loop(Vmh,tpts,ats,nm_penal);
    Vnhp=Non_man_penalty_loop(Vnh,tpts,ats,nm_penal);
    for a = 1:ats
        for z_tilda = 1:tpts
            for a_tilda = 1:ats
                gt_meet_nm(a_tilda,z_tilda,a) = max(max(Vmhp(a,z_tilda),Vnhp(a,z_tilda)) - Veh(a) - Vnh(a_tilda,z_tilda) + Veh(a_tilda), tol);
            end
        end
    end
    %Gains from trade from meeting firm with team
    % %Getting the manager
    %gt_meet_t_m(a_tilda,z_tilda,q_tilda,a)
    gt_meet_t_m= zeros(ats,tpts,tpts,ats);
    for a = 1:ats
        for q_tilda=1:tpts
            for z_tilda = 1:tpts
                for a_tilda = 1:ats
                    gt_meet_t_m(a_tilda,z_tilda,q_tilda,a) = max(max(Vmh(a,z_tilda),Vnh(a,z_tilda)) - Veh(a)...
                        - Vth(a_tilda,z_tilda,q_tilda) + Vnh(a_tilda,q_tilda), tol);
                end
            end
        end
    end
    %Getting the non-manager
    %gt_meet_t_nm(a_tilda,z_tilda,q_tilda,a)
    gt_meet_t_nm= zeros(ats,tpts,tpts,ats);
    for a = 1:ats
        for q_tilda=1:tpts
            for z_tilda = 1:tpts
                for a_tilda = 1:ats
                    gt_meet_t_nm(a_tilda,z_tilda,q_tilda,a) = max(max(Vmhp(a,q_tilda),Vnhp(a,q_tilda)) - Veh(a)...
                        - Vth(a_tilda,z_tilda,q_tilda) + Vmh(a_tilda,z_tilda), tol);
                end
            end
        end
    end     
    % Final Gains from trade from meeting team
    %gt_meet_t(a_tilda,z_tilda,q_tilda,a)
    gt_meet_t= max(gt_meet_t_m,gt_meet_t_nm);
end
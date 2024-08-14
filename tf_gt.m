% Compute all gains from trade from firms with team
function [gt_meet_u, gt_meet_m, gt_meet_nm, gt_meet_t_m, gt_meet_t_nm, gt_meet_t] = tf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,true,cd,cp)
    %%% Firm with team continuation value Vt(a,z,q)
    tol=1e-8;
    %Gains from trade from meeting unemployed
    gt_meet_u = zeros(tpts,ats,tpts,tpts);
    %gt_meet_u(z_tilda,a,z,q)
    for q = 1:tpts
        for z=1:tpts
            for a = 1:ats
                for z_tilda = 1:tpts
                    nu=[Vth(a,z_tilda,q)+Uh(z), Vth(a,z_tilda,z)+Uh(q)-cd, Vth(a,z,z_tilda)+Uh(q),Vth(a,q,z_tilda)+Uh(z)-cp];
                    gt_meet_u(z_tilda,a,z,q) = max(max(nu) - Vth(a,z,q)   -Uh(z_tilda),tol);
                end
            end
        end
    end
    
    %Gains from trade from meeting firm with manager
    gt_meet_m = zeros(ats,tpts,ats,tpts,tpts);
    %gt_meet_m(a_tilda,z_tilda,a,z,q)
    for q = 1:tpts
        for z=1:tpts
            for a = 1:ats
                for z_tilda = 1:tpts
                    for a_tilda = 1:ats
                        nu=[Vth(a,z_tilda,q)+Uh(z), Vth(a,z_tilda,z)+Uh(q)-cd, Vth(a,z,z_tilda)+Uh(q),Vth(a,q,z_tilda)+Uh(z)-cp];
                        gt_meet_m(a_tilda,z_tilda,a,z,q) = max(max(nu) - Vth(a,z,q) - Vmh(a_tilda,z_tilda) + Veh(a_tilda), tol);
                    end
                end
            end
        end
    end
    
    %Gains from trade from meeting firm with non-manager
    gt_meet_nm = zeros(ats,tpts,ats,tpts,tpts);
    %gt_meet_nm(a_tilda,z_tilda,a,z,q)
    %Vmhp=Non_man_penalty_loop(Vmh,tpts,ats,true);
    %Vnhp=Non_man_penalty_loop(Vnh,tpts,ats,true);
    Vthpm=Non_man_penalty_full_m_loop(Vth,tpts,ats,true);
    Vthpn=Non_man_penalty_full_n_loop(Vth,tpts,ats,true);
    for q = 1:tpts
        for z=1:tpts
            for a = 1:ats
                for z_tilda = 1:tpts
                    for a_tilda = 1:ats
                        nu=[Vthpm(a,z_tilda,q)+Uh(z), Vthpn(a,z,z_tilda)+Uh(q), Vthpm(a,z_tilda,z)+Uh(q)-cd, Vthpn(a,q,z_tilda)+Uh(z)-cp];
                        gt_meet_nm(a_tilda,z_tilda,a,z,q) = max(max(nu) - Vth(a,z,q) - Vnh(a_tilda,z_tilda) + Veh(a_tilda), tol);
                    end
                end
            end
        end
    end
    
    %Gains from trade from meeting firm with team
    
    %Grabbing the manager
    gt_meet_t_m= zeros(ats,tpts,tpts,ats,tpts,tpts);
    for q = 1:tpts
        for z=1:tpts
            for a = 1:ats
                for q_tilda = 1:tpts
                    for z_tilda = 1:tpts
                        for a_tilda = 1:ats
                            nu=[Vth(a,z_tilda,q)+Uh(z), Vth(a,z_tilda,z)+Uh(q)-cd, Vth(a,z,z_tilda)+Uh(q),Vth(a,q,z_tilda)+Uh(z)-cp];
                            gt_meet_t_m(a_tilda,z_tilda,q_tilda,a,z,q) = max(max(nu) - Vth(a,z,q) - Vth(a_tilda,z_tilda,q_tilda) + Vnh(a_tilda,q_tilda), tol);
                        end
                    end
                end
            end
        end
    end
    
    %Grabbing the non-manager
    gt_meet_t_nm= zeros(ats,tpts,tpts,ats,tpts,tpts);
    for q = 1:tpts
        for z=1:tpts
            for a = 1:ats
                for q_tilda = 1:tpts
                    for z_tilda = 1:tpts
                        for a_tilda = 1:ats
                            nu=[Vthpm(a,z_tilda,q)+Uh(z), Vthpn(a,z,z_tilda)+Uh(q), Vthpm(a,z_tilda,z)+Uh(q)-cd, Vthpn(a,q,z_tilda)+Uh(z)-cp];
                            gt_meet_t_nm(a_tilda,z_tilda,q_tilda,a,z,q) = max(max(nu) - Vth(a,z,q) - Vth(a_tilda,z_tilda,q_tilda) + Vmh(a_tilda,z_tilda), tol);
                        end
                    end
                end
            end
        end
    end
    
    %Final Gains from trade from meeting team
    gt_meet_t=max(gt_meet_t_m,gt_meet_t_nm);

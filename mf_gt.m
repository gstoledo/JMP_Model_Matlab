% Compute all gains from trade from firms with managers
function [gt_meet_u, gt_meet_m, gt_meet_nm, gt_meet_t_m, gt_meet_t_nm, gt_meet_t] = mf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,true,cd)
    %%% Firm with manager continuation value Vm(a,z)

    %Gains from trade from meeting unemployed
    %gt_meet_u(z_tilda,a,z)
    gt_meet_u= zeros(tpts,ats,tpts);
    for z = 1:tpts
        for a = 1:ats
            for z_tilda = 1:tpts
                nu=[Vmh(a,z_tilda)+Uh(z),Vth(a,z,z_tilda),Vth(a,z_tilda,z)-cd];
                gt_meet_u(z_tilda,a,z) = max(max(nu) - Vmh(a,z)  -Uh(z_tilda),0);  
            end
        end
    end
    
    
    
    %Gains from trade from meeting firm with manager
    %gt_meet_m(a_tilda,z_tilda,a,z)
    gt_meet_m= zeros(ats,tpts,ats,tpts);
    for z = 1:tpts
        for a = 1:ats
            for z_tilda = 1:tpts
                for a_tilda = 1:ats
                    nu=[Vmh(a,z_tilda)+Uh(z),Vth(a,z,z_tilda),Vth(a,z_tilda,z)-cd];
                    gt_meet_m(a_tilda,z_tilda,a,z) = max(max(nu) - Vmh(a,z) - Vmh(a_tilda,z_tilda) + Veh(a_tilda), 0);
                end
            end
        end
    end
    
    
    
    %Gains from trade from meeting firm with non-manager
    %gt_meet_n(a_tilda,z_tilda,a,z)
    Vmhp=Non_man_penalty_loop(Vmh,tpts,ats,true);
    %Vnhp=Non_man_penalty_loop(Vnh,tpts,ats,true);
    Vthpm=Non_man_penalty_full_m_loop(Vth,tpts,ats,true);
    Vthpn=Non_man_penalty_full_n_loop(Vth,tpts,ats,true);
    gt_meet_nm= zeros(ats,tpts,ats,tpts);
    for z = 1:tpts
        for a = 1:ats
            for z_tilda = 1:tpts
                for a_tilda = 1:ats
                    nu=[Vmhp(a,z_tilda)+Uh(z),Vthpn(a,z,z_tilda),Vthpm(a,z_tilda,z)-cd];
                    gt_meet_nm(a_tilda,z_tilda,a,z) = max(max(nu) - Vmh(a,z) - Vmh(a_tilda,z_tilda) + Veh(a_tilda), 0);
                end
            end
        end
    end        
    
    
    %Gains from trade from meeting firm with team
    %gt_meet_t(a_tilda,z_tilda,q_tilda,a,z)
    
    %Grabbing the manager
    gt_meet_t_m= zeros(ats,tpts,tpts,ats,tpts);
    for z = 1:tpts
        for a = 1:ats
            for q_tilda = 1:tpts
                for z_tilda = 1:tpts
                    for a_tilda = 1:ats
                        nu=[Vmh(a,z_tilda)+Uh(z),Vth(a,z,z_tilda),Vth(a,z_tilda,z)-cd];
                        gt_meet_t_m(a_tilda,z_tilda,q_tilda,a,z) = max(max(nu) - Vmh(a,z) - Vth(a_tilda,z_tilda,q_tilda) + Vnh(a_tilda,q_tilda), 0);
                    end
                end
            end
        end
    end
    
    %Grabbing the non-manager
    gt_meet_t_nm= zeros(ats,tpts,tpts,ats,tpts);
    for z = 1:tpts
        for a = 1:ats
            for q_tilda = 1:tpts
                for z_tilda = 1:tpts
                    for a_tilda = 1:ats
                        nu=[Vmhp(a,z_tilda)+Uh(z),Vthpn(a,z,z_tilda),Vthpm(a,z_tilda,z)-cd];
                        gt_meet_t_nm(a_tilda,z_tilda,q_tilda,a,z) = max(max(nu) - Vmh(a,z) - Vth(a_tilda,z_tilda,q_tilda) + Vmh(a_tilda,z_tilda), 0);
                    end
                end
            end
        end
    end
    
    % Final Gains from trade from meeting team
    %gt_meet_t(a_tilda,z_tilda,q_tilda,a,z)
    gt_meet_t= max(gt_meet_t_m,gt_meet_t_nm);
%LOM interation for SS 
function [e_udist, e_edist, e_mdist, e_ndist, e_tdist,failed]=lom_iterationNCC(eini_udist,eini_edist,eini_mdist,eini_ndist,eini_tdist,ats,tpts,phim,phin,Veh,Vmh,Vnh,U,Vth,Ve,Vm,Vn,Vt,cost_d,cost_p,lamu, lam, del, death, n,u_trans, a_trans,q_trans,typebirth,it_joint,speed_dist,display_iter_dist,nm_penal)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Lets do the iteration of the LOMs
    %Initial guesses for the LOM 
    eplus_udist=eini_udist;
    eplus_edist=eini_edist;
    eplus_mdist=eini_mdist;
    eplus_ndist=eini_ndist;
    eplus_tdist=eini_tdist;

    
    %Iterate on masses of workers
    diff_dist=1;
    diff_dist_max=1e-5;
    it_dist=0;
     
    %Speed of convergence for distributions
    % speed_dist=1;
    
    it_dist_max=20000;
    if it_joint>30 && it_joint<50
        min_it=5;
    elseif it_joint>=50
        min_it=1;
    else
        min_it=20;
    end
    failed=0; %Flag for failed convergence

    while (diff_dist>diff_dist_max  | it_dist<min_it) && it_dist<it_dist_max
        it_dist=it_dist+1;
        if it_dist==1
            e_udist=eplus_udist;
            e_edist=eplus_edist;
            e_mdist=eplus_mdist;
            e_ndist=eplus_ndist;
            e_tdist=eplus_tdist;
        end
        if it_dist>1
            e_udist=(1-speed_dist)*e_udist+speed_dist*eplus_udist;
            e_edist=(1-speed_dist)*e_edist+speed_dist*eplus_edist;
            e_mdist=(1-speed_dist)*e_mdist+speed_dist*eplus_mdist;
            e_ndist=(1-speed_dist)*e_ndist+speed_dist*eplus_ndist;
            e_tdist=(1-speed_dist)*e_tdist+speed_dist*eplus_tdist;
    
            %u=sum(e_udist); %mass of unemployed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Distributions and LOMs 
            
            % Distributions after the SM (implied by hiring and allocation policies)   
            % e1 distributions
            [e1_udist,e1_edist, e1_mdist, e1_ndist,e1_tdist]=e1_distNCC_phi(Veh,Vmh,Vnh,U,Vth,ats,tpts,phim,phin,cost_d,cost_p,e_udist,e_edist,e_mdist,e_ndist,e_tdist,lamu,lam,n,del,nm_penal);
            
            %It is not really getting to one...
            % sum(e1_edist)+sum(e1_mdist,"all")+sum(e1_ndist,"all")+sum(e1_tdist,"all"); % We are probabily missing a flow somewhere, that we will check later 
            
            
            % Distributions before productions (reallocation policies)
            %e2 distributions
            [e2_udist,e2_edist,e2_mdist,e2_ndist,e2_tdist]=e2_dist(Ve,Vm,Vn,U,Vt,ats,tpts,cost_d,cost_p,e1_udist,e1_edist,e1_mdist,e1_ndist,e1_tdist);
            
            %Check sum
            % sum(e2_edist,"all")+sum(e2_mdist,"all")+sum(e2_ndist,"all")+sum(e2_tdist,"all"); %This is summing up to one if the fed in dist sums to 1 
            
            %Disrtribtuions after shocks 
            %e3 distributions
            [e3_udist,e3_edist,e3_mdist,e3_ndist,e3_tdist]=e3_dist(ats,tpts,e2_udist,e2_edist,e2_mdist,e2_ndist,e2_tdist,u_trans,a_trans,q_trans);
            
            %Check sum
            % sum(e3_edist)+sum(e3_mdist,"all")+sum(e3_ndist,"all")+sum(e3_tdist,"all"); %This is summing up to one if the fed in dist sums to 1
            
            %Dist after entry and exit
            %eplus Distributions
            [eplus_udist,eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist]=eplus_dist(ats,tpts,e3_udist,e3_edist,e3_mdist,e3_ndist,e3_tdist,death,typebirth);
            
            %Check sum
            % nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all"); %This is summing up to one if the fed in dist sums to 1
            % popplus=sum(eplus_udist)+ sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");
            % store_n(it_dist-1)=nplus;
            % store_p(it_dist-1)=popplus;


            %%Error
            diff_dist=max([max(abs(eplus_udist-e_udist)),max(abs(eplus_edist-e_edist)),max(abs(eplus_mdist-e_mdist),[],[1 2]),max(abs(eplus_ndist-e_ndist),[],[1 2]),max(abs(eplus_tdist-e_tdist),[],[1 2 3])]);
            if display_iter_dist==1
                %Print every 100 iterations
                if mod(it_dist,100)==0
                    fprintf('Dist Iteration %d, error %f \n', it_dist, diff_dist)
                end 
                if diff_dist<diff_dist_max && it_dist>=min_it
                    cprintf('green','Distributions Converged in %d iterations\n',it_dist)
                end
            end
            %In red failed to converge
            if it_dist==it_dist_max
                fprintf(2,'LOM Failed to converge\n')
                failed=1;
            end

        end
    end
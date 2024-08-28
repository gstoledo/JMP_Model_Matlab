function [status, co_type, tenure, move, wage, co_tenure, co_wage]=FirmMeetsWorker( worker_state, worker_type, tie_break,  ...
        cdf_worker_state, cdf_worker_unempl, cdf_solo, cdf_team, ...
        type0, status0, co_type0, tenure0, co_tenure0, wage0, co_wage0,...
        gamma,nh,nquantiles,h0u_k, h1u_ik, h2u_ijk, r_ijk, h0e_k, h1e_ik, h2e_ijk, h0ee_kl,h1ee_ikl, h2ee_ijkl, wage_w,W_wkl,W_wk, U_k,Vmarginal_k, Vmarginal_ki, Vmarginal_kij )
 


        if ( status0 == 1 )   % worker is solo
            % WHO does the firm MEET?
            if ( worker_state <= cdf_worker_state(1) )  % meet an unemployed worker
                i=draw_CDF_1d( cdf_worker_unempl, worker_type );
                if ( h1u_ik(type0, i ) > 0.00 )  % firm hires unemployed worker to join team 
                    % write(*,*) "-firm hires unemployed worker of type", i
                    status = 2;
                    tenure = tenure0 + 1;
                    move = 0;
                    co_type = i;
                    co_tenure = 1;
                    % calculate wage for new hire 
                    target_value = (1.00 - gamma) * U_k(i) + gamma * Vmarginal_ki(i,type0);
                    co_wage=InterpolateWage( W_wkl(:,i,type0), target_value, wage_w );
                    wage = wage0;
                else % firm does not hire contacted worker_type
                    % write(*,*) "-firm does not hire unemployed worker of type ", i
                    status = 1;
                    tenure = tenure0 + 1;
                    move = 0;
                    co_type = 0;
                    co_tenure = 0;
                    co_wage = 0.00;
                    wage = wage0;
                end
            elseif ( (worker_state > cdf_worker_state(1)) & (worker_state <= cdf_worker_state(2)) ) 
                % meet a solo worker 
                i=draw_CDF_1d( cdf_solo, worker_type );
                if ( h1e_ik(type0, i) > 0.00 )  % firm hires solo worker to join team 
                    % write(*,*) "-firm hires solo worker of type", i
                    status = 2;
                    tenure = tenure0 + 1;
                    move = 0;
                    co_type = i;
                    co_tenure = 1;
                    % calculate wage for new hire 
                    target_value = (1.00 - gamma) * Vmarginal_k(i) + gamma * Vmarginal_ki(i,type0);
                    co_wage=InterpolateWage( W_wkl(:,i,type0), target_value, wage_w );
                    % update worker's wage in case it is outside bargaining set (here we use monotonicity of worekr value function in term of the wage, given match type)
                    wage = wage0;
                else % firm does not hire contacted worker_type
                    % write(*,*) "-firm does not hire solo worker of type ", i
                    status = 1;
                    tenure = tenure0 + 1;
                    move = 0;
                    co_type = 0;
                    co_tenure = 0;
                    co_wage = 0.00;
                    wage = wage0;
                end
            else % firm meets a worker in a team
                [i, j] = draw_CDF_2d( cdf_team, nh, worker_type );
                if ( h1ee_ikl(type0, i, j) > 0.00 )  % firm poaches worker from team 
                    % write(*,*) "-firm hires worker from a team: ", i, j
                    status = 2;
                    tenure = tenure0 + 1;
                    move = 0;
                    co_type = i;
                    co_tenure = 1;
                    % calculate wage for new hire 
                    target_value = (1.00 - gamma) * Vmarginal_ki(i,j) + gamma * Vmarginal_ki(i,type0);
                    co_wage=InterpolateWage( W_wkl(:,i,type0), target_value, wage_w );
                    % update worker's wage in case it is outside bargaining set (here we use monotonicity of worekr value function in term of the wage, given match type)
                    wage = wage0;
                else % firm does not hire contacted worker_type
                    % write(*,*) "-firm does not hire worker from a team: ", i, j
                    status = 1;
                    tenure = tenure0 + 1;
                    move = 0;
                    co_type = 0;
                    co_tenure = 0;
                    co_wage = 0.00;
                    wage = wage0;
                end
            end 

        else  % worker is in a team
            % WHO does the firm MEET?
            if ( worker_state <= cdf_worker_state(1) )  % meet an unemployed worker
                i=draw_CDF_1d( cdf_worker_unempl, worker_type );
                if ( h2u_ijk(type0, co_type0, i ) > 0.00 )  % firm hires unemployed worker to join team
                    % write(*,*) "-firm hires unemployed worker of type", i
                    % NOTE: the unifrom random numbers are in [0,1). If r_ijk = 0 this always evaluates to 0; if r_ijk = 1
                    %       this always evaluates to 1; if r_ijk = 0.5 it evaluates to 0 or 1 with equal probability.
                    if ( tie_break <  r_ijk(type0, co_type0, i) )  % worker is replaced
                        % write(*,*) "--worker is replaced"
                        status = 0;
                        tenure = 1;
                        move = 1;
                        co_type = 0;
                        co_tenure = 0;
                        wage = 0.00;
                        co_wage = 0.00;
                    else
                        % write(*,*) "--co-worker is replaced"
                        status = 2;           % co-worker is replaced 
                        tenure = tenure0 + 1;
                        move = 0;
                        co_type = i;
                        co_tenure = 1;
                        % calculate wage for new hire 
                        target_value = (1.00 - gamma) * U_k(i) + gamma * Vmarginal_kij(i,type0,co_type0);
                        co_wage=InterpolateWage( W_wkl(:,i,type0), target_value, wage_w );
                        % update worker's wage in case it is outside bargaining set (here we use monotonicity of worekr value function in term of the wage, given match type)
                        wage = wage0;
                    end
                else % firm does not hire contacted worker_type
                    % write(*,*) "-firm does not hire unemployed worker of type", i
                    status = 2;
                    tenure = tenure0 + 1;
                    move = 0;
                    co_type = co_type0;
                    co_tenure = co_tenure0 + 1;
                    co_wage = co_wage0;
                    wage = wage0;
                end
            elseif ( (worker_state > cdf_worker_state(1)) & (worker_state <= cdf_worker_state(2)) )  % meet a solo worker 
                i=draw_CDF_1d( cdf_solo, worker_type );
                if ( h2e_ijk(type0, co_type0, i ) > 0.00 )  % firm hires solo worker to join team 
                    % write(*,*) "-firm hires solo worker of type", i
                    if ( tie_break <=  r_ijk(type0, co_type0, i) )  % worker is replaced
                        % write(*,*) "--worker is replaced"
                        status = 0;
                        tenure = 1;
                        move = 1;
                        co_type = 0;
                        co_tenure = 0;
                        wage = 0.00;
                        co_wage = 0.00;
                    else
                        % write(*,*) "--co-worker is replaced"
                        status = 2;
                        tenure = tenure0 + 1;
                        move = 0;
                        co_type = i;
                        co_tenure = 1;
                        % calculate wage for new hire 
                        target_value = (1.00 - gamma) * Vmarginal_k(i) + gamma * Vmarginal_kij(i,type0,co_type0);
                        co_wage=InterpolateWage( W_wkl(:,i,type0), target_value, wage_w );
                        wage = wage0;
                    end
                else % firm does not hire contacted worker_type
                    % write(*,*) "-firm does not hire solo worker of type ", i
                    status = 2;
                    tenure = tenure0 + 1;
                    move = 0;
                    co_type = co_type0;
                    co_tenure = co_tenure0 + 1;
                    co_wage = co_wage0;
                    wage = wage0;
                end
            else % firms meets a worker in a team 
                [i, j]=draw_CDF_2d( cdf_team, nh, worker_type ); 
                if ( h2ee_ijkl(type0, co_type0, i, j) > 0.00 )  % firm poaches the i worker from (i,j) team 
                    % write(*,*) "-firm hires worker from a team: ", i, j
                    if ( tie_break <=  r_ijk(type0, co_type0, i) )  % worker is replaced
                        % write(*,*) "--worker is replaced"
                        status = 0;
                        tenure = 1;
                        move = 1;
                        co_type = 0;
                        co_tenure = 0;
                        wage = 0.00;
                        co_wage = 0.00;
                    else
                        % write(*,*) "--co-worker is replaced"
                        status = 2;
                        tenure = tenure0 + 1;
                        move = 0;
                        co_type = i;
                        co_tenure = 1;
                        % calculate wage for new hire 
                        target_value = (1.00 - gamma) * Vmarginal_ki(i,j) + gamma * Vmarginal_kij(i,type0,co_type0);
                        co_wage=InterpolateWage( W_wkl(:,i,type0), target_value, wage_w );
                        wage = wage0;
                    end
                else % firm does not hire contacted worker_type
                    % write(*,*) "-firm does not hire worker from a team: ", i, j
                    status = 2;
                    tenure = tenure0 + 1;
                    move = 0;
                    co_type = co_type0;
                    co_tenure = co_tenure0 + 1;
                    co_wage = co_wage0;
                    wage =  wage0;
                end
            end
        end

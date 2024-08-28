function [status, co_type, tenure, move, wage, co_tenure, co_wage]=WorkerMeetsFirm(  firm_state, firm_type, tie_break, cdf_firm_state, cdf_solo, cdf_team, ...
    wage_target_MP_solo, quantiles, type0, status0, co_type0, tenure0, co_tenure0, wage0, co_wage0,...
    gamma,nh,nquantiles,h0u_k, h1u_ik, h2u_ijk, r_ijk, h0e_k, h1e_ik, h2e_ijk, h0ee_kl,h1ee_ikl, h2ee_ijkl, wage_w,W_wkl,W_wk, U_k,Vmarginal_k, Vmarginal_ki, Vmarginal_kij )
 
if ( status0 == 0 )
    % WHO DO I MEET?
    if ( firm_state <= cdf_firm_state(1) )  % meet a vacant firm
        if ( h0u_k(type0) > 0.00 )  % move to new firm
            % write(*,*) "--worker joins vacant firm"
            status = 1;
            tenure = 1;
            move = 1;
            target_value = (1.00 - gamma) * U_k(type0) + gamma * Vmarginal_k(type0);
            wage=InterpolateWage( W_wk(:,type0) , target_value, wage_w );
            co_type = 0;
            co_tenure = 0;
            co_wage = 0.00 ;
        else % STAY UNEMPLOYED
            % write(*,*) "--worker does not join vacant firm"
            status = 0;
            tenure = tenure0 + 1;
            move = 0;
            wage = 0.00;
            co_type = 0;
            co_tenure = 0;
            co_wage = 0.00;
        end
    elseif (  ( firm_state > cdf_firm_state(1) ) & firm_state <= cdf_firm_state(2) )  % meet a solo firm
        j=draw_CDF_1d( cdf_solo, firm_type);
        if ( h1u_ik(j,type0 ) > 0.00 )  % form new team
            % write(*,*) "--worker joins firm with worker ", j
            status = 2 ;
            co_type = j ;
            tenure = 1;
            move = 1;
            target_value = (1.00-gamma) * U_k(type0) + gamma * Vmarginal_ki(type0,j);
            wage=InterpolateWage( W_wkl(:,type0, co_type), target_value, wage_w);
            co_tenure = 1;  
            k = 1+floor( tie_break*nquantiles ) ;
            co_wage = quantiles( k ,co_type, nh+1 );
            
        else % Stay unemployed
            % write(*,*) "--worker does not join firm with worker ", j
            status = 0;
            tenure = tenure0 + 1;
            move = 0;
            wage = 0.00;
            co_type = 0;
            co_tenure = 0;
            co_wage = 0.00;
        end
    else % meet an existing (i,j) team
        [i, j]=draw_CDF_2d( cdf_team, nh, firm_type );
        if ( h2u_ijk(i,j,type0) > 0.00 )
            % write(*,*) "--worker joins firm with old team ", i, j
            status = 2;
            tenure = 1;
            move = 1;
            co_type = j ;                                  % needed to take care of boundry at 1.
            if ( tie_break < r_ijk(i,j,type0) )  % worker replaces i
                % write(*,*) "--worker new teammate ", j
                co_type = j;
                old_co_type = i ;
            else
                % write(*,*) "--worker new teammate ", i
                co_type = i;
                old_co_type = j;
            end
            
            target_value = (1.00-gamma) * U_k(type0) + gamma * Vmarginal_kij(type0,i,j);
            wage=InterpolateWage( W_wkl(:,type0, co_type), target_value, wage_w );
            co_tenure = 1;   
            k = 1+floor( tie_break*nquantiles );
            co_wage = quantiles( k, co_type, old_co_type );
            
        else % Stay unemployed
            % write(*,*) "--worker does not join firm with old team ", i, j
            status = 0;
            tenure = tenure0 + 1;
            move = 0;
            wage = 0.00;
            co_type = 0;
            co_tenure = 0;
            co_wage = 0.00;
        end
    end
elseif ( status0 == 1 )  % Worker is employed solo
    if ( firm_state <= cdf_firm_state(1) )  % meet a vacant firm
 
        % This event always bids worker to her marginal value and she stays in the current match.
        status = 1;
        tenure = tenure0 + 1;
        wage = wage_target_MP_solo(type0);
        move = 0;
        co_type = 0;
        co_tenure = 0;
        co_wage = 0.00;
        
    elseif ( (firm_state > cdf_firm_state(1)) & (firm_state <= cdf_firm_state(2) ) )  % meet solo-worker firm
        j=draw_CDF_1d( cdf_solo, firm_type);
        if ( h1e_ik(j,type0 ) > 0.00 )  % move and form new team
            % write(*,*) "--worker moves to join firm of type ", j
            status = 2 ;
            co_type = j ;
            tenure = 1;
            move = 1;
            target_value = (1.00-gamma) * Vmarginal_k(type0) + gamma * Vmarginal_ki(type0,j);
            wage=InterpolateWage( W_wkl(:,type0, co_type), target_value, wage_w );
            co_tenure = 1;   
            k = 1+floor( tie_break*nquantiles );
            co_wage = quantiles( k ,co_type, nh+1 );
            
        else % Stay in current firm
            % write(*,*) "--worker does not join firm of type ", j
            status = 1;
            tenure = tenure0 + 1;
            % check if worker uses poaching firm to bid up wage
            target_value = Vmarginal_ki(type0,j);  
            current_value =InterpolateWage( wage_w, wage0, W_wk(:,type0) );
            if ( target_value > current_value )  % worker renegotiates wage with current firm
                % write(*,*) "----worker renegotiates wage with current firm"
                wage=InterpolateWage( W_wk(:,type0), target_value, wage_w );
            else
                wage = wage0;
            end
            move = 0;
            co_type = 0;
            co_tenure = 0;
            co_wage = 0.00;
        end
    else % meet an existing (i,j) team
        [i, j]=draw_CDF_2d( cdf_team, nh, firm_type );
        if ( h2e_ijk(i,j,type0) > 0.00 )
            % write(*,*) "--worker moves to join firm team of type ", i, j
            status = 2;
            tenure = 1;
            move = 1;
            co_type = j;                                   % needed to take care of boundry at 1.
            if ( tie_break < r_ijk(i,j,type0) )  % worker replaces i
                % write(*,*) "----worker has new co-worker ", j
                co_type = j;
                old_co_type = i ;
            else
                % write(*,*) "----worker has new co-worker ", i
                co_type = i;
                old_co_type = j;
            end
            
            target_value = (1.00-gamma) * Vmarginal_k(type0) + gamma * Vmarginal_kij(type0,i,j) ;
            wage=InterpolateWage( W_wkl(:,type0, co_type), target_value, wage_w );
            co_tenure = 1;   
            k = 1+floor( tie_break*nquantiles );
            co_wage = quantiles( k ,co_type, old_co_type );
            
        else % Stay in current solo firm
            % write(*,*) "--worker does not moves to join firm team of type ", i, j
            status = 1;
            tenure = tenure0 + 1;
            move = 0;
            % check if worker uses poaching firm to bid up wage
            target_value = Vmarginal_kij(type0,i,j); 
            current_value=InterpolateWage( wage_w, wage0, W_wk(:,type0) );
            if ( target_value > current_value )  % worker renegotiates wage with current firm
                % write(*,*) "----worker renegotiates wage with current firm"
                wage=InterpolateWage( W_wk(:,type0), target_value, wage_w );
            else
                wage = wage0;
            end
            co_type = 0;
            co_tenure = 0;
            co_wage = 0.00;
        end
        
    end
    
else % worker is in a team with co_type0
    
    if ( firm_state <= cdf_firm_state(1) )  % meet a vacant firm
        if ( h0ee_kl(type0,co_type0) > 0.00 )  % move to new firm
            % write(*,*) "--worker joins vacant firm"
            status = 1;
            tenure = 1;
            move = 1;
            target_value = (1.00 - gamma) * Vmarginal_ki(type0,co_type0) + gamma * Vmarginal_k(type0);
            wage=InterpolateWage( W_wk(:,type0) , target_value, wage_w );
            co_type = 0;
            co_tenure = 0;
            co_wage = 0.00;
        else % Stay in current firm
            % write(*,*) "--worker does not joins vacant firm"
            status = 2;
            tenure = tenure0 + 1;
            move = 0;
            co_type = co_type0;
            co_tenure = co_tenure0 + 1;
            co_wage = co_wage0;
            % check if worker uses poaching firm to bid up wage
            target_value = Vmarginal_k(type0 ); 
            current_value=InterpolateWage( wage_w, wage0, W_wkl(:,type0,co_type0) );
            if ( target_value > current_value )  % worker renegotiates wage with current firm
                % write(*,*) "----worker renegotiates wage with current firm"
                wage=InterpolateWage( W_wkl(:,type0,co_type0), target_value, wage_w );
            else
                wage = wage0;
            end
        end
    elseif ( (firm_state > cdf_firm_state(1)) & (firm_state <= cdf_firm_state(2) ) )  % meet solo-worker firm
        j= draw_CDF_1d( cdf_solo, firm_type );
        if ( h1ee_ikl(j, type0, co_type0 ) > 0.00 )  % move and form new team
            % write(*,*) "--worker moves to join firm of type ", j
            status = 2 ;
            co_type = j ;
            tenure = 1;
            move = 1;
            target_value = (1.00-gamma) * Vmarginal_ki(type0, co_type0) + gamma * Vmarginal_ki(type0, j);
            wage=InterpolateWage( W_wkl(:,type0, j), target_value, wage_w );
            co_tenure = 1;   
            k = 1+floor( tie_break*nquantiles ) ;
            co_wage = quantiles( k ,co_type, nh+1 );
            
        else % Stay in current firm
            % write(*,*) "--worker does not move to join firm of type ", j
            status = 2;
            tenure = tenure0 + 1;
            move = 0;
            co_type = co_type0;
            co_tenure = co_tenure0 + 1;
            co_wage = co_wage0;
            % check if worker uses poaching firm to bid up wage
            target_value = Vmarginal_ki(type0, j );  
            current_value=InterpolateWage( wage_w, wage0, W_wkl(:,type0,co_type0) );
            if ( target_value > current_value )  % worker renegotiates wage with current firm
                % write(*,*) "----worker renegotiates wage with current firm"
                wage=InterpolateWage( W_wkl(:,type0,co_type0), target_value, wage_w );
            else
                wage = wage0;
            end
        end
    else % meet an existing (i,j) team
        [i, j]=draw_CDF_2d( cdf_team, nh, firm_type );
        if ( h2ee_ijkl(i,j,type0,co_type0) > 0.00 )
            % write(*,*) "--worker moves to join firm team of type ", i, j
            status = 2;
            tenure = 1;
            move = 1;
            if ( tie_break < r_ijk(i,j,type0) )  % worker replaces i
                % write(*,*) "----worker has new co-worker ", j
                co_type = j;
                old_co_type = i ;
            else
                % write(*,*) "----worker has new co-worker ", i
                co_type = i;
                old_co_type = j ;
            end
            
            target_value = (1.00-gamma) * Vmarginal_ki(type0,co_type0) + gamma * Vmarginal_kij(type0,i,j) ;
            wage=InterpolateWage( W_wkl(:,type0, co_type), target_value, wage_w );
            co_tenure = 1;    
            k = 1+floor( tie_break*nquantiles ) ;
            co_wage = quantiles( k ,co_type, old_co_type );
            
        else % Stay in current team
            % write(*,*) "--worker does not move to join firm team of type ", i, j
            status = 2;
            tenure = tenure0 + 1;
            move = 0;
            co_type = co_type0;
            co_tenure = co_tenure0 + 1;
            co_wage = co_wage0;
            % check if worker uses poaching firm to renegotiate wage
            target_value = Vmarginal_kij(type0,i,j) ; 
            current_value=InterpolateWage( wage_w, wage0, W_wkl(:,type0,co_type0) );
            if ( target_value > current_value )  % worker renegotiates wage with current firm
                % write(*,*) "----worker renegotiates wage with current firm"
                wage= InterpolateWage( W_wkl(:,type0,co_type0), target_value, wage_w );
            else
                wage = wage0;
            end
        end
        
    end
end

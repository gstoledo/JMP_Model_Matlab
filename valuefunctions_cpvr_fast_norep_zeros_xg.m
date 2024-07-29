function [V0, V1, V2, V1h, V2h, U ]=valuefunctions_cpvr_fast_norep_zeros_xg(udist,solodist,teamdist,V0ini,V1ini,V2ini,Uini,utot,etot,nfirm1,nfirm2,nfirm0,lam0, lam1, del, bt, death, bpf, bpw, nfirm , tpts, b, fsolo, fteam,   uplus, soloplus,   teamplus)
%Evolution of unemployed workers
%udistup is post-type dist of unemployed
%solodistup is post-type dist of solo workers
%matchup is post-type dist of teams

 
%Mass of unemployed, solo, and teams
udistc = udist;
solodistc = solodist;
teamdistc = teamdist;

%Use initial guesses for value functions
V0up = V0ini;
V1up = V1ini;
V2up = V2ini;
Uup = Uini;

diff=100;
diffmax=.5e-4;
 
it=0; 
itmax=1000;
 

zero_tpts=zeros(1,tpts);
zero_tpts_tpts=zeros(tpts,tpts);

Ub=zeros(1,tpts);
V1b=zeros(1,tpts);
V2b=zeros(tpts,tpts);

while diff>diffmax & it<itmax
 it=it+1;
 
    V0 = V0up; %Update step
    V1 = V1up;
    V2 = V2up;
    U = Uup;
 
    V1h =  max( V1, V0 + U ) ;
    
    %This is used to deal with the two possible dimensions of full teams 
    Urep=repmat(U,tpts,1);
    Uprimerep=repmat(U',1,tpts);
        
    V2h=max(V2, max(repmat(V1',1,tpts) + Urep, max(repmat(V1,tpts,1)+ Uprimerep, V0 + Uprimerep+Urep ) ) );
    

    %Unemployed value after type evolution 
    for kup=1:tpts
        
        V2hkup=V2h(kup, :); 
        Ukup=U(kup);
        V1hkup=V1h(kup);
        
        inner1= max( repmat(V2hkup',1,tpts) - V2h + Urep - Ukup, 0);      % Replace first guy
        inner2= max( repmat(V2hkup,tpts,1) - V2h + Uprimerep - Ukup, 0);  % Replace second guy
        
        % The marginal values of hiring from U to different types of firms
        ev_meet_vacant=max(  V1hkup - V0 - Ukup, 0  );
        ev_meet_solo=max(V2hkup - V1h - Ukup, 0); %Dont have the problem of firing this guy here 
        ev_meet_team=max(inner1 ,inner2   );

        Ub(kup)= bt*(death*(0 - Ukup) +...
            lam0*(nfirm0/nfirm)*bpw.*ev_meet_vacant + ...
            lam0/nfirm*bpw.*sum(ev_meet_solo.*solodistc)+...
            lam0/(2*nfirm)*bpw.*sum(sum(ev_meet_team.*teamdistc)) +...
            Ukup);

    end
    

    %Unemployed home production and cont value

    for k=1:tpts 

        ev_u=Ub(kminus(k): kplus(k)).*uplus(k,kminus(k): kplus(k));
        Uup(k) =  b(k) + sum(ev_u);
        
    end
    
    %Idle firm cont value
 
    %Expectations over meeting unempl, solo, ith member of team (i,j)
 
    inner3=repmat(V1h',1,tpts) - V2h + repmat(V1h,tpts,1) - V0; %Idle firm stealing from full firm
    ev_firm_meet_u=max(V1h - U - V0, zero_tpts); 
    ev_firm_meet_team=max(inner3, zero_tpts_tpts);
    
    %Value function of idle firm
    
    V0up =  0 + bt*((lam0/nfirm)*bpf*sum(ev_firm_meet_u.*udistc) + (lam1/nfirm)*bpf*sum(sum(ev_firm_meet_team.*teamdistc)) + V0);
    
    %Value of solo firm
    for kup=1:tpts
        
        V2hkup=V2h(kup, :); %Value of aquiring a type kup
        Ukup=U(kup);
        V1hkup=V1h(kup);        
         
        ev_meet_vacant=0;   
        
        
        inner4=max( repmat(V2hkup',1,tpts) - V2h + Urep  - V1hkup + V0, zero_tpts_tpts);
        inner5=max(repmat(V2hkup,tpts,1) - V2h + Uprimerep - V1hkup + V0, zero_tpts_tpts);
        ev_meet_solo=max(V2hkup - V1h - V1hkup+V0, zero_tpts); %Worker meets solo
        
        ev_firm_meet_u=max( V2hkup - U - V1hkup, zero_tpts ); %Firm meets unempl
        ev_firm_meet_solo=max( V2hkup - V1h + V0 - V1hkup, zero_tpts ); %Firm meets solo
        

        ev_meet_team= max(inner4, inner5 ); %Worker meets team
                
        inner6=repmat(V2hkup',1,tpts) - V2h + repmat(V1h,tpts,1) - V1hkup;
        ev_firm_meet_team= max(inner6, zero_tpts_tpts);        
                
        
        
        
        %Cont value
        V1b(kup) =   bt*(death*(V0 - V1hkup) + del*(Ukup + V0 - V1hkup) +...
            ev_meet_vacant +...  
            lam1/nfirm*bpw.*sum(ev_meet_solo.*solodistc) + ...
            lam1/(2*nfirm)*bpw.*sum(sum(ev_meet_team.*teamdistc)) +...
            (lam0/nfirm)*bpf.*sum(ev_firm_meet_u.*udistc) +...
            (lam1/nfirm)*bpf.*sum(ev_firm_meet_solo.*solodistc) +...
            (lam1/nfirm)*bpf.*sum(sum( ev_firm_meet_team.*teamdistc)) +  V1hkup);
        
    end
    
    
    %Cont value of solo
    for k=1:tpts

        ev_solo=V1b(kminus(k):kplus(k)).*soloplus(k, kminus(k):kplus(k));
        V1up(k) =  fsolo(k) + sum(ev_solo);
        
    end
    
    
    %Cont value of team
    for kup=1:tpts
        
        V2hkup=V2h(kup, :);
        Ukup=U(kup);
        V1hkup=V1h(kup);        
        for lup=1:tpts
            
            
        	V2hlup=V2h(lup, :);
            V2hkuplup=V2h(kup, lup);
            V1hlup=V1h(lup);
            Ulup=U(lup); 
            
            
            ev_meet_solo=max(V2hkup - V2hkuplup - V1h + V1hlup ,  zero_tpts);
            ev_co_meet_solo=max(V2hlup - V2hkuplup - V1h + V1hkup, zero_tpts);
            
            inner7=max(V2hlup - V2hkuplup + Ukup - U ,   zero_tpts) ;
            ev_firm_meet_u=  max(V2hkup - V2hkuplup + Ulup - U , inner7 );
            
            
            inner8=max(V2hlup -  V2hkuplup + Ukup - V1h + V0 , zero_tpts);
            ev_firm_meet_solo=max(V2hkup - V2hkuplup + Ulup - V1h + V0, inner8   );
                
                
            inner9=repmat(V2hkup',1,tpts) - V2hkuplup - V2h + Urep  + V1hlup;
            inner10=max(repmat(V2hkup,tpts,1) -  V2hkuplup - V2h + Uprimerep + V1hlup, zero_tpts_tpts);
            ev_meet_team= max(inner9,inner10 );
            
            inner11=repmat(V2hlup',1,tpts) - V2hkuplup - V2h + Urep  + V1hkup ;
            inner12=max(repmat(V2hlup,tpts,1) - V2hkuplup - V2h + Uprimerep + V1hkup, zero_tpts_tpts);
            ev_co_meet_team= max(inner11, inner12 );
            
            
            inner13=repmat(V2hkup',1,tpts) -  V2hkuplup - V2h  + repmat(V1h,tpts,1) + Ulup ;
            inner14=max(repmat(V2h(:, lup),1,tpts) -  V2hkuplup - V2h + repmat(V1h,tpts,1) + Ukup , zero_tpts_tpts);
            ev_firm_meet_team=  max(inner13,inner14 );            
                    

            ev_meet_vacant=max(V1hkup-V0  - V2hkuplup + V1hlup , 0);
            ev_co_meet_vacant=max(V1hlup-V0  - V2hkuplup + V1hkup, 0);
            
            V2b(kup,lup) =  bt*( death*(V1hlup - V2hkuplup) +...
                del*(V1hlup + Ukup - V2hkuplup) +...
                death*(V1hkup - V2hkuplup) +...
                del*(V1hkup + Ulup - V2hkuplup) +...
                lam1*nfirm0/nfirm*bpw.*ev_meet_vacant +  ...
                lam1/nfirm*bpw.*sum(ev_meet_solo.*solodistc) +   ...
                lam1/(2*nfirm)*bpw.*sum(sum(ev_meet_team.*teamdistc))+...
                lam1*nfirm0/nfirm*bpw.*ev_co_meet_vacant +   ...
                lam1/nfirm*bpw.*sum(ev_co_meet_solo.*solodistc)+   ...
                lam1/(2*nfirm)*bpw.*sum(sum(ev_co_meet_team.*teamdistc)) +  ...
                ((lam0)/nfirm)*bpf.*sum(ev_firm_meet_u.*udistc) +...
                ((lam1)/nfirm)*bpf.*sum(ev_firm_meet_solo.*solodistc) +...
                ( lam1/nfirm)*bpf.*sum(sum(ev_firm_meet_team.*teamdistc)) +...
                V2hkuplup );
        end
    end
    
    
    %Add in production
    for k=1:tpts
        for l=1:tpts
            
            ev_team=zeros(1,tpts);
      
            for kup=kminus(k):kplus(k)
                
                ev_team(kup)= teamplus(k, kup,l)*sum(V2b(kup,kminus(l):kplus(l)).*teamplus(l,kminus(l):kplus(l),k));
            end
            
            
            
            V2up(k,l) =    fteam(k, l) + sum(ev_team) ;
        end
    end
    
    
    diff=max([max(abs(V0 - V0up)),max(max(abs(V1 - V1up))),max(max(abs(V2 - V2up))),max(abs(U-Uup))]);
    
end %end loop over iterations
 

V0 = V0up;
V1 = V1up;
V2 = V2up;
U = Uup;

for i=1:tpts-1
    for j=i+1:tpts
        V2(j, i) = V2(i, j);
    end
end

 % Final hat VFs that the function will also spit out
V1h=max(V1,V0+U);

 
V2h=max( max(max(V2, repmat(V1',1,tpts) + Urep), repmat(V1,tpts,1) + Uprimerep),...
        V0 + Uprimerep+ Urep  );


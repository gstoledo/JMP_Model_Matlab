clear all
close all
%load('xmin_results_combined')
load('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit/xmin_results_combined.mat')
x=xmin;

%Creating all global variables
global lam0 lam1 del bt death bpf bpw nfirm  tpts spts typemin typemax...
    type typebirth fcomp homeprod soloprod b fsolo fteam ulose ugain ustay uplus...
    sololose sologain solostay soloplus   teamplus...
    zero_tol wgrid ftypetot ftype1 ftype2    wtypetot  wtypeu  wtypesolo  wtypeteam ...
    wtypetot_pdf team_learn_param_up_sym team_learn_param_up team_learn_param_down wtypetot_initial_pdf

%% Worker for storage of results
task = getCurrentTask;
if isempty(task)
    xxx=1;
else
    xxx=task.ID;
end

tic

%% Unpack
lam0                    =x(1);
lam1                    =x(2);
fcomp                   =x(3); %Complementarity in F
mubar                   =x(4);
del                     =x(5);
ulose                   =x(6);
sologain                =x(7);
team_learn_param_up  	=x(8);
mnew	                =x(9);
team_learn_param_up_sym	=x(10);
bpf                     =x(11);
homeprod                =x(12);
sololose                =0;
A                       =x(14);
team_learn_param_down   =0;
mnew_high               =floor(x(15)); 

%Toggles
dayy          ='1'   ; %Version of the code    
zero_tol      =1e-10 ; %Tolerance for zero
use_guess     ='n'   ;
 
%Fundamentals
bt   =1/(1.1^(1/12))   ; %beta, discount factor
death=1/(35*12)        ; %Probability agent dies  
bpw  =1-bpf            ; %bpw is bargaining power of worker
nfirm=1                ; %Measure of firms 
tpts =7                ; %type point space
spts =tpts+2           ; %state S for updating -- {u,0,{j}} -- dimension of that space is type+2

typemin=1; %Lowest type
typemax=mubar ; %Highest type

type=[typemin:(typemax-typemin)/(tpts-1):typemax]; %Types

typebirth=zeros(1,tpts);
for i=1:mnew_high
    typebirth(i)= mnew*exp(-type(i)*mnew)./sum(mnew*exp(-type(1:mnew_high)*mnew)); %Birth rate
end
 
soloprod=1    ;
pnew    =fcomp;
pold    =1    ;
 
fteam=zeros(tpts,tpts);
fsolo=zeros(1,tpts);
for i=1:tpts
    fsolo(i)= soloprod*((type(i).^(fcomp)+type(1).^(fcomp))^(1/fcomp))./(2^(1/fcomp)) ;
    for j=1:tpts
        fteam(i,j)= A*((type(i).^(fcomp)+type(j).^(fcomp))^(1/fcomp))./(2^(1/fcomp))  ; %Team production (CES)
    end
end

b=homeprod*fsolo ; %Type home production vector
 
ugain=0.00           ; %Probability unemployed move up

ustay=1-ulose-ugain  ; %Probability unemployed stay

uplus=create_trans(ulose,ustay,ugain,tpts);

solostay=1-sololose-sologain; %Probability unemployed stay

soloplus=create_trans(sololose,solostay,sologain,tpts);

%teamplus(i,i+,j)--> prob of going from i to i+ given coworker j
for j=1:tpts %Loop over coworker type
    
    for i=1  
        teamplus(i,:,j)=[teamlose(type(i),type(j))+teamstay(type(i),type(j)),teamgain(type(i),type(j)),zeros(1,tpts-(i+1))];
    end
    
    for i=2:tpts-1
        teamplus(i,:,j)=[zeros(1,i-2),teamlose(type(i),type(j)),teamstay(type(i),type(j)),teamgain(type(i),type(j)),zeros(1,tpts-(i+1))];
    end
    
    for i=tpts
        teamplus(i,:,j)=[zeros(1,i-2),teamlose(type(i),type(j)),teamstay(type(i),type(j))+teamgain(type(i),type(j))];
    end
end
  
%%
%VFI initialization
V0ini = 0; %Guess value functinos -- idle firm
V1ini = fsolo /(1 - bt) ; %Single worker firm
V2ini = fteam /(1 - bt) ; %Two worker firm
Uini =  b /(1 - bt) ; %Unemployed worker
 
udist    = repmat((1/spts)*1/tpts, 1, tpts) ; %Unempl dist across types, guess
solodist = repmat((1/spts)*1/tpts,1,tpts)   ; %Emply dist across types, guess
teamdist = repmat((1/spts)*1/tpts,tpts,tpts); %Team dist, guess

utot   = sum(udist)                         ; %Guess unempl distribution
etot   = sum(solodist) + sum(sum(teamdist)) ; %Guess empl dist
nfirm1 = sum(solodist)                      ; %Guess mass of single worker firms
nfirm2 = sum(sum(teamdist))/2               ; %Guess mass of teams
nfirm0 = nfirm - nfirm1 - nfirm2            ; %Guess mass of idle firms

%Load existing v-func guesses to save on time
if exist('vval.mat')==2 && use_guess(1)=='y'
    load('vval.mat','V0up','V1up','V2up','Uup','V1hup','V2hup','udist','solodist','teamdist','utot','etot','nfirm1','nfirm2','nfirm0')
else
    [V0up, V1up, V2up, V1hup, V2hup, Uup] = valuefunctions_cpvr_fast_norep_zeros_xg(udist,solodist,teamdist,V0ini,V1ini,V2ini,Uini,utot,etot,nfirm1,nfirm2,nfirm0,lam0, lam1, del, bt, death, bpf, bpw, nfirm , tpts, b, fsolo, fteam,   uplus, soloplus,   teamplus);
end
  
%%
%VFI

%Iterate on value functions and type distribution
diff_joint    =100;
diff_joint_max=1e-4; %Value func/distribution max tolerance
  
it_joint      =0;
it_joint_max  =2000;
 
V0=V0up; 
V1=V1up;
V2=V2up;
V1h=V1hup;
V2h=V2hup;
U=Uup;
 
update_speed_v=1;
update_speed=1;
 
teamdist_outer=teamdist;
solodist_outer=solodist;
udist_outer=udist;
 
diff_joint_lag1=0;
diff_joint_lag2=0;


 while (diff_joint>diff_joint_max | diff_joint_lag1>diff_joint_max  | diff_joint_lag2>diff_joint_max | it_joint<250 ) &&  it_joint<it_joint_max
     
    it_joint=it_joint+1;
     
    udist    =update_speed*udist   +(1-update_speed)*udist_outer; %Update guesses of value functions, slowly update to avoid cycling
    solodist =update_speed*solodist+(1-update_speed)*solodist_outer;
    teamdist =update_speed*teamdist+(1-update_speed)*teamdist_outer;
    typedist(1,:)= sum(teamdist.') + solodist + udist; %Final distribution
    
    utot = sum(udist);
    etot = sum(solodist) + sum(sum(teamdist));
    nfirm1 = sum(solodist);
    nfirm2 = sum(sum(teamdist))/2;
    nfirm0 = nfirm - nfirm1 - nfirm2;
    
    %Given distributions, update value functions
    [V0up, V1up, V2up, V1hup, V2hup, Uup] = valuefunctions_cpvr_fast_norep_zeros_xg(udist,solodist,teamdist,V0,V1,V2,U,utot,etot,nfirm1,nfirm2,nfirm0,lam0, lam1, del, bt, death, bpf, bpw, nfirm , tpts, b, fsolo, fteam,   uplus, soloplus,   teamplus);
    
 
    V0=update_speed_v*V0up+(1-update_speed_v)*V0; %Update 
    V1=update_speed_v*V1up+(1-update_speed_v)*V1;
    V2=update_speed_v*V2up+(1-update_speed_v)*V2;
    V1h=update_speed_v*V1hup+(1-update_speed_v)*V1h;
    V2h=update_speed_v*V2hup+(1-update_speed_v)*V2h;
    U=update_speed_v*Uup+(1-update_speed_v)*U;
    
    
    %Iterate on masses of workers
    diff_dist=1;
    diff_dist_max=.5e-4;
    it_dist=0;
    
 
    it_dist_max=50;
    
    typedist_outer=typedist; %Initial distribution
    teamdist_outer=teamdist;
    solodist_outer=solodist;
    udist_outer=udist;
    
 
    if it_joint>100 & it_joint<200
        min_it=5;
    elseif it_joint>=200
        min_it=1;
    else
        min_it=25;
    end
    
    while (diff_dist>diff_dist_max  | it_dist<min_it) && it_dist<it_dist_max
        
        it_dist=it_dist+1;
        typedist_old=typedist; %Store older type distribution
        
        udistup = zeros(1, tpts);
        solodistup = zeros( 1, tpts);
        matchup = zeros(tpts,tpts) ;
        
        %Solve for resulting distribution given guesses
        [udistup, solodistup, matchup,h0u_k,h1u_ik, h2u_ijk] =  uevolution_vec_norep_zeros_xg(udistup,solodistup,matchup,...
            U,V0,V1h,V2h,V1,V2,udist,solodist,teamdist,utot,etot,nfirm1,nfirm2,nfirm0,lam0, lam1, del,  nfirm,  tpts, zero_tol);
        [udistup, solodistup, matchup,h0e_k, h1e_ik, h2e_ijk] = soloevolution_vec_norep_zeros_xg(udistup,solodistup,matchup,...
            U,V0,V1h,V2h,V1,V2,udist,solodist,teamdist,utot,etot,nfirm1,nfirm2,nfirm0,lam0, lam1, del,  nfirm,  tpts, zero_tol);
        [udistup, solodistup, matchup,h0ee_kl,  h1ee_ikl,  h2ee_ijkl, r_ijk] = teamevolution_vec_norep_zeros_xg(udistup,solodistup,matchup,...
            U,V0,V1h,V2h,V1,V2,udist,solodist,teamdist,utot,etot,nfirm1,nfirm2,nfirm0,lam0, lam1, del,  nfirm,  tpts, zero_tol);
        [udistup, solodistup, matchup] = dismissals_zeros_xg(udistup,solodistup,matchup,...
            U,V0,V1h,V2h,V1,V2,udist,solodist,teamdist,utot,etot,nfirm1,nfirm2,nfirm0,lam0, lam1, del,  nfirm,  tpts, zero_tol);
        udistprod=udistup; %Production stage
        solodistprod=solodistup;
        teamdistprod=matchup;
        [udistup, solodistup, matchup] = learning_xg(udistup,solodistup,matchup,...
            U,V0,V1h,V2h,V1,V2,udist,solodist,teamdist,utot,etot,nfirm1,nfirm2,nfirm0, tpts , uplus ,soloplus ,teamplus);
        [udistup, solodistup, matchup] = birth_xg(udistup,solodistup,matchup,...
            U,V0,V1h,V2h,V1,V2,udist,solodist,teamdist,utot,etot,nfirm1,nfirm2,nfirm0,death,tpts,typebirth);
         
        %Update
        udist=udistup;
        solodist=solodistup;
        teamdist=matchup;
        
        %Total employment
        utot = sum(udist);
        etot = sum(solodist) + sum(sum(teamdist));
        nfirm1 = sum(solodist);
        nfirm2 = sum(sum(teamdist))/2;
        nfirm0 = nfirm - nfirm1 - nfirm2;
         
        typedist(1,:)= sum(teamdist.') + solodist + udist; %Final distribution
         
        diff_dist=max(abs(typedist-typedist_old));
         
    end
     
    overall_dist_diff=max(abs(typedist_outer-typedist)); %Initial vs. final
    
    diff_joint_lag2=diff_joint_lag1    ; %Lags
    diff_joint_lag1=diff_joint     ; %Lags
    diff_joint= max(  [overall_dist_diff ] );
     
    display(['joint iteration: ' num2str(it_joint) ': joint diff: ' num2str(diff_joint) ' : type vector: ' num2str(typedist)])
     
    %Identify cycling distributions
    
    if overall_dist_diff>.01 && it_joint>100
        cycle=1;
        break
    else
        cycle=0;
    end
    
 end
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lets undertand until here 
 %Wages
 
% wmin  =  -6*min(fteam(:)); 
% wmax = max(fteam(:));
% wpts = 100;
 
% wgrid = [wmin:(wmax - wmin)/(wpts - 1):wmax];
 
% W1up = zeros(tpts,wpts);
% W2up = zeros(tpts,tpts,wpts);
 
% diff_w    =100;
% diff_max_w=1e-4;
% it_w      =0;
% it_max_w  =3000;
 
% Urep=repmat(U,tpts,1);
% Uprimerep=repmat(U',1,tpts);
% V1hprimerep=repmat(V1h',1,tpts);
% V1hrep=repmat(V1h,tpts,1);

% %VFI for employed value functions with wages
% while diff_w>diff_max_w & it_w<it_max_w
    
%     it_w=it_w+1;
    
%     %Update value functions
%     W1 = W1up;
%     W2 = W2up;
    
%     W1b =zeros(tpts,wpts);
%     W1up =zeros(tpts,wpts);
%     W2b = zeros(tpts,tpts,wpts);
%     W2up =zeros(tpts,tpts,wpts);
    
 
%     for w=1:wpts
%         for kup=1:tpts
%             V2hkup=V2h(kup, :);
%             V1hkup=V1h(kup);
%             Ukup=U(kup);
%             W2kupw=W2(kup, :, w);
            
%             inner1=contv(W1(kup, w), V1hkup - V0, Ukup);
%             inner2=aucv_zeros( V1hkup - V0, V2hkup - V1h) ;
            
%             ev_meet_solo= (solodist/nfirm).* max( inner1, inner2 ) ;
            
%             inner3=max( repmat(V2hkup',1,tpts) + Urep - V2h,   repmat(V2hkup,tpts,1) + Uprimerep - V2h  );
%             ev_meet_team=(teamdist/(2*nfirm)).*max( inner1, ...
%                 aucv_zeros( V1hkup - V0,  inner3 )  ) ;
            
%             I1=(V2hkup - V1hkup -U>zero_tol );
%             inner4=contv(W2kupw, V2hkup - V1h  , Ukup );
%             ev_firm_meets_u= (udist/nfirm).*(I1.*inner4 ...
%                 + (~I1).*inner1 ) ;
            
%             I2=(V2hkup - V1hkup -(V1h - V0)> zero_tol);
            
%             ev_firm_meets_solo= (solodist/nfirm).*( I2.*inner4 ...
%                 + (~I2)*inner1 ) ;
            
%             I3=(repmat(V2hkup',1,tpts) - V1hkup -(V2h - V1hrep)> zero_tol );
%             inner5=contv(repmat(W2kupw',1,tpts), repmat(V2hkup',1,tpts) - V1hprimerep , Ukup);
%             ev_firm_meets_team= (teamdist/nfirm).*( ...
%                 I3.*inner5 + ...
%                 (~I3).*inner1 ) ;
            
%             ev_meet_vacant=lam1*nfirm0/nfirm*max([ inner1, ...
%                 aucv_zeros(V1hkup - V0, V1hkup - V0) ]);
            
%             W1b(kup, w) =  bt*del*Ukup + ...
%                 bt*ev_meet_vacant  +  ...
%                 bt*lam1*sum(ev_meet_solo) + ...
%                 bt*lam1*sum(sum(ev_meet_team)) + ...  
%                 bt*lam0*sum(ev_firm_meets_u) + ...
%                 bt*lam1*sum(ev_firm_meets_solo) + ...
%                 bt*lam1*sum(sum(ev_firm_meets_team)) + ...  
%                 bt*(1 - del - death - lam1 - (lam0*utot + lam1*etot)/nfirm)*contv(W1(kup, w), V1hkup - V0 , Ukup) ;
%         end
%     end
    
 
%     for w=1:wpts
%         for k=1:tpts
            
%             W1up(k, w) =  wgrid(w) + sum(W1b(kminus(k):kplus(k), w).*(soloplus(k,kminus(k):kplus(k))'));
%         end
%     end
    
%     %Team worker, with wage w, own type kup, and coworker type lup
%     for w=1:wpts
%         for lup=1:tpts
%             for kup=1:tpts
%                 V2hkup=V2h(kup, :);
%                 V1hkup=V1h(kup);
%                 Ukup=U(kup);
%                 V2hlup=V2h(:, lup);
%                 V2hluprep=repmat(V2hlup,1,tpts);
%                 Ulup=U(lup);
%                 W2kupw=W2(kup, :, w);
%                 V2hkuplup=V2h(kup, lup);
                
%                 inner6=contv(W2(kup, lup, w), V2hkuplup - V1h(lup) , Ukup);
%                 inner7=aucv_zeros(V2hkuplup - V1h(lup),  V2hkup - V1h);
%                 ev_meet_solo= (solodist/nfirm) .* max(  inner6,  inner7 );
                
%                 inner8 =aucv_zeros(V2hkuplup - V1h(lup),    max(repmat(V2hkup',1,tpts) + Urep, repmat(V2hkup,tpts,1) + Uprimerep) -  V2h) ;
%                 ev_meet_team=(teamdist/(2*nfirm)).* max( inner6,   inner8    )   ;
                
%                 inner9=contv(W1(kup, w), V1hkup - V0 , Ukup);
%                 I4=(V2hlup' - V1h -(V2hkuplup - V1hkup)>zero_tol );
%                 ev_co_meets_solo= (solodist/nfirm).*(...
%                     I4.*inner9 ...
%                     + (~I4)*inner6  ) ;
                
%                 inner10=contv(W1(kup, w), V1hkup - V0 , Ukup);
%                 I5=(max(V2hluprep + Urep, V2hluprep' + Uprimerep ) -  V2h -(V2hkuplup - V1hkup)> zero_tol );
%                 ev_co_meets_team=(teamdist/(2*nfirm)).*(...
%                     I5.*inner10 ...
%                     +(~I5 )*inner6  ) ;
                
                
%                 I6=((max( V2hkup + Ulup, V2hlup' + Ukup  )  - U )-V2hkuplup>zero_tol );
%                 ev_firm_meets_u=(udist/nfirm).*( ...
%                     (~I6).*inner6 ...
%                     + (I6).*(  ...
%                     (V2hkup + Ulup -(V2hlup' + Ukup) > zero_tol).*contv(W2kupw, V2hkup - V1h , Ukup)...
%                     +(zero_tol< V2hlup' + Ukup -(V2hkup + Ulup)).* Ukup ...
%                     +(abs(V2hkup + Ulup -(V2hlup' + Ukup))<=zero_tol).* (.5*Ukup + .5*contv(W2kupw, V2hkup - V1h , Ukup) ) ...
%                     )  ) ;
                
                
%                 I7=((max( V2hkup + Ulup,  V2hlup' + Ukup ) - (V1h - V0)) -V2hkuplup>zero_tol   );
%                 ev_firm_meets_solo=(solodist/nfirm).*( ...
%                     (~I7).*inner6...
%                     + (I7).* ( ...
%                     (V2hkup + Ulup -(V2hlup' + Ukup)>zero_tol).*contv(W2kupw, V2hkup - V1h  , Ukup) ...
%                     + (zero_tol< V2hlup' + Ukup-(V2hkup + Ulup)).* Ukup ...
%                     + (abs(V2hkup + Ulup -(V2hlup' + Ukup))<=zero_tol).*( .5*Ukup + .5*contv(W2kupw, V2hkup - V1h , Ukup )  ) ...
%                     ) ) ;
                
                
                
%                 I8=((max(repmat(V2hkup',1,tpts) + Ulup, V2hluprep + Ukup ) - (V2h - V1hrep))-V2hkuplup >zero_tol  );
%                 ev_firm_meets_team=  teamdist/nfirm.* ( ...
%                     (~I8).*inner6 ...
%                     +   (I8).* ( ...
%                     (repmat(V2hkup',1,tpts) + Ulup -(V2hluprep + Ukup) >zero_tol ).*contv(repmat(W2kupw',1,tpts), repmat(V2hkup',1,tpts) - V1hprimerep , Ukup)  ...
%                     +( zero_tol< V2hluprep + Ukup-(repmat(V2hkup',1,tpts) + Ulup)).*Ukup  ...
%                     +( abs(repmat(V2hkup',1,tpts) + Ulup -( V2hluprep + Ukup))<=zero_tol).*(.5*Ukup + .5*contv(repmat(W2kupw',1,tpts), repmat(V2hkup',1,tpts) - V1hprimerep , Ukup) )   ...
%                     ) ) ;
                
%                 ev_meet_vacant=lam1*nfirm0/nfirm*max([ inner6, aucv_zeros(V2hkuplup - V1h(lup), V1hkup - V0 ) ]);
                
%                 I9=(V1h(lup) - V0 -(V2hkuplup - V1hkup)>zero_tol );
%                 ev_co_meets_vacant=lam1*nfirm0/nfirm*(  I9*contv(W1(kup, w), V1hkup - V0 , Ukup)...
%                     +(~I9)*inner6 );
                
%                 W2b(kup, lup, w) = bt*del*Ukup +...
%                     bt*ev_meet_vacant + ...
%                     bt*lam1*sum(ev_meet_solo) + ...
%                     bt*lam1*sum(sum( ev_meet_team )) + ...
%                     bt*(del + death)*contv(W1(kup, w), V1hkup - V0 , Ukup) + ...
%                     bt*ev_co_meets_vacant + ...
%                     bt*lam1*sum(ev_co_meets_solo) + ...
%                     bt*lam1*sum(sum( ev_co_meets_team )) + ...
%                     bt*lam0*sum(ev_firm_meets_u) + ...
%                     bt*lam1*sum(ev_firm_meets_solo) + ...
%                     bt*lam1*sum(sum( ev_firm_meets_team ))+ ...
%                     bt*(1 - 2*del - 2*death - 2*lam1 - (lam0*utot + lam1*etot)/nfirm)*inner6 ;
                
%             end
%         end
%     end
    
%     %Learning and collect wage
%     for w=1:wpts
%         for l=1:tpts
%             for k=1:tpts
                
%                 ev_w=zeros(tpts,1);
%                 for kup=kminus(k):kplus(k)
                    
%                     ev_w(kup)= sum((W2b(kup, kminus(l):kplus(l), w).*teamplus(l, kminus(l):kplus(l),k))*teamplus(k, kup,l));
                    
%                 end
%                 W2up(k, l, w) =   wgrid(w) + sum(ev_w);
%             end
%         end
%     end
    
%     %Check if v-funs are within tol
%     diff_w=max([max(max(abs(    W1-W1up))),max(max(max(max(abs(    W2-W2up)))))]);
%     display(['w iter: ' num2str(it_w) ' : vfunc diff w : ' num2str(diff_w)])
    
    
    
% end
% %final iteration stats
% display(['w iter: ' num2str(it_w) ' : vfunc diff w : ' num2str(diff_w)])

% toc

% save('wval.mat','W1up','W2up')
 
% %Distributions
% ftypetot   = [nfirm0, nfirm1, nfirm2]/nfirm ; %Firm type distribution
% ftype1     = solodist/nfirm1                ; %Solo firm pdf
% ftype2_mat = teamdist/(2*nfirm2)            ; %Dual firm pdf 
% ftype2     = reshape(ftype2_mat',1,tpts*tpts);
  
% %Needs to sum to precise 1, to e-14 precision
% ftypetot(end)=1-sum(ftypetot(1:end-1))  ;
% ftype1(end)  =1-sum(ftype1(1:end-1))    ;
% ftype2(end)  =1-sum(ftype2(1:end-1))    ;
 
% wtypetot      = [lam0*sum(udist),  lam1*sum(solodist),   lam1*sum(sum(teamdist))]; %worker count by job status in sampling distribution
% wtypetot_pdf  = wtypetot./sum(wtypetot)                                          ; %worker count by job status

% wtypetot_initial_pdf      = [sum(udist),  sum(solodist),   sum(sum(teamdist))]; %worker count by job status in cross section
% wtypetot_initial_pdf = wtypetot_initial_pdf./sum(wtypetot_initial_pdf);

% wtypeu        = lam0*udist/wtypetot(1)          ; %Unempl pdf
% wtypesolo     = lam1*solodist/wtypetot(2)       ; %Solo pdf
% wtypeteam_mat = lam1*teamdist/wtypetot(3)       ; %team pdf, not flattened yet
% wtypeteam     =reshape(wtypeteam_mat',1,tpts*tpts); %Flattened type distribution,

% %Needs to sum to precise 1, to e-14 precision
% wtypeu(end)      =1-sum(wtypeu(1:end-1))   ;
% wtypesolo(end)   =1-sum(wtypesolo(1:end-1));
% all_but_last=sum(sum(wtypeteam))-wtypeteam(end,end);
% wtypeteam(end,end)=1-all_but_last          ;

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Hiring rules 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for k=1:tpts
%     h0u_k(1,k) = (V1h(k) - V0 - U(k) > zero_tol); %Pr vacant firm hires unemployed type k
% end

% for k=1:tpts
%     for i=1:tpts
%         h1u_ik(i,k) = (V2h(k, i) - V1h(i) - U(k) > zero_tol); % Pr firm with type i hires unemployed type k
%     end
% end

% for i=1:tpts %Loop over firm type i
%     for j=1:tpts %Loop over firm type j
%         for k=1:tpts %Loop over my type k
%             choices = [V2h(k, j) + U(i) - V2h(i, j) - U(k),   V2h(k, i) + U(j) - V2h(i, j) - U(k), 0];
%             hire = (max(choices) > zero_tol);
%             replace =(choices(1)-choices(2) > zero_tol) +.5*(abs(choices(1)-choices(2))<=zero_tol); %1=replace i
%             h2u_ijk(i,j,k)=hire;
%             r_ijk_old(i,j,k)=replace; %probability (i,j) replaces i with k
%         end
%     end
% end

% r_ijk=zeros(tpts,tpts,tpts);
% for i=1:tpts
%     for j=1:tpts
%         for k=1:tpts
%             if V2h(k, j)-V2h(k, i) + U(i) -U(j)>zero_tol
%                 r_ijk(i,j,k)=1;
%             elseif abs(V2h(k, j)-V2h(k, i) + U(i) -U(j))<=zero_tol
%                 r_ijk(i,j,k)=.5;
%                 %(V2h(k, j)-V2h(k, i) + U(i) -U(j)>zero_tol)+.5*(abs(V2h(k, j)-V2h(k, i) + U(i) -U(j))<zero_tol);
%             end
%         end
%     end
% end

 
% %Who hires solo worker
% for k=1:tpts
%     h0e_k(1,k) = (V1h(k) - V0 - (V1h(k) - V0) > zero_tol);
% end


% for k=1:tpts %Your type
%     for i=1:tpts %Meet another solo
%         h1e_ik(i,k) = (V2h(k, i) - V1h(i) - (V1h(k) - V0) > zero_tol);
%     end
% end
 
% for i=1:tpts %Loop over firm type i
%     for j=1:tpts %Loop over firm type j
%         for k=1:tpts %Loop over my type k
%             choices = [V2h(k, j) + U(i) - V2h(i, j) - (V1h(k) - V0), V2h(k, i) + U(j) - V2h(i, j) - (V1h(k) - V0), 0];
%             hire = (max(choices) > zero_tol);
%             replace =(choices(1)-choices(2) > zero_tol) +.5*(abs(choices(1)-choices(2))<=zero_tol); %1=replace i, 2=replace j
%             h2e_ijk(i,j,k)=hire;
%             r_ijk_check(i,j,k)=replace; %Check this is the same
%         end
%     end
% end

% %Who hires team worker
% for k=1:tpts
%     for l=1:tpts %Loop over coworker type
%         h0ee_kl(k,l) = (V1h(k) - V0 - (V2h(k, l) - V1h(l)) > zero_tol);
%     end
% end

% for k=1:tpts %Loop over your type
%     for i=1:tpts %Who you meet
%         for l=1:tpts %Loop over coworker type
%             h1ee_ikl(i,k,l) = (V2h(k, i) - V1h(i) - (V2h(k, l) - V1h(l))> zero_tol);
%         end
%     end
% end

% for i=1:tpts %Loop over firm type i
%     for j=1:tpts %Loop over firm type j
%         for k=1:tpts %Loop over my type k
%             for l=1:tpts %Loop over coworker type
%                 choices =[V2h(k, j) + U(i) - V2h(i, j) - (V2h(k, l) - V1h(l)), V2h(k, i) + U(j) - V2h(i, j) - (V2h(k, l) - V1h(l)), 0];
%                 hire = (max(choices) > zero_tol);
%                 h2ee_ijkl(i,j,k,l)=hire;
%             end
%         end
%     end
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Dismissal probabilities 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dismissk_k   = zeros(tpts,1);
% dismissk_kl  = zeros(tpts,tpts);
% dismisskl_kl = zeros(tpts,tpts);
% keepkl_kl    = zeros(tpts,tpts);

% for k=1:tpts
%     if (U(k) + V0 -  V1(k) > zero_tol )
%         dismissk_k(k) = 1.00;
%     end
% end

% for l=1:tpts
%     for k=1:tpts
%         % check if firm should keep both workers
%         if ( max( [ U(k) + U(l) + V0 - V2(k,l),  V1(l) + U(k) - V2(k,l),   V1(k) + U(l) - V2(k,l)] ) < zero_tol )
%             keepkl_kl(k,l) = 1.00 ;
%         else
%             % if the firm does not keep both workers,  check if firm should dismiss both workers
%             if ( U(k) + U(l) + V0 - V2(k,l) ...
%                     - max( V1(l) + U(k) - V2(k,l),  V1(k) + U(l) - V2(k,l) ) > zero_tol )
%                 dismisskl_kl(k,l) = 1.00 ;% both separated
%             end
            
%             % if the firm does not keep both workers,  check if firm should dismiss k
%             if ( V1(l) + U(k) - V2(k,l) - (V1(k) + U(l) - V2(k,l)) > zero_tol)
%                 dismissk_kl(k,l) = 1.00;
%             elseif ( abs( V1(l) + U(k) - V2(k,l) - (V1(k) + U(l) - V2(k,l)) ) <= zero_tol )        % need to account for indifference in dismissal of k or l
%                 dismissk_kl(k,l) = 0.50  ;
%             end
%         end
%     end
% end
 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %CDFs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% max_age=12*200;
% for i=1:12*200
%     cdf_age(i) = 1.0  - (1.0 - death)^i;
% end
% cdf_age(max_age) = 1.0; %age distribution

% cdf_firm_state(1)=ftypetot(1); %{vacant,solo,team}
% cdf_firm_state(2)=sum(ftypetot(1:2));
% cdf_firm_state(3)=1;

% cdf_worker_initial_state(1)=wtypetot_initial_pdf(1); %conditional on meeting worker, what is their status {u,solo,team}
% cdf_worker_initial_state(2)=sum(wtypetot_initial_pdf(1:2));
% cdf_worker_initial_state(3)=1;

% prob_firm_contact_worker=sum(wtypetot); % 

% cdf_worker_state(1)=wtypetot(1)/prob_firm_contact_worker; % conditional on meeting worker, what is their type {u,solo,team}
% cdf_worker_state(2)=sum(wtypetot(1:2))/prob_firm_contact_worker;
% cdf_worker_state(3)=1;

% Pi_k=typebirth;
% cdf_newborn=measure_to_cdf(Pi_k);
% gamma=bpw;
% delta=del;
% sigma=death;
% lambda1=lam1;
% lambda0=lam0;
 
% p1_i= solodist ; %solo distribution of types
% p2_ij= teamdist/2; %team distribution of types
% p0=1-sum(p1_i)-sum(sum(p2_ij));

% cdf_solo=measure_to_cdf(solodist); % 
% cdf_team=measure2d_to_cdf(p2_ij);
% cdf_worker_unempl=measure_to_cdf(udist);


% wage_w=wgrid; %Wage grid
% W_wk=permute(W1,[2,1]); %Worker values -- wage is first dimension
% W_wkl=permute(W2,[3, 1 , 2]); %Worker team value -- wage is first dimension, then you, then co
% U_k=U; 
% Vhat_k=V1h;
% Vhat_kl=V2h;


% nh = tpts;
% gu          = ulose;  
% g0          = sologain;  
% g1plus      = team_learn_param_up; 
% g1minus     = team_learn_param_up_sym;  

% pr.gu_down_k = repmat(gu,nh,1) ;
% pr.gu_down_k(1) = 0.00;
% pr.gu_stay_k = 1.00 - pr.gu_down_k;

% pr.g_up_k = repmat(g0,nh,1);
% pr.g_up_k(nh) = 0.00;
% pr.g_stay_k = 1.00 - pr.g_up_k;

% pr.g_up_kl(nh,:) = 0.00;
% tmp = 1.00 / ( type(nh) - type(1) ) ;  
% for k=1:(nh-1)
%     for l=1:nh
%         pr.g_up_kl(k,l) = g0 + tmp * ( g1plus * max(type(l)-type(k), 0.00 ) ...
%             + g1minus * max(type(k)-type(l), 0.00 ) ) ;
%     end
% end
% pr.g_stay_kl = 1.00 - pr.g_up_kl;
 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Compute marginals
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Vmarginal_k = Vhat_k - V0;

% for k=1:nh
%     for i=1:nh
%         Vmarginal_ki(k,i) = Vhat_kl(k,i) - Vhat_k(i);
%     end
% end

% for k=1:nh
%     for i=1:nh
%         for j=1:nh
%             Vmarginal_kij(k,i,j) = max( Vhat_kl(k,j) - Vhat_kl(i,j) + U_k(i), ...
%                 Vhat_kl(i,k) - Vhat_kl(i,j) + U_k(j) );
%         end
%     end
% end
 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Calculate wage bounds for all matches
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k=1:tpts
%     wage_target_U_solo(k)=InterpolateWage( W_wk(:,k), U_k(k), wage_w );
% end

% for k=1:tpts
%     wage_target_MP_solo(k)=InterpolateWage( W_wk(:,k), Vmarginal_k(k), wage_w );
% end

% for l=1:tpts
%     for k=1:tpts
%         wage_target_U_team(k,l)=InterpolateWage( W_wkl(:,k,l), U_k(k), wage_w );
%     end
% end

% for l=1:tpts
%     for k=1:tpts
%         wage_target_MP_team(k,l)=InterpolateWage( W_wkl(:,k,l), Vmarginal_ki(k,l), wage_w );
%     end
% end
 
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%
% %Monte Carlo
% %%%%%%%%%%%%%%%%%%%%%%%%%
% n_months=240;
% n_workers=200000;
% n_years=n_months/12;

% STAT_store=[];
% Nrep=5; %Repeat simulations
% for jjj=1:Nrep
    
%     jjj
%     rng(jjj,'twister'); %Set the seed
    
    
%     r.state = rand(n_months,n_workers);  
%     r.age =rand(n_months,n_workers);  
%     r.type = rand(n_months,n_workers);  
%     r.learn =rand(n_months,n_workers);  
%     r.event =rand(n_months,n_workers);  
%     r.tie_break = rand(n_months,n_workers);  
%     r.co_learn = rand(n_months,n_workers);  
%     r.death =rand(n_months,n_workers); 
%     r.age  =rand(n_months,n_workers);  
     
%     %%%%%%%%%%%%%%%
%     %Sim params
%     %%%%%%%%%%%%%%%
%     burn_years=10;
%     t_burn=12*burn_years;
%     stable_years = n_years - burn_years;
%     nquantiles=10;
%     t = 1; %
%     n_firms = n_workers; % 
    
%     %Initialize storage of firm variables
%     firm_status_mi=zeros(n_months,n_firms);
%     firm_t1_mi =zeros(n_months,n_firms);
%     firm_t2_mi =zeros(n_months,n_firms);
%     firm_w1_mi =zeros(n_months,n_firms);
%     firm_w2_mi =zeros(n_months,n_firms);
    
    
%     for i=1:n_firms 
%         firm_status_mi(t,i) = draw_CDF_1d( cdf_firm_state, r.state(t,i));
%         firm_status_mi(t,i) = firm_status_mi(t,i) - 1;  
%         if ( firm_status_mi(t,i) == 0 )% Firm is vacant
            
 
%             firm_t1_mi(t,i) = 0;
%             firm_t2_mi(t,i) = 0;
%             firm_w1_mi(t,i) = 0.0;
%             firm_w2_mi(t,i) = 0.0;
%         elseif ( firm_status_mi(t,i) == 1 )% Firm has one worker
       
%             if (r.tie_break(t,i) < 0.5)
%                 firm_t1_mi(t,i)=draw_CDF_1d( cdf_solo, r.type(t,i) );
                
         
%                 firm_t2_mi(t,i) = 0 ;
%                 target_value = (1.0-gamma) * U_k( firm_t1_mi(t,i) ) + gamma * Vmarginal_k( firm_t1_mi(t,i) );
%                 firm_w1_mi(t,i)=InterpolateWage( W_wk(:,firm_t1_mi(t,i)), target_value, wage_w );
%                 firm_w2_mi(t,i) = 0.0;
%             else
%                 firm_t2_mi(t,i)=draw_CDF_1d( cdf_solo, r.type(t,i)  );
                
             
%                 firm_t1_mi(t,i) = 0;
%                 target_value = (1.0-gamma) * U_k( firm_t2_mi(t,i) ) + gamma * Vmarginal_k( firm_t2_mi(t,i) );
%                 firm_w2_mi(t,i)=InterpolateWage( W_wk(:,firm_t2_mi(t,i)), target_value, wage_w );
%                 firm_w1_mi(t,i) = 0.0;
%             end
%         else  % Firm has a team
%             [firm_t1_mi(t,i), firm_t2_mi(t,i) ]=draw_CDF_2d( cdf_team, nh, r.type(t,i) );
            
      
%             target_value = (1.0-gamma) * U_k( firm_t1_mi(t,i) ) ...
%                 + gamma * Vmarginal_ki( firm_t1_mi(t,i), firm_t2_mi(t,i) );
%             firm_w1_mi(t,i) =InterpolateWage( W_wkl(:,firm_t1_mi(t,i), firm_t2_mi(t,i)), target_value, wage_w);
%             target_value = (1.0-gamma) * U_k( firm_t2_mi(t,i) ) + gamma * Vmarginal_ki( firm_t2_mi(t,i), firm_t1_mi(t,i) );
%             firm_w2_mi(t,i)=InterpolateWage( W_wkl(:,firm_t2_mi(t,i), firm_t1_mi(t,i)), target_value, wage_w);
%         end
%     end
      
%     for i=1:n_firms
        
%         for t=2:n_months  
     
%             % Vacant Firm
%             if ( firm_status_mi(t-1,i) == 0 )% vacant firm 
%                 if (r.event(t,i) <= prob_firm_contact_worker )% firm contacts a worker_type
%                     if ( r.state(t,i) <= cdf_worker_state(1) )% meet an unemployed worker
                
%                         ta=draw_CDF_1d( cdf_worker_unempl, r.type(t,i));
%                         if ( h0u_k(ta) > 0.0 )% firm hires unemployed worker
%                             if (r.tie_break(t,i) < 0.5) % flip a coin for which slot to use
 
%                                 firm_status_mi(t,i) = 1;
%                                 firm_t1_mi(t,i) = ta;
%                                 target_value = (1.0 - gamma) * U_k(ta) + gamma * Vmarginal_k(ta);
%                                 firm_w1_mi(t,i)=InterpolateWage( W_wk(:,ta), target_value, wage_w );
%                                 firm_t2_mi(t,i) = 0;
%                                 firm_w2_mi(t,i) = 0.0;
%                             else
                     
%                                 firm_status_mi(t,i) = 1;
%                                 firm_t2_mi(t,i) = ta;
%                                 target_value = (1.0 - gamma) * U_k(ta) + gamma * Vmarginal_k(ta);
%                                 firm_w2_mi(t,i)=InterpolateWage( W_wk(:,ta), target_value, wage_w);
%                                 firm_t1_mi(t,i) = 0;
%                                 firm_w1_mi(t,i) = 0.0;
%                             end
%                         else  % firm does not hire unemployed worker
          
%                             firm_status_mi(t,i) = 0;
%                             firm_t1_mi(t,i) = 0;
%                             firm_t2_mi(t,i) = 0;
%                             firm_w1_mi(t,i) = 0.0;
%                             firm_w2_mi(t,i) = 0.0;
%                         end
%                     elseif ( r.state(t,i) > cdf_worker_state(1) ...
%                             && r.state(t,i) <= cdf_worker_state(2) ) % meet solo worker
 
%                         firm_status_mi(t,i) = 0;
%                         firm_t1_mi(t,i) = 0;
%                         firm_t2_mi(t,i) = 0;
%                         firm_w1_mi(t,i) = 0.0;
%                         firm_w2_mi(t,i) = 0.0;
%                     else  % meet a worker on a team 
%                         [ta, tb]=draw_CDF_2d( cdf_team, nh, r.type(t,i));
%                         if ( h0ee_kl(ta,tb) > 0.0 )% firm hires worker
%                             if (r.tie_break(t,i) < 0.5)
 
%                                 firm_status_mi(t,i) = 1;
%                                 firm_t1_mi(t,i) = ta;
%                                 target_value = (1.0 - gamma) * Vmarginal_ki(ta,tb) + gamma * Vmarginal_k(ta);
%                                 firm_w1_mi(t,i) = InterpolateWage( W_wk(:,ta), target_value, wage_w);
%                                 firm_t2_mi(t,i) = 0;
%                                 firm_w2_mi(t,i) = 0.0;
%                             else
 
%                                 firm_status_mi(t,i) = 1;
%                                 firm_t2_mi(t,i) = ta;
%                                 target_value = (1.0 - gamma) * Vmarginal_ki(ta,tb) + gamma * Vmarginal_k(ta);
%                                 firm_w2_mi(t,i) = InterpolateWage( W_wk(:,ta), target_value, wage_w );
%                                 firm_t1_mi(t,i) = 0;
%                                 firm_w1_mi(t,i) = 0.0;
                                
%                             end
%                         else  % firm does not hire  worker
 
%                             firm_status_mi(t,i) = 0;
%                             firm_t1_mi(t,i) = 0;
%                             firm_t2_mi(t,i) = 0;
%                             firm_w1_mi(t,i) = 0.0;
%                             firm_w2_mi(t,i) = 0.0;
%                         end
%                     end
%                 else  % firm does not contact any worker
 
%                     firm_status_mi(t,i) = 0;
%                     firm_t1_mi(t,i) = 0;
%                     firm_t2_mi(t,i) = 0;
%                     firm_w1_mi(t,i) = 0.0;
%                     firm_w2_mi(t,i) = 0.0;
%                 end
%                 % Firm with one worker
%             elseif ( firm_status_mi(t-1,i) == 1 )% firm with one worker_type
 
%                 % check if worker learns
%                 if ( firm_t1_mi(t-1,i) > 0)% worker in slot 1
%                     if ( r.learn(t,i) >= pr.g_stay_k(firm_t1_mi(t-1,i)) )% gain human capital
%                         firm_t1_mi(t,i) = firm_t1_mi(t-1,i) + 1;
%                         firm_t2_mi(t,i) = 0;
  
%                     else  % human capital stays the same
%                         firm_t1_mi(t,i) = firm_t1_mi(t-1,i);
%                         firm_t2_mi(t,i) = 0;
       
%                     end
%                 else  % worker in slot 2
%                     if ( r.learn(t,i) >= pr.g_stay_k(firm_t2_mi(t-1,i)) )% gain human capital
%                         firm_t1_mi(t,i) = 0;
%                         firm_t2_mi(t,i) = firm_t2_mi(t-1,i) + 1;
               
%                     else  % human capital stays the same
%                         firm_t1_mi(t,i) = 0;
%                         firm_t2_mi(t,i) = firm_t2_mi(t-1,i);
          
%                     end
%                 end
                
%                 if ( r.event(t,i) <= delta + sigma ) % worker leaves firm
          
%                     firm_status_mi(t,i) = 0;
%                     firm_t1_mi(t,i) = 0;
%                     firm_t2_mi(t,i) = 0;
%                     firm_w1_mi(t,i) = 0.0;
%                     firm_w2_mi(t,i) = 0.0;
%                 elseif ( r.event(t,i) > delta + sigma & ...
%                         r.event(t,i) <= delta + sigma + prob_firm_contact_worker )% firm contacts a worker
%                     if ( r.state(t,i) <= cdf_worker_state(1) )% meet an unemployed worker
        
%                         ta=draw_CDF_1d( cdf_worker_unempl, r.type(t,i) );
%                         if ( firm_t1_mi(t,i) > 0 )% there is a worker in slot 1, fill slot 2
%                             if ( h1u_ik(firm_t1_mi(t,i),ta) > 0.0 )% firm hires unemployed worker
 
%                                 firm_status_mi(t,i) = 2;
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                                 firm_t2_mi(t,i) = ta;
%                                 target_value = (1.0 - gamma) * U_k(ta) + gamma * Vmarginal_ki(ta,firm_t1_mi(t,i));
%                                 firm_w2_mi(t,i)=InterpolateWage( W_wkl(:,ta,firm_t1_mi(t,i)), ...
%                                     target_value, wage_w);
%                             else  % firm does not hire unemployed worker
                      
%                                 firm_status_mi(t,i) = 1;
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                                 firm_w2_mi(t,i) = 0.0;
%                             end
                            
%                         else  % There is a worker in slot 2, fill slot 1
%                             if ( h1u_ik(firm_t2_mi(t,i),ta) > 0.0 )% firm hires unemployed worker
                 
%                                 firm_status_mi(t,i) = 2;
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                                 firm_t1_mi(t,i) = ta;
%                                 target_value = (1.0 - gamma) * U_k(ta) + gamma * Vmarginal_ki(ta,firm_t2_mi(t,i));
%                                 firm_w1_mi(t,i) =InterpolateWage( W_wkl(:,ta,firm_t2_mi(t,i)), ...
%                                     target_value, wage_w );
%                             else  % firm does not hire unemployed worker
                     
%                                 firm_status_mi(t,i) = 1;
%                                 firm_w1_mi(t,i) = 0.0;
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                             end
%                         end
%                     elseif ( r.state(t,i) > cdf_worker_state(1) ...
%                             & r.state(t,i) <= cdf_worker_state(2) ) % meet solo worker
            
%                         ta=draw_CDF_1d( cdf_solo, r.type(t,i) );
%                         if ( firm_t1_mi(t,i) > 0 )% there is a worker in slot 1, fill slot 2
%                             if (h1e_ik(firm_t1_mi(t,i),ta) > 0.0 )% firm hires worker
                
%                                 firm_status_mi(t,i) = 2;
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                                 firm_t2_mi(t,i) = ta;
%                                 target_value = (1.0 - gamma) * Vmarginal_k(ta) ...
%                                     + gamma * Vmarginal_ki(ta,firm_t1_mi(t,i));
%                                 firm_w2_mi(t,i)=InterpolateWage( W_wkl(:,ta,firm_t1_mi(t,i)), target_value, ...
%                                     wage_w);
%                             else  % firm does not hire solo  worker
                        
%                                 firm_status_mi(t,i) = 1;
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                                 firm_w2_mi(t,i) = 0.0;
%                             end
%                         else  % there is a worker in slot 2, fill slot 1
%                             if (h1e_ik(firm_t2_mi(t,i),ta) > 0.0 )% firm hires worker
                               
%                                 firm_status_mi(t,i) = 2;
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                                 firm_t1_mi(t,i) = ta;
%                                 target_value = (1.0 - gamma) * Vmarginal_k(ta) ...
%                                     + gamma * Vmarginal_ki(ta,firm_t2_mi(t,i));
%                                 firm_w1_mi(t,i)=InterpolateWage( W_wkl(:,ta,firm_t2_mi(t,i)), ...
%                                     target_value, wage_w);
%                             else  % firm does not hire solo  worker
                          
%                                 firm_status_mi(t,i) = 1;
%                                 firm_w1_mi(t,i) = 0.0;
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                             end
%                         end
%                     else  % meet a worker on a team
               
%                         [ta, tb]=draw_CDF_2d( cdf_team, nh, r.type(t,i));
%                         if ( firm_t1_mi(t,i) > 0 )% there is a worker in slot 1, fill slot 2
%                             if ( h1ee_ikl(firm_t1_mi(t,i),ta,tb) > 0.0 )% firm hires worker
                 
%                                 firm_status_mi(t,i) = 2;
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                                 firm_t2_mi(t,i) = ta;
%                                 target_value = (1.0 - gamma) * Vmarginal_ki(ta,tb) ...
%                                     + gamma * Vmarginal_ki(ta,firm_t1_mi(t,i));
%                                 firm_w2_mi(t,i)=InterpolateWage( W_wkl(:,ta,firm_t1_mi(t,i)), ...
%                                     target_value, wage_w );
%                             else  % firm does not hire  worker
                        
%                                 firm_status_mi(t,i) = 1;
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                                 firm_w2_mi(t,i) = 0.0;
%                             end
%                         else  % there is a worker in slot 2, fill slot 1
%                             if ( h1ee_ikl(firm_t2_mi(t,i),ta,tb) > 0.0 )% firm hires worker
                 
%                                 firm_status_mi(t,i) = 2;
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                                 firm_t1_mi(t,i) = ta;
%                                 target_value = (1.0 - gamma) * Vmarginal_ki(ta,tb) ...
%                                     + gamma * Vmarginal_ki(ta,firm_t2_mi(t,i));
%                                 firm_w1_mi(t,i)=InterpolateWage( W_wkl(:,ta,firm_t2_mi(t,i)), ...
%                                     target_value, wage_w);
%                             else
                      
%                                 firm_status_mi(t,i) = 1;
%                                 firm_w1_mi(t,i) = 0.0;
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                             end
%                         end
%                     end
%                 elseif ( r.event(t,i) > delta + sigma + prob_firm_contact_worker & ...
%                         r.event(t,i) <= delta + sigma + prob_firm_contact_worker + lambda1 )% Worker meets a firm
      
%                     if ( r.state(t,i) <= cdf_firm_state(1) )% meet vacant firm 
     
%                         if ( firm_t1_mi(t,i) > 0 )% there is a worker in slot 1
%                             firm_status_mi(t,i) = 1;
%                             firm_w1_mi(t,i) = wage_target_MP_solo( firm_t1_mi(t,i) );
%                             firm_w2_mi(t,i) = 0.0;
%                         else  % there is a worker in slot 2
%                             firm_status_mi(t,i) = 1;
%                             firm_w1_mi(t,i) = 0.0;
%                             firm_w2_mi(t,i) = wage_target_MP_solo( firm_t2_mi(t,i) );
%                         end
%                     elseif ( r.state(t,i) > cdf_firm_state(1) & ...
%                             r.state(t,i) <= cdf_firm_state(2) )% meet a solo firm
 
%                         ta=draw_CDF_1d( cdf_solo, r.type(t,i));
%                         if ( firm_t1_mi(t,i) > 0 )% there is a worker in slot 1
%                             if ( h1e_ik(ta, firm_t1_mi(t,i)) > 0.0 )% worker leaves
                
%                                 firm_status_mi(t,i) = 0;
%                                 firm_t1_mi(t,i) = 0;
%                                 firm_w1_mi(t,i) = 0.0;
%                                 firm_t2_mi(t,i) = 0;
%                                 firm_w2_mi(t,i) = 0.0;
%                             else  % worker stays
                           
%                                 firm_status_mi(t,i) = 1;
%                                 % check if worker bids up her wage
%                                 target_value = Vmarginal_ki(firm_t1_mi(t,i),ta) ;    
%                                 current_value=InterpolateWage( wage_w, firm_w1_mi(t-1,i), ...
%                                     W_wk(:,firm_t1_mi(t,i)));
%                                 if ( target_value > current_value )% worker renegotiates wage with current firm 
%                                     firm_w1_mi(t,i)=InterpolateWage( W_wk(:,firm_t1_mi(t,i)), ...
%                                         target_value, wage_w);
%                                 else
%                                     firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                                 end
%                                 firm_t2_mi(t,i) = 0;
%                                 firm_w2_mi(t,i) = 0.0;
%                             end
%                         else  % Worker is in slot 2
%                             if ( h1e_ik(ta, firm_t2_mi(t,i)) > 0.0 )% worker leaves 
%                                 firm_status_mi(t,i) = 0;
%                                 firm_t1_mi(t,i) = 0;
%                                 firm_w1_mi(t,i) = 0.0;
%                                 firm_t2_mi(t,i) = 0;
%                                 firm_w2_mi(t,i) = 0.0;
%                             else  % worker stays 
%                                 firm_status_mi(t,i) = 1;
%                                 % check if worker bids up her wage
%                                 target_value = Vmarginal_ki(firm_t2_mi(t,i),ta) ;  
%                                 current_value =InterpolateWage( wage_w, firm_w2_mi(t-1,i), ...
%                                     W_wk(:,firm_t2_mi(t,i)));
%                                 if ( target_value > current_value )% worker renegotiates wage with current firm 
%                                     firm_w2_mi(t,i)=InterpolateWage( W_wk(:,firm_t2_mi(t,i)), ...
%                                         target_value, wage_w );
%                                 else
%                                     firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                                 end
%                                 firm_t1_mi(t,i) = 0;
%                                 firm_w1_mi(t,i) = 0.0;
%                             end
%                         end
%                     else  % worker meets team
%                         [ta, tb]=draw_CDF_2d( cdf_team, nh, r.type(t,i)); 
%                         if ( firm_t1_mi(t,i) > 0 )% there is a worker in slot 1
%                             if ( h2e_ijk(ta, tb, firm_t1_mi(t,i)) > 0.0 )% worker leaves 
%                                 firm_status_mi(t,i) = 0;
%                                 firm_t1_mi(t,i) = 0;
%                                 firm_w1_mi(t,i) = 0.0;
%                                 firm_t2_mi(t,i) = 0;
%                                 firm_w2_mi(t,i) = 0.0;
%                             else  % worker stays 
%                                 firm_status_mi(t,i) = 1;
%                                 % check if worker bids up her wage
%                                 target_value = Vmarginal_kij(firm_t1_mi(t,i),ta,tb) ; 
%                                 current_value =InterpolateWage( wage_w, firm_w1_mi(t-1,i), ...
%                                     W_wk(:,firm_t1_mi(t,i)));
%                                 if ( target_value > current_value )% worker renegotiates wage with current firm 
%                                     firm_w1_mi(t,i)=InterpolateWage( W_wk(:,firm_t1_mi(t,i)), ...
%                                         target_value, wage_w );
%                                 else
%                                     firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                                 end
%                                 firm_t2_mi(t,i) = 0;
%                                 firm_w2_mi(t,i) = 0.0;
%                             end
%                         else  % the worker is in slot 2
%                             if ( h2e_ijk(ta, tb, firm_t2_mi(t,i)) > 0.0 )% worker leaves 
%                                 firm_status_mi(t,i) = 0;
%                                 firm_t1_mi(t,i) = 0;
%                                 firm_w1_mi(t,i) = 0.0;
%                                 firm_t2_mi(t,i) = 0;
%                                 firm_w2_mi(t,i) = 0.0;
%                             else  % the worker stays 
%                                 firm_status_mi(t,i) = 1;
%                                 % check if worker bids up her wage
%                                 target_value = Vmarginal_kij(firm_t2_mi(t,i),ta,tb); 
%                                 current_value=InterpolateWage( wage_w, firm_w2_mi(t-1,i), ...
%                                     W_wk(:,firm_t2_mi(t,i)));
%                                 if ( target_value > current_value )% worker renegotiates wage with current firm
                                 
%                                     firm_w2_mi(t,i)=InterpolateWage( W_wk(:,firm_t2_mi(t,i)), ...
%                                         target_value, wage_w );
%                                 else
%                                     firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                                 end
%                                 firm_t1_mi(t,i) = 0;
%                                 firm_w1_mi(t,i) = 0.0;
%                             end
%                         end
%                     end
%                 else  % no separations or contacts
                   
%                     firm_status_mi(t,i) = 1;
%                     firm_w1_mi(t,i) = firm_w1_mi(t-1,i); % one of these will be zero
%                     firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                 end
                
%             else  % firm has a team
           
%                 % check if worker 1 learns
%                 if ( r.learn(t,i) >= pr.g_stay_kl(firm_t1_mi(t-1,i), firm_t2_mi(t-1,i)) )% gain human capital
%                     firm_t1_mi(t,i) = firm_t1_mi(t-1,i) + 1;
                  
%                 else  % human capital stays the same
%                     firm_t1_mi(t,i) = firm_t1_mi(t-1,i);
                    
%                 end
%                 % check if worker 2 learns
%                 if ( r.co_learn(t,i) >= pr.g_stay_kl(firm_t2_mi(t-1,i), firm_t1_mi(t-1,i)) )% gain human capital
%                     firm_t2_mi(t,i) = firm_t2_mi(t-1,i) + 1;
                   
%                 else  % human capital stays the same
%                     firm_t2_mi(t,i) = firm_t2_mi(t-1,i);
                 
%                 end
                
%                 if ( r.event(t,i) <= delta + sigma )% worker 1 exits
      
%                     firm_status_mi(t,i) = 1;
%                     firm_t1_mi(t,i) = 0;
%                     firm_w1_mi(t,i) = 0.0;
%                     firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                 elseif (r.event(t,i) > delta + sigma & ...
%                         r.event(t,i) <= delta + sigma + delta + sigma )% worker 2 exits
         
%                     firm_status_mi(t,i) = 1;
%                     firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                     firm_t2_mi(t,i) = 0;
%                     firm_w2_mi(t,i) = 0.0;
%                 elseif (r.event(t,i) >  delta + sigma + delta + sigma & ...
%                         r.event(t,i) <= delta + sigma + delta + sigma + prob_firm_contact_worker ) % firm contacts a new worker
          
%                     if ( r.state(t,i) <= cdf_worker_state(1) )% meet an unemployed worker
%                         ta=draw_CDF_1d( cdf_worker_unempl, r.type(t,i) );
%                         if ( h2u_ijk(firm_t1_mi(t,i),firm_t2_mi(t,i),ta) > 0.0 )% firm hires unemployed worker
                
%                             firm_status_mi(t,i) = 2;
%                             % check which worker is replaced
%                             if ( r.tie_break(t,i) < r_ijk(firm_t1_mi(t,i),firm_t2_mi(t,i),ta) )% replace worker 1
                      
%                                 target_value = (1.0 - gamma) * U_k(ta) ...
%                                     + gamma * Vmarginal_kij(ta,firm_t1_mi(t,i),firm_t2_mi(t,i));
%                                 firm_w1_mi(t,i)=InterpolateWage( W_wkl(:,ta,firm_t2_mi(t,i)), ...
%                                     target_value, wage_w );
%                                 firm_t1_mi(t,i) = ta;
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                             else  % replace worker 2
                           
%                                 target_value = (1.0 - gamma) * U_k(ta) ...
%                                     + gamma * Vmarginal_kij(ta,firm_t1_mi(t,i),firm_t2_mi(t,i));
%                                 firm_w2_mi(t,i)=InterpolateWage( W_wkl(:,ta,firm_t1_mi(t,i)), ...
%                                     target_value, wage_w );
%                                 firm_t2_mi(t,i) = ta;
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             end
%                         else  % firm does not hire unemployed worker
                        
%                             firm_status_mi(t,i) = 2;
%                             firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                         end
%                     elseif ( r.state(t,i) > cdf_worker_state(1) ...
%                             & r.state(t,i) <= cdf_worker_state(2) ) % meet solo worker
%                         ta=draw_CDF_1d( cdf_solo, r.type(t,i));
%                         if (h2e_ijk(firm_t1_mi(t,i),firm_t2_mi(t,i),ta) > 0.0 )% firm hires worker
                   
%                             firm_status_mi(t,i) = 2;
%                             % check which worker is replaced
%                             if ( r.tie_break(t,i) < r_ijk(firm_t1_mi(t,i),firm_t2_mi(t,i),ta) )% replace worker 1
                            
%                                 target_value = (1.0 - gamma) * Vmarginal_k(ta) ...
%                                     + gamma * Vmarginal_kij(ta,firm_t1_mi(t,i),firm_t2_mi(t,i));
%                                 firm_w1_mi(t,i)=InterpolateWage(W_wkl(:,ta,firm_t2_mi(t,i)), ...
%                                     target_value, wage_w) ;
%                                 firm_t1_mi(t,i) = ta;
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                             else  % replace worker 2
                    
%                                 target_value = (1.0 - gamma) * Vmarginal_k(ta) ...
%                                     + gamma * Vmarginal_kij(ta,firm_t1_mi(t,i),firm_t2_mi(t,i));
%                                 firm_w2_mi(t,i)=InterpolateWage( W_wkl(:,ta,firm_t1_mi(t,i)), ...
%                                     target_value, wage_w );
%                                 firm_t2_mi(t,i) = ta;
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             end
%                         else  % firm does not hire solo  worker
                          
%                             firm_status_mi(t,i) = 2;
%                             firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                         end
%                     else  % meet a worker on a team
%                         [ta, tb]=draw_CDF_2d( cdf_team, nh, r.type(t,i) );
%                         if ( h2ee_ijkl(firm_t1_mi(t,i),firm_t2_mi(t,i),ta,tb) > 0.0 )% firm hires worker
                     
%                             firm_status_mi(t,i) = 2;
%                             if ( r.tie_break(t,i) < r_ijk(firm_t1_mi(t,i),firm_t2_mi(t,i),ta) )% replace worker 1
                          
%                                 target_value = (1.0 - gamma) * Vmarginal_ki(ta,tb) ...
%                                     + gamma * Vmarginal_kij(ta,firm_t1_mi(t,i),firm_t2_mi(t,i));
%                                 firm_w1_mi(t,i)=InterpolateWage( W_wkl(:,ta,firm_t2_mi(t,i)), ...
%                                     target_value, wage_w );
%                                 firm_t1_mi(t,i) = ta;
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                             else  % replace worker 2
                              
%                                 target_value = (1.0 - gamma) * Vmarginal_ki(ta,tb) ...
%                                     + gamma * Vmarginal_kij(ta,firm_t1_mi(t,i),firm_t2_mi(t,i));
%                                 firm_w2_mi(t,i)=InterpolateWage( W_wkl(:,ta,firm_t1_mi(t,i)), ...
%                                     target_value, wage_w );
%                                 firm_t2_mi(t,i) = ta;
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             end
%                         else  % firm does not hire team  worker
                            
%                             firm_status_mi(t,i) = 2;
%                             firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                         end
%                     end
%                 elseif ( r.event(t,i) > delta + sigma + delta + sigma + prob_firm_contact_worker & ...
%                         r.event(t,i) <= delta + sigma + delta + sigma ...
%                         + prob_firm_contact_worker + lambda1 )% worker 1 contacts a firm
                    
%                     if ( r.state(t,i) <= cdf_firm_state(1) )% meet vacant firm
                     
%                         if ( h0ee_kl(firm_t1_mi(t,i),firm_t2_mi(t,i)) > 0.0 )% worker 1 leaves
                         
%                             firm_status_mi(t,i) = 1;
%                             firm_t1_mi(t,i) = 0;
%                             firm_w1_mi(t,i) = 0.0;
%                             firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                         else  % worker stays
                          
%                             firm_status_mi(t,i) = 2;
%                             target_value = Vmarginal_k(firm_t1_mi(t,i))  ;
%                             current_value = InterpolateWage( wage_w, firm_w1_mi(t-1,i), ...
%                                 W_wkl(:,firm_t1_mi(t,i),firm_t2_mi(t,i)));
%                             if ( target_value > current_value )% worker renegotiates wage with current firm
                               
%                                 firm_w1_mi(t,i)=InterpolateWage( W_wkl(:,firm_t1_mi(t,i),firm_t2_mi(t,i)), ...
%                                     target_value, wage_w );
%                             else
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             end
%                             firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                         end
%                     elseif ( r.state(t,i) > cdf_firm_state(1) & ...
%                             r.state(t,i) <= cdf_firm_state(2) )% meet a solo firm
                         
%                         ta=draw_CDF_1d( cdf_solo, r.type(t,i) );
%                         if ( h1ee_ikl(ta, firm_t1_mi(t,i), firm_t2_mi(t,i)) > 0.0 )% worker leaves
                             
%                             firm_status_mi(t,i) = 1;
%                             firm_t1_mi(t,i) = 0;
%                             firm_w1_mi(t,i) = 0.0;
%                             firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                         else  % worker stays
                             
%                             firm_status_mi(t,i) = 2;
%                             % check if worker bids up her wage
%                             target_value = Vmarginal_ki(firm_t1_mi(t,i),ta); 
%                             current_value=InterpolateWage( wage_w, firm_w1_mi(t-1,i), ...
%                                 W_wkl(:,firm_t1_mi(t,i),firm_t2_mi(t,i)) );
%                             if ( target_value > current_value )% worker renegotiates wage with current firm
                                
%                                 firm_w1_mi(t,i)=InterpolateWage( W_wkl(:,firm_t1_mi(t,i),firm_t2_mi(t,i)), ...
%                                     target_value, wage_w );
%                             else
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             end
%                             firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                         end
%                     else  % worker meets team
                 
%                         [ta, tb]=draw_CDF_2d( cdf_team, nh, r.type(t,i) );
%                         if ( h2ee_ijkl(ta, tb, firm_t1_mi(t,i), firm_t2_mi(t,i)) > 0.0 )% worker leaves
                            
%                             firm_status_mi(t,i) = 1;
%                             firm_t1_mi(t,i) = 0;
%                             firm_w1_mi(t,i) = 0.0;
%                             firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                         else  % worker stays
                         
%                             firm_status_mi(t,i) = 2;
%                             % check if worker bids up her wage
%                             target_value = Vmarginal_kij(firm_t1_mi(t,i),ta,tb);  
%                             current_value=InterpolateWage( wage_w, firm_w1_mi(t-1,i), ...
%                                 W_wkl(:,firm_t1_mi(t,i),firm_t2_mi(t,i)) );
%                             if ( target_value > current_value )% worker renegotiates wage with current firm
                             
%                                 firm_w1_mi(t,i)=InterpolateWage( W_wkl(:,firm_t1_mi(t,i),firm_t2_mi(t,i)), ...
%                                     target_value, wage_w );
%                             else
%                                 firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             end
%                             firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                         end
%                     end
                    
%                 elseif ( r.event(t,i) > delta + sigma + delta + sigma ...
%                         + prob_firm_contact_worker + lambda1 ...
%                         & r.event(t,i) <= delta + sigma + delta + sigma ...
%                         + prob_firm_contact_worker + lambda1 + lambda1 )% worker 2 contacts a firm
                    
%                     if ( r.state(t,i) <= cdf_firm_state(1) )% meet vacant firm
            
%                         if ( h0ee_kl(firm_t2_mi(t,i),firm_t1_mi(t,i)) > 0.0 )% worker 2 leaves
                  
%                             firm_status_mi(t,i) = 1;
%                             firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             firm_t2_mi(t,i) = 0;
%                             firm_w2_mi(t,i) = 0.0;
%                         else  % worker stays 
%                             firm_status_mi(t,i) = 2;
%                             target_value = Vmarginal_k(firm_t2_mi(t,i)) ; 
%                             current_value=InterpolateWage( wage_w, firm_w2_mi(t-1,i), ...
%                                 W_wkl(:,firm_t2_mi(t,i),firm_t1_mi(t,i)) );
%                             if ( target_value > current_value )% worker renegotiates wage with current firm
                        
%                                 firm_w2_mi(t,i)=InterpolateWage( W_wkl(:,firm_t2_mi(t,i),firm_t1_mi(t,i)), ...
%                                     target_value, wage_w );
%                             else
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                             end
%                             firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                         end
%                     elseif ( r.state(t,i) > cdf_firm_state(1) & ...
%                             r.state(t,i) <= cdf_firm_state(2) )% meet a solo firm
                   
%                         ta=draw_CDF_1d( cdf_solo, r.type(t,i) );
%                         if ( h1ee_ikl(ta, firm_t2_mi(t,i), firm_t1_mi(t,i)) > 0.0 )% worker 2 leaves
                     
%                             firm_status_mi(t,i) = 1;
%                             firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             firm_t2_mi(t,i) = 0;
%                             firm_w2_mi(t,i) = 0.0;
%                         else  % worker stays
                      
%                             firm_status_mi(t,i) = 2;
%                             % check if worker bids up her wage
%                             target_value = Vmarginal_ki(firm_t2_mi(t,i),ta); 
%                             current_value =InterpolateWage( wage_w, firm_w2_mi(t-1,i), ...
%                                 W_wkl(:,firm_t2_mi(t,i),firm_t1_mi(t,i)));
%                             if ( target_value > current_value )% worker renegotiates wage with current firm 
%                                 firm_w2_mi(t,i)=InterpolateWage( W_wkl(:,firm_t2_mi(t,i),firm_t1_mi(t,i)), ...
%                                     target_value, wage_w );
%                             else
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                             end
%                             firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                         end
                        
                        
%                     else  % worker meets team
          
%                         [ta, tb ]=draw_CDF_2d( cdf_team, nh, r.type(t,i));
%                         if ( h2ee_ijkl(ta, tb, firm_t2_mi(t,i), firm_t1_mi(t,i)) > 0.0 )% worker 2 leaves
                    
%                             firm_status_mi(t,i) = 1;
%                             firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                             firm_t2_mi(t,i) = 0;
%                             firm_w2_mi(t,i) = 0.0;
%                         else  % worker stays
                          
%                             firm_status_mi(t,i) = 2;
%                             % check if worker bids up her wage
%                             target_value = Vmarginal_kij(firm_t2_mi(t,i),ta,tb) ; 
%                             current_value=InterpolateWage( wage_w, firm_w2_mi(t-1,i), ...
%                                 W_wkl(:,firm_t2_mi(t,i),firm_t1_mi(t,i)) );
%                             if ( target_value > current_value )% worker renegotiates wage with current firm
                        
%                                 firm_w2_mi(t,i)=InterpolateWage( W_wkl(:,firm_t2_mi(t,i),firm_t1_mi(t,i)), ...
%                                     target_value, wage_w  );
%                             else
%                                 firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                             end
%                             firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                         end
%                     end
                    
%                 else  % no contacts or separations
       
%                     firm_status_mi(t,i) = 2;
%                     firm_w1_mi(t,i) = firm_w1_mi(t-1,i);
%                     firm_w2_mi(t,i) = firm_w2_mi(t-1,i);
%                 end
%             end
            
%             % Check wages are in the bargaining set
%             % make sure wages are inside the bargaining set.
%             if ( firm_status_mi(t,i) == 0 )
 
%             elseif ( firm_status_mi(t,i) == 1 )
             
%                 if (firm_t1_mi(t,i) > 0)% worker is in slot 1
                    
%                    firm_w1_mi(t,i) = min( max(wage_target_U_solo(firm_t1_mi(t,i)),  firm_w1_mi(t,i)), ...
%                         wage_target_MP_solo(firm_t1_mi(t,i)) );
                     
%                 else   % worker is in slot 2
                     
%                     firm_w2_mi(t,i) = min( max(wage_target_U_solo(firm_t2_mi(t,i)), firm_w2_mi(t,i)), ...
%                         wage_target_MP_solo(firm_t2_mi(t,i)) );
                    
%                 end
%             elseif ( firm_status_mi(t,i) == 2 )
                
%                 firm_w1_mi(t,i) = min( max(wage_target_U_team(firm_t1_mi(t,i), firm_t2_mi(t,i)), ...
%                     firm_w1_mi(t,i)), wage_target_MP_team(firm_t1_mi(t,i), firm_t2_mi(t,i) ) );
%                 firm_w2_mi(t,i) = min( max(wage_target_U_team(firm_t2_mi(t,i), firm_t1_mi(t,i)), ...
%                     firm_w2_mi(t,i)), wage_target_MP_team(firm_t2_mi(t,i), firm_t1_mi(t,i) ) );
       
%             end
             
%         end
         
%         if r.state(t,i)>.999
%             i
%         end
%     end
    
     
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Worker sim
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%     quantiles=zeros(nquantiles,nh,nh);
    
%     for l=1:nh
%         for k=1:nh
%             I=find(firm_status_mi(n_months,:) == 2 & firm_t1_mi(n_months,:) == k  & firm_t2_mi(n_months,:) == l );
            
%             if isempty(I)==0
%                 quantiles(:,k,l)=quantile(firm_w1_mi(n_months,I),nquantiles);
%             end
             
%         end
%     end
    
    
    
%     for k=1:nh
%         I=find(firm_status_mi(n_months,:) == 1 & (firm_t1_mi(n_months,:) == k  | firm_t2_mi(n_months,:) == k) );
         
%         if isempty(I)==0
%             quantiles(:,k,nh+1)=quantile(firm_w1_mi(n_months,I)+firm_w2_mi(n_months,I),nquantiles);
%         end
         
%     end
    
    
%     %Initialize worker states
%     type_mi=zeros(n_months,n_workers);
%     status_mi=zeros(n_months,n_workers);
%     age_mi=zeros(n_months,n_workers);
%     tenure_mi=zeros(n_months,n_workers);
%     move_mi=zeros(n_months,n_workers);
%     wage_mi=zeros(n_months,n_workers);
%     co_type_mi=zeros(n_months,n_workers);
%     co_tenure_mi=zeros(n_months,n_workers);
%     co_wage_mi=zeros(n_months,n_workers);
    
    
%     t=1;
%     % Initial conditions
%     for i=1:n_workers
        
%         status_mi(t,i)=draw_CDF_1d( cdf_worker_initial_state, r.state(t,i) );
%         status_mi(t,i) = status_mi(t,i) - 1;  
%         age_mi(t,i)=draw_CDF_1d( cdf_age, r.age(i));
       
%         if ( status_mi(t,i) == 0 )   % Worker is Unemployed
    
%             type_mi(t,i)=draw_CDF_1d( cdf_worker_unempl, r.type(t,i) );
%             tenure_mi(t,i) = 1 ; 
%             move_mi(t,i)   = 1 ;  
%             wage_mi(t,i)   = 0.00;
%             co_type_mi(t,i) = 0;
%             co_tenure_mi(t,i) = 0;
%             co_wage_mi(t,i) = 0.00;
            
%         elseif ( status_mi(t,i) == 1 )   % worker is solo employed
        
%             type_mi(t,i)=draw_CDF_1d( cdf_solo, r.type(t,i) ); 
%             k = 1 + floor( r.tie_break(t,i)*nquantiles );
%             wage_mi(t,i) = quantiles( k ,type_mi(t,i), nh+1 );
%             tenure_mi(t,i) = age_mi(t,i) / 10 ;    
%             move_mi(t,i)   = 0       ;  
%             co_type_mi(t,i) = 0;
%             co_tenure_mi(t,i) = 0;
%             co_wage_mi(t,i) = 0.00;
            
%         else   % Worker is part of a team
         
%             [type_mi(t,i), co_type_mi(t,i)]=draw_CDF_2d( cdf_team, nh, r.type(t,i) ); 
%             k = 1 + floor( r.tie_break(t,i)*nquantiles );
%             wage_mi(t,i) = quantiles( k ,type_mi(t,i), co_type_mi(t,i) );
%             tenure_mi(t,i) = age_mi(t,i) / 10 ;    
%             move_mi(t,i)   = 0     ;      
%             target_value = (1.00-gamma) * U_k( co_type_mi(t,i) ) ...
%                 + gamma * Vmarginal_ki( co_type_mi(t,i), type_mi(t,i) ) ; 
%             k = 1 + floor( r.death(t,i)*nquantiles );
%             co_wage_mi(t,i) = quantiles( k, co_type_mi(t,i) ,type_mi(t,i) ) ; 
%             co_tenure_mi(t,i) = age_mi(t,i) / 10 ;   
%         end
%     end
     
    
%     for i=1:n_workers
%         for t=2:n_months
             
%             if ( (r.death(t,i) <= sigma) | age_mi(t-1,i) >= max_age )   % worker dies and is replaced
          
%                 type_mi(t,i)=draw_CDF_1d( cdf_newborn, r.type(t,i) );
%                 status_mi(t,i) = 0;
%                 age_mi(t,i) = 0;
%                 tenure_mi(t,i) = 0;
%                 move_mi(t,i) = 0;
%                 wage_mi(t,i) = 0.00;
%                 co_type_mi(t,i) = 0;
%                 co_tenure_mi(t,i) = 0;
%                 co_wage_mi(t,i) = 0.00;
%             else % Worker did not die.
                
%                 age_mi(t,i) = age_mi(t-1,i) + 1;  % update worker age
                
%                 if ( status_mi(t-1,i) == 0 )   % unemployed worker
    
%                     if ( r.learn(t,i) < pr.gu_down_k(type_mi(t-1,i)) )   % lose human capital
                        
%                         type_mi(t,i) = type_mi(t-1,i) - 1;
%                     else % human capital stays the same
                     
%                         type_mi(t,i) = type_mi(t-1,i);
%                     end
                    
%                     % unemployed worker evernts:
%                     if ( r.event(t,i) < lambda0 )   % unemployed worker contacts a firm
                   
%                         [status_mi(t,i), co_type_mi(t,i), tenure_mi(t,i), move_mi(t,i),  wage_mi(t,i), co_tenure_mi(t,i), co_wage_mi(t,i)]=...
%                             WorkerMeetsFirm( r.state(t,i), r.type(t,i), r.tie_break(t,i), cdf_firm_state, cdf_solo, cdf_team,...
%                             wage_target_MP_solo, quantiles, type_mi(t,i), 0, 0, tenure_mi(t-1,i),0, 0.00, 0.00,...
%                             gamma,nh,nquantiles,h0u_k, h1u_ik, h2u_ijk, r_ijk, h0e_k, h1e_ik, h2e_ijk, h0ee_kl,h1ee_ikl, h2ee_ijkl, wage_w,W_wkl,W_wk, U_k,Vmarginal_k, Vmarginal_ki, Vmarginal_kij );
                        
%                     else
                     
%                         status_mi(t,i) = 0;
%                         tenure_mi(t,i) = tenure_mi(t-1,i) + 1;
%                         move_mi(t,i) = 0;
%                         wage_mi(t,i) = 0.00;
%                         co_type_mi(t,i) = 0;
%                         co_tenure_mi(t,i) = 0;
%                         co_wage_mi(t,i) = 0.00;
%                     end
                    
%                 elseif ( status_mi(t-1,i) == 1 )   % solo employed worker
                    
%                     if ( r.learn(t,i) >= pr.g_stay_k(type_mi(t-1,i)) )   % gain human capital
                      
%                         type_mi(t,i) = type_mi(t-1,i) + 1;
%                     else % human capital stays the same
               
%                         type_mi(t,i) = type_mi(t-1,i);
%                     end
                    
%                     % solo worker events
%                     if ( r.event(t,i) <= delta )   % separated to unemployment
                        
%                         status_mi(t,i) = 0;
%                         tenure_mi(t,i) = 1;
%                         move_mi(t,i) = 1;
%                         wage_mi(t,i) = 0.00;
%                         co_type_mi(t,i) = 0;
%                         co_tenure_mi(t,i) = 0;
%                         co_wage_mi(t,i) = 0.00;
%                     elseif ( (r.event(t,i) > delta) & ( r.event(t,i) <= delta + lambda1  ) )
                       
%                         [status_mi(t,i), co_type_mi(t,i), tenure_mi(t,i),    move_mi(t,i), wage_mi(t,i), co_tenure_mi(t,i), co_wage_mi(t,i)]= ...
%                             WorkerMeetsFirm( r.state(t,i), r.type(t,i), r.tie_break(t,i), ...
%                             cdf_firm_state, cdf_solo, cdf_team, wage_target_MP_solo, quantiles, ...
%                             type_mi(t,i), 1, 0, tenure_mi(t-1,i), ...
%                             0, wage_mi(t-1,i), 0.00,...
%                             gamma,nh,nquantiles,h0u_k, h1u_ik, h2u_ijk, r_ijk, h0e_k, h1e_ik, h2e_ijk, h0ee_kl,h1ee_ikl, h2ee_ijkl, wage_w,W_wkl,W_wk, U_k,Vmarginal_k, Vmarginal_ki, Vmarginal_kij );
                        
                        
%                     elseif ( (r.event(t,i) > delta + lambda1) & (r.event(t,i) ...
%                             <= delta + lambda1 + prob_firm_contact_worker ) )
                       
%                         [status_mi(t,i), co_type_mi(t,i), tenure_mi(t,i), move_mi(t,i), ...
%                             wage_mi(t,i), co_tenure_mi(t,i), co_wage_mi(t,i)]=FirmMeetsWorker(  r.state(t,i), r.type(t,i), r.tie_break(t,i),  ...
%                             cdf_worker_state, cdf_worker_unempl, cdf_solo, cdf_team, ...
%                             type_mi(t,i), 1, 0, tenure_mi(t-1,i), 0, wage_mi(t-1,i), co_wage_mi(t-1,i), ...
%                             gamma,nh,nquantiles,h0u_k, h1u_ik, h2u_ijk, r_ijk, h0e_k, h1e_ik, h2e_ijk, h0ee_kl,h1ee_ikl, h2ee_ijkl, wage_w,W_wkl,W_wk, U_k,Vmarginal_k, Vmarginal_ki, Vmarginal_kij );
%                     else % No contact events this period
                  
%                         status_mi(t,i) = 1;
%                         tenure_mi(t,i) = tenure_mi(t-1,i) + 1;
%                         move_mi(t,i) = 0;
%                         wage_mi(t,i) = wage_mi(t-1,i);
%                         co_type_mi(t,i) = 0;
%                         co_tenure_mi(t,i) = 0;
%                         co_wage_mi(t,i) = 0.00;
%                     end
                    
                    
%                 else % team worker events
%                    if ( r.learn(t,i) >= pr.g_stay_kl(type_mi(t-1,i), co_type_mi(t-1,i) ) )   % gain human capital
                         
%                         type_mi(t,i) = type_mi(t-1,i) + 1;
%                     else % human capital stays the same 
%                         type_mi(t,i) = type_mi(t-1,i);
%                     end
%                     % update co-worker human capital
                    
%                     if ( r.co_learn(t,i) >= pr.g_stay_kl(co_type_mi(t-1,i), type_mi(t-1,i) ) )   % co worker gain human capital
                        
%                         co_type_mi(t,i) = co_type_mi(t-1,i) + 1;
%                     else % human capital stays the same
                     
%                         co_type_mi(t,i) = co_type_mi(t-1,i);
%                     end
%                     co_type0 = co_type_mi(t,i);
                    
%                     % team worker events
%                     if ( r.event(t,i) <= delta )   % separated to unemployment
                      
%                         status_mi(t,i) = 0;
%                         tenure_mi(t,i) = 1;
%                         move_mi(t,i) = 1;
%                         wage_mi(t,i) = 0.00;
%                         co_type_mi(t,i) = 0;
%                         co_tenure_mi(t,i) = 0;
%                         co_wage_mi(t,i) = 0.00;
%                     elseif ( (r.event(t,i) > delta ) ...
%                             & (r.event(t,i) <= delta+delta ) )
                 
%                         status_mi(t,i) = 1;
%                         tenure_mi(t,i) = tenure_mi(t-1,i) + 1;
%                         move_mi(t,i) = 0;
%                         wage_mi(t,i) = wage_mi(t-1,i);
%                         co_type_mi(t,i) = 0;
%                         co_tenure_mi(t,i) = 0;
%                         co_wage_mi(t,i) = 0.00;
%                     elseif ( (r.event(t,i) > delta + delta ) ...
%                             & (r.event(t,i) <= delta + delta + lambda1 ) )
             
%                         [status_mi(t,i), co_type_mi(t,i), tenure_mi(t,i), ...
%                             move_mi(t,i),  wage_mi(t,i), co_tenure_mi(t,i), co_wage_mi(t,i)]=coWorkerMeetsFirm( r.state(t,i), r.type(t,i), cdf_firm_state, ...
%                             cdf_solo, cdf_team, type_mi(t,i), co_type0, tenure_mi(t-1,i), co_tenure_mi(t-1,i), ...
%                             wage_mi(t-1,i), co_wage_mi(t-1,i), ...
%                             gamma,nh,nquantiles,h0u_k, h1u_ik, h2u_ijk, r_ijk, h0e_k, h1e_ik, h2e_ijk, h0ee_kl,h1ee_ikl, h2ee_ijkl, wage_w,W_wkl,W_wk, U_k,Vmarginal_k, Vmarginal_ki, Vmarginal_kij);
%                     elseif ( (r.event(t,i) > delta + delta + lambda1 ) ...
%                             & (r.event(t,i) <= delta + delta + lambda1 + lambda1 ) )
%                         % Worker receives an offer
                      
%                         [status_mi(t,i), co_type_mi(t,i), tenure_mi(t,i), move_mi(t,i),   wage_mi(t,i), co_tenure_mi(t,i), co_wage_mi(t,i)]=...
%                             WorkerMeetsFirm( r.state(t,i), r.type(t,i), r.tie_break(t,i), ...
%                             cdf_firm_state, cdf_solo, cdf_team, wage_target_MP_solo,  quantiles, ...
%                             type_mi(t,i), 2, co_type0, tenure_mi(t-1,i), ...
%                             co_tenure_mi(t-1,i), wage_mi(t-1,i), co_wage_mi(t-1,i),...
%                             gamma,nh,nquantiles,h0u_k, h1u_ik, h2u_ijk, r_ijk, h0e_k, h1e_ik, h2e_ijk, h0ee_kl,h1ee_ikl, h2ee_ijkl, wage_w,W_wkl,W_wk, U_k,Vmarginal_k, Vmarginal_ki, Vmarginal_kij );
                        
%                     elseif ( (r.event(t,i) > delta + delta + lambda1 + lambda1 ) ...
%                             & (r.event(t,i) <= delta+ delta + lambda1 ...
%                             + lambda1 + prob_firm_contact_worker) )
                         
%                         % Frim contacts another worker.
%                         [status_mi(t,i), co_type_mi(t,i), tenure_mi(t,i), move_mi(t,i), ...
%                             wage_mi(t,i), co_tenure_mi(t,i), co_wage_mi(t,i)]= FirmMeetsWorker(  r.state(t,i), r.type(t,i), r.tie_break(t,i),  ...
%                             cdf_worker_state, cdf_worker_unempl, cdf_solo, cdf_team, ...
%                             type_mi(t,i), 2, co_type0, tenure_mi(t-1,i), co_tenure_mi(t-1,i), ...
%                             wage_mi(t-1,i), co_wage_mi(t-1,i),...
%                             gamma,nh,nquantiles,h0u_k, h1u_ik, h2u_ijk, r_ijk, h0e_k, h1e_ik, h2e_ijk, h0ee_kl,h1ee_ikl, h2ee_ijkl, wage_w,W_wkl,W_wk, U_k,Vmarginal_k, Vmarginal_ki, Vmarginal_kij );
                        
%                     else % No contact events this period
                         
%                         status_mi(t,i) = 2;
%                         tenure_mi(t,i) = tenure_mi(t-1,i) + 1;
%                         move_mi(t,i) = 0;
%                         wage_mi(t,i) = wage_mi(t-1,i);
%                         co_type_mi(t,i) = co_type0;
%                         co_tenure_mi(t,i) = co_tenure_mi(t-1,i) + 1;
%                         co_wage_mi(t,i) = co_wage_mi(t-1,i);
%                     end
%                 end
                
%                 % Dismissal Stage
%                 % Check if still employed at dismissal stage
%                 if ( status_mi(t,i) == 0 )
%                     status_mi(t,i) = 0;
%                 elseif ( status_mi(t,i) == 1 )
%                     if ( dismissk_k(type_mi(t,i)) > r.tie_break(t,i) )
%                         status_mi(t,i) = 0;
%                         tenure_mi(t,i) = 0;
%                         move_mi(t,i) = 1;
%                         wage_mi(t,i) = 0.00;
%                         co_type_mi(t,i) = 0;
%                         co_tenure_mi(t,i) = 0;
%                         co_wage_mi(t,i) = 0.00;
%                     end
%                     % check whether or not the firm keeps both workers and, if not, who is terminated
%                 else
%                     if (keepkl_kl(type_mi(t,i),co_type_mi(t,i)) > 0.00)
%                         status_mi(t,i) = 2;
%                     elseif ( (dismisskl_kl(type_mi(t,i),co_type_mi(t,i)) >0 ))
%                         status_mi(t,i) = 0;
%                         tenure_mi(t,i) = 0;
%                         move_mi(t,i) = 1;
%                         wage_mi(t,i) = 0.00;
%                         co_type_mi(t,i) = 0;
%                         co_tenure_mi(t,i) = 0;
%                         co_wage_mi(t,i) = 0.00;
%                     elseif ( (dismissk_kl(type_mi(t,i),co_type_mi(t,i)) > r.tie_break(t,i) ))
%                         status_mi(t,i) = 0;
%                         tenure_mi(t,i) = 0;
%                         move_mi(t,i) = 1;
%                         wage_mi(t,i) = 0.00;
%                         co_type_mi(t,i) = 0;
%                         co_tenure_mi(t,i) = 0;
%                         co_wage_mi(t,i) = 0.00;
%                     else  
%                         status_mi(t,i) = 1;
%                         tenure_mi(t,i) = tenure_mi(t,i);
%                         move_mi(t,i) = 0;
%                         wage_mi(t,i) = wage_mi(t,i);
%                         co_type_mi(t,i) = 0;
%                         co_tenure_mi(t,i) = 0;
%                         co_wage_mi(t,i) = 0.00;
%                     end
%                 end
%             end
            
%             % make sure wages are inside the bargaining set.
%             if ( status_mi(t,i) == 0 )
               
%             elseif ( status_mi(t,i) == 1 )
               
%                 wage_mi(t,i) = min( max(wage_target_U_solo(type_mi(t,i)), wage_mi(t,i)), ...
%                     wage_target_MP_solo(type_mi(t,i)) );
%             elseif ( status_mi(t,i) == 2 )
                 
%                 wage_mi(t,i) = min( max(wage_target_U_team(type_mi(t,i), co_type_mi(t,i)), wage_mi(t,i)), ...
%                     wage_target_MP_team(type_mi(t,i),co_type_mi(t,i) ) );
%                 co_wage_mi(t,i) = min( max(wage_target_U_team(co_type_mi(t,i), type_mi(t,i)), co_wage_mi(t,i)), ...
%                     wage_target_MP_team(co_type_mi(t,i),type_mi(t,i) ) );
%             end
             
%         end
%     end
    
     
%     age1c=age_mi;
%     type1c=type_mi;
%     status1c=(status_mi>0);
%     wage1c=wage_mi;
%     coworker1c=co_type_mi;
%     tenure1c=tenure_mi.*(status_mi>0);
%     dur1c=tenure_mi.*(status_mi==0);
%     death1c=(age1c==0);
%     move1c=move_mi;
%     cotenure1c=co_tenure_mi;
%     cowage1c=co_wage_mi;
    
%     age1c(1:t_burn,:)      =[]; %start age 1
%     type1c(1:t_burn,:)     =[]; %Draw type
%     status1c(1:t_burn,:)   =[]; %Start unemployed
%     wage1c(1:t_burn,:)     =[]; %Wage zero for unempl
%     coworker1c(1:t_burn,:) =[]; %Coworker type
%     tenure1c(1:t_burn,:)   =[]; %Tenure
%     death1c(1:t_burn,:)    =[];
%     move1c(1:t_burn,:)     =[];
%     cotenure1c(1:t_burn,:) =[]; %Co Tenure
%     cowage1c(1:t_burn,:)   =[]; %Wage zero for unempl
%     dur1c(1:t_burn,:)      =[]; %duration    
     
 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Wide format data for LEHD regs
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     N=n_workers;
%     T=n_months;
    
%     death1c_annual= death1c 	+ lag(death1c,1) + lag(death1c,2) + lag(death1c,3)+ lag(death1c,4)+ lag(death1c,5)+ lag(death1c,6) ...
%         + lag(death1c,7)+lag(death1c,8)+lag(death1c,9)+lag(death1c,10)+lag(death1c,11)  ;
    
%     cowage1c_annual= cowage1c + lag(cowage1c,1) + lag(cowage1c,2) + lag(cowage1c,3)+ lag(cowage1c,4)+ lag(cowage1c,5)+ lag(cowage1c,6) ...
%         + lag(cowage1c,7)+lag(cowage1c,8)+lag(cowage1c,9)+lag(cowage1c,10)+lag(cowage1c,11) ;
    
%     wage1c_annual= wage1c + lag(wage1c,1) + lag(wage1c,2) + lag(wage1c,3)+ lag(wage1c,4)+ lag(wage1c,5)+ lag(wage1c,6) ...
%         + lag(wage1c,7)+lag(wage1c,8)+lag(wage1c,9)+lag(wage1c,10)+lag(wage1c,11)	 ;
 
%     cotenure1c_f1=lag(cotenure1c,-1); 
    
%     %Quarterly wage
%     wage1c_q4=    wage1c       + lag(wage1c,1) + lag(wage1c,2)  ;
%     wage1c_q3=   lag(wage1c,3) + lag(wage1c,4) + lag(wage1c,5)  ;
%     wage1c_q2=   lag(wage1c,6) + lag(wage1c,7) + lag(wage1c,8)  ;
%     wage1c_q1=   lag(wage1c,9) + lag(wage1c,10)+ lag(wage1c,11) ;

%     cowage1c_q4=    cowage1c       + lag(cowage1c,1) + lag(cowage1c,2)  ;
%     cowage1c_q3=   lag(cowage1c,3) + lag(cowage1c,4) + lag(cowage1c,5)  ;
%     cowage1c_q2=   lag(cowage1c,6) + lag(cowage1c,7) + lag(cowage1c,8)  ;
%     cowage1c_q1=   lag(cowage1c,9) + lag(cowage1c,10)+ lag(cowage1c,11) ;
        
    
%     %Move
%     move1c_q4=    move1c       + lag(move1c,1) + lag(move1c,2)  ;
%     move1c_q3=   lag(move1c,3) + lag(move1c,4) + lag(move1c,5)  ;
%     move1c_q2=   lag(move1c,6) + lag(move1c,7) + lag(move1c,8)  ;
%     move1c_q1=   lag(move1c,9) + lag(move1c,10)+ lag(move1c,11) ;
    
%     %Status
%     status1c_q4=    status1c       + lag(status1c,1) + lag(status1c,2)  ;
%     status1c_q3=   lag(status1c,3) + lag(status1c,4) + lag(status1c,5)  ;
%     status1c_q2=   lag(status1c,6) + lag(status1c,7) + lag(status1c,8)  ;
%     status1c_q1=   lag(status1c,9) + lag(status1c,10)+ lag(status1c,11) ;
    
    
%     id=0;
%     period=0;
    
%     data=zeros(floor(N*(T-t_burn)/12),28); %Wide format, store 23 variables
%     data_row=0;
    
%     for j=1:N
%         period=0   ; %Reset period
%         id    =id+1;
%         for i=12:12:T-t_burn %Step forward annually
            
            
%             if death1c_annual(i,j)==0
%                 data_row=data_row+1; %Advance data row
%                 period=period+1; %Advance period
                
%                 %Stack data in annual wide format
%                 data(data_row,:)=[ id,period, wage1c_annual(i,j),wage1c_q1(i,j),wage1c_q2(i,j),wage1c_q3(i,j), wage1c_q4(i,j),...
%                     cowage1c_annual(i,j), tenure1c(i,j),cotenure1c(i,j),type1c(i,j),coworker1c(i,j),...
%                     move1c_q1(i,j),move1c_q2(i,j),move1c_q3(i,j), move1c_q4(i,j),status1c_q1(i,j),status1c_q2(i,j),status1c_q3(i,j), status1c_q4(i,j),...
%                     age1c(i,j),cotenure1c_f1(i,j),cowage1c_q1(i,j),cowage1c_q2(i,j),cowage1c_q3(i,j), cowage1c_q4(i,j),i,j ];
%             else
%                 %Skip deaths
%                 period=0;
%                 id=id+1;
                
%             end
%         end
%     end
    
%     I_remove=find(sum(data==0,2)==28);
%     data(I_remove,:)=[];
    
%     %Unpack data
    
%     id                  =data(:,1);
%     period              =data(:,2);
%     wage1c_annual       =data(:,3);
%     wage1c_q1           =data(:,4);
%     wage1c_q2           =data(:,5);
%     wage1c_q3           =data(:,6);
%     wage1c_q4           =data(:,7);
%     cowage1c_annual     =data(:,8);
%     tenure1c_annual     =data(:,9);
%     cotenure1c_annual   =data(:,10);
%     type1c_annual       =data(:,11);
%     coworker1c_annual   =data(:,12);
%     move1c_q1           =data(:,13);
%     move1c_q2           =data(:,14);
%     move1c_q3           =data(:,15);
%     move1c_q4           =data(:,16);
%     status1c_q1         =data(:,17);
%     status1c_q2         =data(:,18);
%     status1c_q3         =data(:,19);
%     status1c_q4         =data(:,20);
%     age1c_annual        =data(:,21);
%     cotenure1c_f1_annual=data(:,22);
%     cowage1c_q1         =data(:,23);
%     cowage1c_q2         =data(:,24);
%     cowage1c_q3         =data(:,25);
%     cowage1c_q4         =data(:,26);    
%     i_annual            =data(:,27);
%     j_annual            =data(:,28);
    
%     clearvars data
    
%     id_l1=lag(id,1);
%     id_f1=lag(id,-1);
    
%     tenure1c_annual_l1=lag(tenure1c_annual,1);
%     tenure1c_annual_f1=lag(tenure1c_annual,-1);
    
%     cotenure1c_annual_l1=lag(cotenure1c_annual,1);
 
%     wage1c_annual_l1=lag(wage1c_annual,1);
%     wage1c_annual_f1=lag(wage1c_annual,-1);
    
%     cowage1c_annual_l1=lag(cowage1c_annual,1);
%     cowage1c_annual_f1=lag(cowage1c_annual,-1);
%     cowage1c_annual_f2=lag(cowage1c_annual,-2);
        
%     status1c_annual=status1c_q1+status1c_q2+status1c_q3+status1c_q4;
%     move1c_annual=move1c_q1+move1c_q2+move1c_q3+move1c_q4;
%     move1c_annual_f1=lag(move1c_annual,-1);
%     move1c_annual_l1=lag(move1c_annual,1);
        
%     %Winsorize variables
%     wins_lb_cut=.5;
%     wins_ub_cut=99.5;
    
    
%     %Scaling
%     ucutoff= mean(wage1c_annual)*1000/45600; %1k cutoff in model
%     lp1    = 0;% 
    
%     i_zero_q=(wage1c_q1<=ucutoff | wage1c_q2<=ucutoff | wage1c_q3<=ucutoff | wage1c_q4<=ucutoff); %Job loss in census
%     i_zero_q_f1=lag(i_zero_q,-1);
    
%     full_yr =(wage1c_q1>ucutoff  & wage1c_q2>ucutoff  & wage1c_q3>ucutoff  & wage1c_q4>ucutoff); %Full year employment in census
%     full_yr_f1=lag(full_yr,-1);
%     full_yr_l1=lag(full_yr,1);
    
%     %Stable worker 
%     i_stable=(id==id_f1 & cotenure1c_annual>=13 ) ...
%                 & (cowage1c_q1>ucutoff  & cowage1c_q2>ucutoff  & cowage1c_q3>ucutoff  & cowage1c_q4>ucutoff); %Stable coworker defintion in census
                
%     i_stable_l1=lag(i_stable,1); %Year before layoff         
     
      
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % w_{i,t+2} = b0+b1*w_{i,t}+b2*w_{-i,t}+eps for EUE at date t
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%     [I_samp]=find(id==id_l1 & id==id_f1 & tenure1c_annual_l1>=12 & i_stable_l1==1 & full_yr_l1==1  & cowage1c_annual_l1>0 &  wage1c_annual_l1>0  ...
%         & i_zero_q==1 & full_yr_f1==1 & cowage1c_annual_f1>0 &  wage1c_annual_f1>0 & age1c_annual>=3*12 & age1c_annual<45*12);
      
%     %Pooled
%     lw=winsorize(log(wage1c_annual_l1(I_samp)+lp1),wins_lb_cut,wins_ub_cut);
%     lc=winsorize(log(cowage1c_annual_l1(I_samp)+lp1),wins_lb_cut,wins_ub_cut);
%     lf=winsorize(log(wage1c_annual_f1(I_samp)+lp1),wins_lb_cut,wins_ub_cut);
%     lfc=winsorize(log(cowage1c_annual_f1(I_samp)+lp1),wins_lb_cut,wins_ub_cut);
      
%     X=[ones(length(I_samp),1),lw,lc];
%     Y=lf;
    
%     [~,nX]=size(X);
%     B_own_wage=(X'*X)\X'*Y;   
%     eps=Y-X*B_own_wage;
%     sig_sq=sum(eps.^2)/length(I_samp);
%     std_error_own_wage=(sig_sq*eye(nX)/(X'*X))^.5;
    
%     Y=lfc;
%     B_co_wage=(X'*X)\X'*Y;  
%     eps=Y-X*B_co_wage;
%     sig_sq=sum(eps.^2)/length(I_samp);
%     std_error_co_wage=(sig_sq*eye(nX)/(X'*X))^.5;
     
%     %Below the mean
%     I_samp_below=I_samp(wage1c_annual_l1(I_samp)<cowage1c_annual_l1(I_samp));
    
%     lw=winsorize(log(wage1c_annual_l1(I_samp_below)+lp1),wins_lb_cut,wins_ub_cut);
%     lc=winsorize(log(cowage1c_annual_l1(I_samp_below)+lp1),wins_lb_cut,wins_ub_cut);
%     lf=winsorize(log(wage1c_annual_f1(I_samp_below)+lp1),wins_lb_cut,wins_ub_cut);
%     lfc=winsorize(log(cowage1c_annual_f1(I_samp_below)+lp1),wins_lb_cut,wins_ub_cut);
    
%     X=[ones(length(I_samp_below),1),lw,lc];
%     Y=lf;
    
%     [~,nX]=size(X);
%     B_own_wage_below=(X'*X)\X'*Y; 
%     eps=Y-X*B_own_wage_below;
%     sig_sq=sum(eps.^2)/length(I_samp_below);
%     std_error_own_wage_below=(sig_sq*eye(nX)/(X'*X))^.5;
    
%     Y=lfc;
%     B_co_wage_below=(X'*X)\X'*Y;  
%     eps=Y-X*B_co_wage_below;
%     sig_sq=sum(eps.^2)/length(I_samp_below);
%     std_error_co_wage_below=(sig_sq*eye(nX)/(X'*X))^.5;
     
%     %Above the mean
%     I_samp_above=I_samp(wage1c_annual_l1(I_samp)>=cowage1c_annual_l1(I_samp));
    
%     lw=winsorize(log(wage1c_annual_l1(I_samp_above)+lp1),wins_lb_cut,wins_ub_cut);
%     lc=winsorize(log(cowage1c_annual_l1(I_samp_above)+lp1),wins_lb_cut,wins_ub_cut);
%     lf=winsorize(log(wage1c_annual_f1(I_samp_above)+lp1),wins_lb_cut,wins_ub_cut);
%     lfc=winsorize(log(cowage1c_annual_f1(I_samp_above)+lp1),wins_lb_cut,wins_ub_cut);
    
%     X=[ones(length(I_samp_above),1),lw,lc];
%     Y=lf;
    
%     [~,nX]=size(X);
%     B_own_wage_above=(X'*X)\X'*Y;  
%     eps=Y-X*B_own_wage_above;
%     sig_sq=sum(eps.^2)/length(I_samp_above);
%     std_error_own_wage_above=(sig_sq*eye(nX)/(X'*X))^.5;
    
%     Y=lfc;
%     B_co_wage_above=(X'*X)\X'*Y; 
%     eps=Y-X*B_co_wage_above;
%     sig_sq=sum(eps.^2)/length(I_samp_above);
%     std_error_co_wage_above=(sig_sq*eye(nX)/(X'*X))^.5;
      
%     N_below=length(I_samp_below);
%     N_above=length(I_samp_above);
     
%     log_wage_growth_l1_f1=winsorize(log(wage1c_annual_f1(I_samp )+lp1),wins_lb_cut,wins_ub_cut)-winsorize(log(wage1c_annual_l1(I_samp )+lp1),wins_lb_cut,wins_ub_cut);
  
%     eue_mean_log_wage_growth_l1_f1=mean(log_wage_growth_l1_f1);
 
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % EU_{i,t+1} = b0+b1*|w_{i,t}-w_{-i,t}|+b1*w_{i,t}+eps for E at date t
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
%     I_eu=         (id==id_f1 & i_stable==1 & full_yr==1 & cowage1c_annual>0 & wage1c_annual>0 ...
%         &   i_zero_q_f1==1 &  move1c_annual_f1>=1 & age1c_annual>=3*12 & age1c_annual<45*12);
    
%     [I_samp]=find(id==id_f1 & i_stable==1 & full_yr==1 & cowage1c_annual>0  & wage1c_annual>0 ...
%                                                   &  age1c_annual>=3*12 & age1c_annual<45*12);
     
%     %Pooled
%     lw=winsorize(log(wage1c_annual(I_samp)+lp1),wins_lb_cut,wins_ub_cut);
%     lg=winsorize(abs(log(cowage1c_annual(I_samp)+lp1)-log(wage1c_annual(I_samp))+lp1),wins_lb_cut,wins_ub_cut);
     
%     %Compute regression coefficients
%     X=[ones(length(I_samp),1),lg,lw];
%     Y= I_eu(I_samp) ;
    
%     [~,nX]=size(X);
%     B_eu=(X'*X)\X'*Y ;  
%     eps=Y-X*B_eu;
%     sig_sq=sum(eps.^2)/length(I_samp);
%     std_error_eu=(sig_sq*eye(nX)/(X'*X))^.5;
     
%     %Below
%     I_samp_below=I_samp(wage1c_annual(I_samp)<cowage1c_annual(I_samp));
    
%     lw=winsorize(log(wage1c_annual(I_samp_below)+lp1),wins_lb_cut,wins_ub_cut);
%     lg=winsorize(abs(log(cowage1c_annual(I_samp_below)+lp1)-log(wage1c_annual(I_samp_below))+lp1),wins_lb_cut,wins_ub_cut);
     
%     %Compute regression coefficients
%     X=[ones(length(I_samp_below),1),lg,lw];
%     Y= I_eu(I_samp_below) ;
    
%     [~,nX]=size(X);
%     B_eu_below=(X'*X)\X'*Y ; %For earnings below their coworker at t, eu_{i,t+1} = b0+b1*|w_{i,t}-w_{-i,t}|+b2*w_{i,t}+eps
%     eps=Y-X*B_eu_below;
%     sig_sq=sum(eps.^2)/length(I_samp_below);
%     std_error_eu_below=(sig_sq*eye(nX)/(X'*X))^.5;
    
%     %Above
%     I_samp_above=I_samp(wage1c_annual(I_samp)>=cowage1c_annual(I_samp));
%     lw=winsorize(log(wage1c_annual(I_samp_above)+lp1),wins_lb_cut,wins_ub_cut);
%     lg=winsorize(abs(log(cowage1c_annual(I_samp_above)+lp1)-log(wage1c_annual(I_samp_above))+lp1),wins_lb_cut,wins_ub_cut);
    
%     %Compute regression coefficients
%     X=[ones(length(I_samp_above),1),lg,lw];
%     Y= I_eu(I_samp_above) ;
    
%     [~,nX]=size(X);
%     B_eu_above=(X'*X)\X'*Y ; %For those earning more than coworkers at t, eu_{i,t+1} = b0+b1*|w_{i,t}-w_{-i,t}|+b2*w_{i,t}+eps
%     eps=Y-X*B_eu_above;
%     sig_sq=sum(eps.^2)/length(I_samp_above);
%     std_error_eu_above=(sig_sq*eye(nX)/(X'*X))^.5;
     
%     N_below_eu=length(I_samp_below);
%     N_above_eu=length(I_samp_above);
    
%     %Heat map -- bins of E
%     focal_bin=quantile( wage1c_annual (I_samp),[0:.1:1]') ; 
%     focal_bin(1)=0; 
%     co_bin=quantile( cowage1c_annual (I_samp),[0:.1:1]')  ;
%     co_bin(1)=0; 
    
 
%     eu_n_obs_heatmap=zeros(10,10);
%     eu_heatmap=zeros(10,10);
%     for i=1:10
%         for j=1:10
%         eu_idx=(id==id_f1 & i_stable==1 & full_yr==1 & cowage1c_annual>0  & wage1c_annual>0 ...
%                                                   &  age1c_annual>=3*12 & age1c_annual<45*12) & ...
%                     (wage1c_annual  >focal_bin(i) & wage1c_annual  <=focal_bin(i+1) ...
%                             & cowage1c_annual  >co_bin(j) & cowage1c_annual  <=co_bin(j+1)); %Everyone
%         eu_n_obs_heatmap(i,j) = sum(eu_idx) ;
%         eu_heatmap(i,j)       = mean(I_eu(eu_idx)) ; %Pull elements where this is true
%         end
%     end    
    
%     tot_e_n_obs=length(I_samp);
%     e_density_heatmap=eu_n_obs_heatmap/tot_e_n_obs;
    
%     eu_avg=sum(sum(I_eu))/length(I_samp);
          
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % JJ_{i,t+1} = b0+b1*|w_{i,t}-w_{-i,t}|+b1*w_{i,t}+eps for E at date t
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%     I_j_to_j=     (id==id_f1 &   i_stable==1 & full_yr==1  & cowage1c_annual>0  & wage1c_annual>0 ...
%         &  move1c_annual_f1>=1 & i_zero_q_f1==0  & age1c_annual>=3*12 & age1c_annual<45*12 & cowage1c_annual_f2>0);
    
%     [I_samp]=find(id==id_f1 &   i_stable==1 & full_yr==1  & cowage1c_annual>0  & wage1c_annual>0 ...
%                                                 & age1c_annual>=3*12 & age1c_annual<45*12 & cowage1c_annual_f2>0); %Everyone included in regression, jjs and stayers
       
%     %Pooled
%     lw=winsorize(log(wage1c_annual(I_samp)+lp1),wins_lb_cut,wins_ub_cut);
%     lg=winsorize(abs(log(cowage1c_annual(I_samp)+lp1)-log(wage1c_annual(I_samp))+lp1),wins_lb_cut,wins_ub_cut);
     
%     %Compute regression coefficients
%     X=[ones(length(I_samp),1),lg,lw];
%     Y= I_j_to_j(I_samp) ;
    
%     [~,nX]=size(X);
%     B_j_to_j=(X'*X)\X'*Y ;
%     eps=Y-X*B_j_to_j;
%     sig_sq=sum(eps.^2)/length(I_samp);
%     std_error_j_to_j=(sig_sq*eye(nX)/(X'*X))^.5;
     
%     %Below
%     I_samp_below=I_samp(wage1c_annual(I_samp)<cowage1c_annual(I_samp));
    
%     lw=winsorize(log(wage1c_annual(I_samp_below)+lp1),wins_lb_cut,wins_ub_cut);
%     lg=winsorize(abs(log(cowage1c_annual(I_samp_below)+lp1)-log(wage1c_annual(I_samp_below))+lp1),wins_lb_cut,wins_ub_cut);
     
%     %Compute regression coefficients
%     X=[ones(length(I_samp_below),1),lg,lw];
%     Y= I_j_to_j(I_samp_below) ;
    
%     [~,nX]=size(X);
%     B_j_to_j_below=(X'*X)\X'*Y ; 
%     eps=Y-X*B_j_to_j_below;
%     sig_sq=sum(eps.^2)/length(I_samp_below);
%     std_error_j_to_j_below=(sig_sq*eye(nX)/(X'*X))^.5;
    
%     %Above
%     I_samp_above=I_samp(wage1c_annual(I_samp)>=cowage1c_annual(I_samp));
%     lw=winsorize(log(wage1c_annual(I_samp_above)+lp1),wins_lb_cut,wins_ub_cut);
%     lg=winsorize(abs(log(cowage1c_annual(I_samp_above)+lp1)-log(wage1c_annual(I_samp_above))+lp1),wins_lb_cut,wins_ub_cut);
     
%     %Compute regression coefficients
%     X=[ones(length(I_samp_above),1),lg,lw];
%     Y= I_j_to_j(I_samp_above) ;
    
%     [~,nX]=size(X);
%     B_j_to_j_above=(X'*X)\X'*Y ;  
%     eps=Y-X*B_j_to_j_above;
%     sig_sq=sum(eps.^2)/length(I_samp_above);
%     std_error_j_to_j_above=(sig_sq*eye(nX)/(X'*X))^.5;
    
%     N_below_j_to_j=length(I_samp_below);
%     N_above_j_to_j=length(I_samp_above);
     
%     jj_n_obs_heatmap=zeros(10,10);
%     jj_heatmap=zeros(10,10);
%     for i=1:10
%         for j=1:10
%         j_to_j_idx=(id==id_f1 &   i_stable==1 & full_yr==1  & cowage1c_annual>0  & wage1c_annual>0 ...
%                                                 & age1c_annual>=3*12 & age1c_annual<45*12 & cowage1c_annual_f2>0) & ...
%                     (wage1c_annual  >focal_bin(i) & wage1c_annual  <=focal_bin(i+1) ...
%                             & cowage1c_annual  >co_bin(j) & cowage1c_annual  <=co_bin(j+1));
%         jj_n_obs_heatmap(i,j) = sum(j_to_j_idx) ;
%         jj_heatmap(i,j)= mean(I_j_to_j(j_to_j_idx)) ; %Pull elements where this is true
%         end
%     end
  
%     wage1c_q1_f1=lag(wage1c_q1,-1); %Year before layoff
%     wage1c_q2_f1=lag(wage1c_q2,-1); %Year before layoff
     
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Wage distribution percentiles, annual
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     I_pos_annual=find(id==id_f1 &   i_stable==1 & full_yr==1 & cowage1c_annual>0 ...
%         &  age1c_annual>=3*12 & age1c_annual<45*12);
    
%     p10_wage_all_annual=prctile(wage1c_annual(I_pos_annual),10);
%     p90_wage_all_annual=prctile(wage1c_annual(I_pos_annual),90);
    
%     wage_p90p10_all_a     =p90_wage_all_annual/p10_wage_all_annual;
%     log_wage_p90p10_all_a =log(p90_wage_all_annual)-log(p10_wage_all_annual);    
%     mean_wage_a           =mean(wage1c_annual(I_pos_annual)) ;
    
%     mean_log_wage_a=mean(log(wage1c_annual(I_pos_annual)));
%     var_log_wage_a=var(log(wage1c_annual(I_pos_annual)));
    
%     I_pos_young_annual=find(id==id_f1 &   i_stable==1 & full_yr==1 & cowage1c_annual>0 ...
%         & age1c_annual>=3*12 & age1c_annual<8*12); %24 to 28
    
%     p10_wage_young_annual=prctile(wage1c_annual(I_pos_young_annual),10);
%     p90_wage_young_annual=prctile(wage1c_annual(I_pos_young_annual),90);
    
%     wage_p90p10_young_a     =p90_wage_young_annual/p10_wage_young_annual ;
%     log_wage_p90p10_young_a =log(p90_wage_young_annual)-log(p10_wage_young_annual) ;%%
%     mean_wage_young_a       =mean(wage1c_annual(I_pos_young_annual));
    
%     [I_w_age_y_annual]=find( id==id_f1 &   i_stable==1 & full_yr==1 & cowage1c_annual>0 ...
%           & age1c_annual>=3*12 & age1c_annual<8*12   ); % 24-28
%     [I_w_age_o_annual]=find( id==id_f1 &   i_stable==1 & full_yr==1 & cowage1c_annual>0 ...
%         & age1c_annual>=29*12 & age1c_annual<34*12   ); % 50-54
    
%     mean_wage_y_a=mean(wage1c_annual(I_w_age_y_annual));
%     mean_wage_o_a=mean(wage1c_annual(I_w_age_o_annual));
    
%     mean_log_wage_y_a=mean(winsorize(log(wage1c_annual(I_w_age_y_annual)),wins_lb_cut,wins_ub_cut)); %Mean of log wage, young
%     mean_log_wage_o_a=mean(winsorize(log(wage1c_annual(I_w_age_o_annual)),wins_lb_cut,wins_ub_cut)); % " old
        
%     var_log_wage_y_a=var(winsorize(log(wage1c_annual(I_w_age_y_annual)),wins_lb_cut,wins_ub_cut)); %Var of log wage, young
%     var_log_wage_o_a=var(winsorize(log(wage1c_annual(I_w_age_o_annual)),wins_lb_cut,wins_ub_cut)); %" old
           
%     wage_ratio_old_young_a =mean_wage_o_a/mean_wage_y_a;
%     log_wage_ratio_old_young_a =mean_log_wage_o_a-mean_log_wage_y_a; %Log wage ratio
    
%     [I_pos_growth]=find( wage1c_annual>0 & wage1c_annual_l1>0  & age1c_annual>=3*12  & age1c_annual<45*12); %Stayer, growth
    
%     ln_wage_growth_wins1_pos_annual=winsorize(log(wage1c_annual(I_pos_growth)+lp1)-log(wage1c_annual_l1(I_pos_growth)+lp1),wins_lb_cut,wins_ub_cut);
    
%     mean_ln_wage_change_a =mean(ln_wage_growth_wins1_pos_annual);
%     std_dev_ln_wage_change_a =std(ln_wage_growth_wins1_pos_annual);
    
%     ur_young=1-mean(status1c(  age1c_annual>=3*12 & age1c_annual<4*12)); %Unemployment rate of young
%     ur_old=1-mean(status1c(  age1c_annual>=33*12 & age1c_annual<34*12)); %Unemployment rate of old
  
%     clearvars lw lc lf lfc lw lg I_samp_below I_samp_above I_samp Y X eps *_annual* *_q* full_* id* period
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Correlation of EU and tenure
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %Correlation of wage and human capital
%     wage1c_l1   =lag(wage1c,1);
%     death1c_l1=lag(death1c,1);
%     status1c_l1=lag(status1c,1);
%     tenure1c_l1=lag(tenure1c,1);
%     I_eu=(status1c_l1==1 & status1c ==0 & death1c_l1 ==0 & death1c ==0  & age1c>=3*12 & age1c<45*12);
%     [I_samp_e]=find( wage1c_l1 >0 & status1c_l1==1 &  death1c_l1 ==0 & death1c ==0   & age1c>=3*12 & age1c<45*12 );
    
%     I_corr_eu_tenure=corr(I_eu(I_samp_e),tenure1c_l1(I_samp_e)/12);
    
%     I_corr_eu_hc=corr(I_eu(I_samp_e),type1c(I_samp_e) );
    
%     clearvars I_eu I_samp_e
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Duration rep rate hazard
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%     dur1c_l1=lag(dur1c,1)	;
%     [I_haz]=find(wage1c>0 & status1c==1   &   status1c_l1==0  &   death1c_l1==0 & death1c==0  & age1c>=3*12  & age1c<45*12);
    
%     B_haz=corr(log(wage1c(I_haz)+lp1),4*dur1c_l1(I_haz));
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Transition rates
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [I_empl]=find(status1c_l1==1  & age1c>=3*12 & age1c<45*12);
    
%     [I_eu]=find(status1c_l1==1 & status1c ==0 & death1c_l1 ==0 & death1c ==0  & age1c>=3*12 & age1c<45*12);
    
%     eu_rate=length(I_eu)./length(I_empl);
    
%     [I_unempl]=find(status1c_l1==0  & age1c>=3*12 & age1c<45*12);
    
%     [I_ue]=find(status1c_l1==0 & status1c ==1 & death1c_l1 ==0 & death1c ==0  & age1c>=3*12  & age1c<45*12);
    
%     ue_rate=length(I_ue)./length(I_unempl);
    
%     [I_ee]=find(status1c_l1==1 & status1c ==1 & move1c==1 & death1c_l1 ==0 & death1c ==0  & age1c>=3*12 & age1c<45*12);
    
%     ee_rate=length(I_ee)./length(I_empl);
    
%     unempl_rate=sum(sum(udist));
    
%     sole_prop_rate=sum(sum(solodist));
    
%     rep_rate= mean(b(type1c(I_eu)))./mean(wage1c_l1(I_eu))  ;
    
%     clearvars I_empl I_eu I_ee
    
%     %Pre/Post layoff wag ratio
%     ecutoff= mean(mean(wage1c))*12*1600/49340; %1600k cutoff for MONTHLY earnings (400 per week)
    
%     status1c_l2  =lag(status1c,2);
%     status1c_l3  =lag(status1c,3);
%     wage1c_l12   =lag(wage1c,12);
%     status1c_l12 =lag(status1c,12);
%     death1c_l12ma=(1/12)*(lag(death1c,1) + lag(death1c,2) + lag(death1c,3)+ lag(death1c,4)+ lag(death1c,5)+ lag(death1c,6) ...
%         + lag(death1c,7)+lag(death1c,8)+lag(death1c,9)+lag(death1c,10)+lag(death1c,11)+lag(death1c,12)	);
    
%     [I_disp]=find( wage1c_l12>0 & status1c_l12==1 & wage1c_l12>ecutoff & (status1c_l1==0 | status1c_l2==0 | status1c_l3==0) ...
%         &   death1c_l12ma==0  & age1c>=3*12  & age1c<45*12 ); %Displaced
    
    
%     wratio=winsorize( wage1c(I_disp)./wage1c_l12(I_disp) ,wins_lb_cut,wins_ub_cut);
    
%     wage_drop=mean(wratio); %Wage drop after displacement
    
%     clearvars I_disp wratio
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Wage distribution percentiles
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     status1c_long=reshape(status1c,1,N*(T-t_burn   ))  ;
%     wage1c_long  =reshape(wage1c  ,1,N*(T-t_burn   ))  ;
%     cowage1c_long=reshape(cowage1c,1,N*(T-t_burn   ))  ;
%     age1c_long   =reshape(age1c   ,1,N*(T-t_burn   ))  ;
    
%     I_emp=find(status1c_long==1 & age1c_long>=3*12 & age1c_long<45*12);
    
%     p10_wage_all=prctile(wage1c_long(I_emp),10);
%     p90_wage_all=prctile(wage1c_long(I_emp),90);
    
%     wage_p90p10_all=p90_wage_all/p10_wage_all;
%     mean_wage      =mean(wage1c_long(I_emp)) ;
    
%     I_emp_young=find(status1c_long==1  & age1c_long>=3*12 & age1c_long<4*12);
    
%     p10_wage_young=prctile(wage1c_long(I_emp_young),10);
%     p90_wage_young=prctile(wage1c_long(I_emp_young),90);
    
%     wage_p90p10_young=p90_wage_young/p10_wage_young ;
%     mean_wage_young  =mean(wage1c_long(I_emp_young));
    
    
%     %JF Rate, Job Loss Rate, and EE Rate by Prior Wage
    
%     percentiles=[10:10:100];
%     eu_rate_by_pctile=zeros(length(percentiles),1);
%     ue_rate_by_pctile=zeros(length(percentiles),1);
%     ee_rate_by_pctile=zeros(length(percentiles),1);
%     for j=1:10
%         i=percentiles(j);
%         pi_wage_all_h=prctile(wage1c_long(I_emp),i); %Isolate wage percentile
%         pi_wage_all_l=prctile(wage1c_long(I_emp),i-10); %Isolate
        
%         %EU
%         [I_empl_2]=find(status1c_l1==1 & wage1c_l1>pi_wage_all_l & wage1c_l1<pi_wage_all_h);
        
%         [I_eu_2]=find(status1c_l1==1 & status1c ==0 & death1c_l1 ==0 & death1c ==0  & wage1c_l1>pi_wage_all_l & wage1c_l1<pi_wage_all_h );
        
%         eu_rate_by_pctile(j)=length(I_eu_2)./length(I_empl_2);
        
%         %EE
%         [I_ee_2]=find(status1c_l1==1 & status1c ==1 & move1c==1 & death1c_l1 ==0 & death1c ==0 & wage1c_l1>pi_wage_all_l & wage1c_l1<pi_wage_all_h );
        
%         ee_rate_by_pctile(j)=length(I_ee_2)./length(I_empl_2);
        
%         %UE
%         [I_unempl_2]=find(status1c_l1==0 & wage1c_l12>pi_wage_all_l  & wage1c_l12<pi_wage_all_h );
        
%         [I_ue_2]=find(status1c_l1==0 & status1c ==1 & death1c_l1 ==0 & death1c ==0 & wage1c_l12>pi_wage_all_l & wage1c_l12<pi_wage_all_h );
        
%         ue_rate_by_pctile(j)=length(I_ue_2)./length(I_unempl_2);
        
%     end
    
    
%     clearvars I_emp  I_emp_young age1c_long status1c_long
    
%     %Monthly moment
%     mean_wage_jf=mean(wage1c(I_ue)) ; %Avg wage of job finders
    
%     %Annual moment
%     death1c_annual_mo= death1c 	+ lag(death1c,1) + lag(death1c,2) + lag(death1c,3)+ lag(death1c,4)+ lag(death1c,5)+ lag(death1c,6) ...
%         + lag(death1c,7)+lag(death1c,8)+lag(death1c,9)+lag(death1c,10)+lag(death1c,11);
%     wage1c_annual_mo= wage1c + lag(wage1c,1) + lag(wage1c,2) + lag(wage1c,3)+ lag(wage1c,4)+ lag(wage1c,5)+ lag(wage1c,6) ...
%         + lag(wage1c,7)+lag(wage1c,8)+lag(wage1c,9)+lag(wage1c,10)+lag(wage1c,11);
    
%     [I_ue_ann]=find(status1c_l1==0 & status1c ==1 & death1c_annual_mo==0  & age1c>=3*12  & age1c<45*12);
    
%     mean_wage_jf_a =mean(wage1c_annual_mo(I_ue_ann)) ; %Avg wage of job finders, annual
    
%     clearvars wage1c_annual_mo death1c_annual_mo
    
%     %Employed, do not die, switch employers in between %CPS
%     [I_switch]=find( status1c_l12==1 & wage1c_l12>0 ...
%         &	 death1c_l12ma==0   & death1c ==0  & status1c==1  & status1c_l1==1   & status1c_l2==1  & tenure1c<12  & age1c>=3*12   & age1c<45*12); %EE switcher
%     frac_wage_drop_switch=mean((wage1c(I_switch)./wage1c_l12(I_switch)<1)) ;
    
%     %Switch employers, stay continuously employed
%     status1c_l12ma=(1/12)*(lag(status1c,1) + lag(status1c,2) + lag(status1c,3)+ lag(status1c,4)+ lag(status1c,5)+ lag(status1c,6) ...
%         + lag(status1c,7)+lag(status1c,8)+lag(status1c,9)+lag(status1c,10)+lag(status1c,11)+lag(status1c,12)	);
%     [I_switch_cont]=find(status1c_l12ma==1 ...
%         &	 death1c_l12ma==0 &  death1c ==0  & status1c==1 & tenure1c<12  & age1c>=3*12  & age1c<45*12 ); %EE switcher
%     frac_wage_drop_switch_cont=mean((wage1c(I_switch_cont)./wage1c_l12(I_switch_cont)<1)) ;
    
%     clearvars  I_switch I_switch_cont
    
%     %Wage Age Elasticity
%     [I_w_age]=find(status1c==1 &	wage1c>0   & age1c>=3*12   & age1c<45*12 ); %EE switcher
    
%     X=[ones(length(I_w_age),1),log( age1c(I_w_age)) ];
%     Y=winsorize(log(wage1c(I_w_age)+lp1),wins_lb_cut,wins_ub_cut) ;
    
%     [~,nX]=size(X);
%     B_w_age=(X'*X)\X'*Y;  %For w_{i,t}>=w_{-i,t}, regress w_{i,t+2} = b0+b1*w_{i,t}+b2*w_{-i,t}+eps
%     eps=Y-X*B_w_age;
%     sig_sq=sum(eps.^2)/length(I_w_age);
%     std_error_w_age=(sig_sq*eye(nX)/(X'*X))^.5;
    
%     clearvars I_w_age X Y eps
    
%     [I_w_age_y2]=find(status1c==1 &	wage1c>0 & age1c>=3*12 & age1c<4*12   ); % 24
%     [I_w_age_o2]=find(status1c==1 &	wage1c>0 & age1c>=33*12 & age1c<34*12   ); % 54
    
%     mean_wage_y2=mean(wage1c(I_w_age_y2));
%     mean_wage_o2=mean(wage1c(I_w_age_o2));
    
%     wage_ratio_old_young2=mean_wage_o2/mean_wage_y2;
    
%     clearvars I_w_age_y2  I_w_age_o2
    
%     %%%%%%%%%%%%%%%%%
%     %Wage growth--CPS
%     %%%%%%%%%%%%%%%%%
    
%     [I_switch_g]=find(status1c_l12==1 ...
%         &	 death1c_l12ma==0 &  death1c ==0  & status1c==1 & status1c_l1==1 & status1c_l2==1 & status1c_l3==1 &  tenure1c<12 & wage1c>0 & wage1c_l12>0  & age1c>=3*12  & age1c<45*12 ); %Switcher, growth
%     [I_stay_g]=find(status1c_l12==1 ...
%         &	 death1c_l12ma==0 &  death1c ==0  & status1c==1 &  tenure1c>=12 & wage1c>0 & wage1c_l12>0  & age1c>=3*12  & age1c<45*12); %Stayer, growth
    
%     ln_wage_growth_wins1_switch=winsorize(log(wage1c(I_switch_g) +lp1 )-log(wage1c_l12(I_switch_g)+lp1 ),wins_lb_cut,wins_ub_cut);
    
%     ln_wage_growth_wins1_stay=winsorize(log(wage1c(I_stay_g)+lp1 )-log(wage1c_l12(I_stay_g)+lp1  ),wins_lb_cut,wins_ub_cut);
    
%     wage_growth_switchers=mean(ln_wage_growth_wins1_switch );
    
%     wage_growth_stayer=mean(ln_wage_growth_wins1_stay) ;
    
%     frac_switch=length(I_switch_g)/(length(I_switch_g)+length(I_stay_g));
    
%     wage_growth=frac_switch*wage_growth_switchers+(1-frac_switch)*wage_growth_stayer;
    
%     frac_wage_growth_switch=frac_switch*wage_growth_switchers/wage_growth;
    
%     clearvars I_switch_g I_stay_g ln_wage_growth_wins1_switch ln_wage_growth_wins1_stay
    
%     [I_pos]=find( wage1c>0 & wage1c_l12>0  & age1c>=3*12  & age1c<45*12); %Stayer, growth
    
%     ln_wage_growth_wins1_pos=winsorize(log(wage1c(I_pos)+lp1)-log(wage1c_l12(I_pos)+lp1),wins_lb_cut,wins_ub_cut);
    
%     mean_ln_wage_change=mean(ln_wage_growth_wins1_pos);
    
%     std_dev_ln_wage_change=std(ln_wage_growth_wins1_pos);
    
%     clearvars  I_pos ln_wage_growth_wins1_pos
    
%     %JF Hazard by duration
%     for i=1:12
%         [I_unempl_1]=find(status1c_l1==0 & dur1c_l1==i & death1c_l1 ==0 & death1c ==0 & age1c>=3*12  & age1c<45*12);
        
%         [I_ue_1]=find(status1c_l1==0 & dur1c_l1==i & status1c ==1 & death1c_l1 ==0 & death1c ==0  & age1c>=3*12  & age1c<45*12);
        
%         ue_rate_by_dur(i)=length(I_ue_1)./length(I_unempl_1);
        
%         n_by_dur(i)=length(I_unempl_1);
%     end
    
%     clearvars  I_unempl_1 I_ue_1
    
%     %Total output
%     ttl_output=sum(udistprod.*b)+sum(solodistprod.*fsolo)+sum(sum(teamdistprod.*fteam))/2;
     
%     I_team_long=find(wage1c_long>0 & cowage1c_long>0); %All operating teams
    
%     economy_wide_avg=mean(log(wage1c_long(I_team_long)+lp1)); %economy wide avg
%     team_avg= (1/2)*(log(wage1c_long(I_team_long)+lp1)+log(cowage1c_long(I_team_long)+lp1)); %team by team avg
%     indiv_wage=log(wage1c_long(I_team_long)+lp1);
    
%     var_within=sum((indiv_wage-team_avg).^2)/length(I_team_long); %within
%     var_between=sum((team_avg-economy_wide_avg).^2)/length(I_team_long); %between
%     tot_var=sum((indiv_wage-economy_wide_avg).^2)/length(I_team_long); %Total variance
    
%     frac_between_tot=var_between/tot_var; %Share of total varaince between

%     clearvars  I_team_long indiv_wage
    
%     %Check monotonicity of wage among UE transitioners
%     for k=1:tpts
%         [I_samp]=find( status1c==1 & lag(status1c,1)==0 ...
%             &  death1c==0 & lag(death1c,1)==0 & type1c==k  );
        
%         ue_wage(k)=mean(wage1c(I_samp));
%     end
    
%     %Fraction with negative wages
%     neg_wage_share=sum(sum(wage1c<0))/sum(sum(status1c));
    
%     toc
    
%     clearvars I_empl_2  I_j_to_j I_unempl  I_ue  I_haz  I_unempl_2  I_ue_2  I_etype  I_ee_2  I_eu_2  I_ij_etype
    
%     %Store Stats
%     STAT_store(jjj,:)  = [ue_rate
%         ee_rate
%         eu_rate
%         wage_p90p10_all
%         wage_p90p10_young
%         mean_wage_young/mean_wage
%         mean_wage_jf/mean_wage
%         wage_ratio_old_young2
%         frac_wage_drop_switch_cont
%         frac_between_tot
%         B_own_wage_below(3)
%         B_own_wage_above(3)
%         ue_rate_by_dur(3)/ue_rate_by_dur(1)
%         B_own_wage_below(2)
%         B_own_wage_above(2)
%         B_co_wage(2)
%         B_co_wage(3)
%         rep_rate
%         B_j_to_j_below(2)
%         B_j_to_j(2)
%         B_j_to_j_above(2)
%         std_dev_ln_wage_change
%         frac_wage_growth_switch
%         wage_drop
%         mean_ln_wage_change
%         wage_growth_switchers
%         B_haz
%         1
%         1
%         1
%         1
%         ue_wage(1)
%         ue_wage(7)
%         B_co_wage_below(3)
%         B_co_wage_above(3)
%         B_co_wage_below(2)
%         B_co_wage_above(2)
%         B_eu(2)
%         B_eu_below(2)
%         B_eu_above(2)
%         I_corr_eu_tenure
%         I_corr_eu_hc
%         wage_growth_switchers
%         frac_wage_drop_switch
%         frac_wage_drop_switch_cont
%         neg_wage_share
%         wage_p90p10_all_a
%         wage_p90p10_young_a
%         mean_wage_jf_a/mean_wage_a
%         mean_wage_young_a/mean_wage_a
%         wage_ratio_old_young_a
%         mean_ln_wage_change_a
%         std_dev_ln_wage_change_a
%         mean_log_wage_a
%         var_log_wage_a
%         mean_log_wage_y_a
%         mean_log_wage_o_a
%         var_log_wage_y_a
%         var_log_wage_o_a
%         1
%         1
%         1
%         1
%         1
%         1
%         1
%         eue_mean_log_wage_growth_l1_f1   ]';
    
% end


% %% MOMENTS OUTPUT
% STAT=mean(STAT_store,1);%

% save main_results 








 

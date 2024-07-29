function [eplus_udist,eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist]=eplus_dist(ats,tpts,e3_udist,e3_edist,e3_mdist,e3_ndist,e3_tdist,death,typebirth)
    %eplus Distributions
    
    %Unemployed
    eplus_udist=(1-death)*e3_udist + death*typebirth;
    
    %Empty firms
    eplus_edist=e3_edist+ death*sum(e3_mdist,2)' + death*sum(e3_ndist,2)' + (death^2)*sum(e3_tdist,[2,3])';
    
    %Firms with manager
    eplus_mdist= (1-death)*e3_mdist + death*(1-death)*sum(e3_tdist,3);
    
    %Firms with no manager
    eplus_ndist= (1-death)*e3_ndist + death*(1-death)*reshape(sum(e3_tdist,2),ats,tpts);
            
    %Firms with team
    eplus_tdist= ((1-death)^2)*e3_tdist;
    
    %Check sum
    sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all"); %This is summing up to one if the fed in dist sums to 1
    
    

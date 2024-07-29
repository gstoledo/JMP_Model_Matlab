%Matrix of outside options for firm with one worker
function [u]=out_opt_lose1(Vmh,Veh,tpts,ats)
    d= reshape(permute(Vmh, [2 1 3]), size(Vmh, 2), [])';
    d=repmat(d,1,1,ats);
    u=  d- repmat(Veh',1,tpts,ats);
end

function [u,v,p] = MAC(u,v,p,f,g,D,m1,m2,L,N)
if N==L  % coarsest level: solve by DGS iterations
    for i = 1:100
        [u,v,p] = DGS(u,v,p,f,g,D,N);
    end
else
    % Presmoothing
    for k = 1:m1
        [u,v,p] = DGS(u,v,p,f,g,D,N);
    end
    % Form residual
    [ru,rv,rp]=Residual(u,v,p,f,g,D,N);
    % Restriction
    [ruc,rvc,rpc] = Res(ru,rv,rp,N);
    % Coarse grid correction
    [u1,v1,p1] = form(N/2);
    [u1,v1,p1] = MAC(u1,v1,p1,ruc,rvc,rpc,m1,m2,L,N/2);
    % Prolongation
    [eu,ev,ep] = Pro(u1,v1,p1,N);
    u=u+eu;v=v+ev;p=p+ep;
    % Postsmoothing
    for i=1:m2
        [u,v,p]=DGS(u,v,p,f,g,D,N);
    end
end
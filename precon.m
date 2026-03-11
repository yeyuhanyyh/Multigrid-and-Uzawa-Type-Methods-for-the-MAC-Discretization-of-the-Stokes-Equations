function [u,v]=precon(u,v,f,g,m1,m2,L,N)
if N==L  % coarsest level
    for i = 1:100
        [u,v] = SGS(u,v,f,g,N);
    end
else
    % Presmoothing
    for k = 1:m1
        [u,v] = SGS(u,v,f,g,N);
    end
    % Form residual
    [ru,rv]=Residual3(u,v,f,g,N);
    % Restriction
    [ruc,rvc] = Res2(ru,rv,N);
    % Coarse grid correction
    [u1,v1] = form(N/2);
    [u1,v1] = precon(u1,v1,ruc,rvc,m1,m2,L,N/2);
    % Prolongation
    [eu,ev] = Pro2(u1,v1,N);
    u=u+eu;v=v+ev;
    % Postsmoothing
    for i=1:m2
        [u,v]=SGS(u,v,f,g,N);
    end
end
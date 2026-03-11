function [u,v,p] = uzawa(u,v,p,f,g,N,e,a)
[u,v] = uzawa1(u,v,p,f,g,N,e);
p=uzawa2(u,v,p,N,a);
end
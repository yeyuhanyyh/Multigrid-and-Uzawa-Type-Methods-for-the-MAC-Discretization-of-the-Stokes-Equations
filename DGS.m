function [u,v,p] = DGS(u,v,p,f,g,D,n)
%DGS iterations process
[u,v] = DGS1(u,v,p,f,g,n);
[u,v,p] = DGS2(u,v,p,D,n);
end

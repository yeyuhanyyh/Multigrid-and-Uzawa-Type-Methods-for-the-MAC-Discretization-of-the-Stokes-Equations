function [u,v,p,b2] = inuzawa(u,v,p,f,g,a,m1,m2,L,N,e,b2)
[f1,g1]=setpre(p,f,g,N);
[u,v] = inuzawa1(u,v,f1,g1,m1,m2,L,N,e,b2);
[p,b2]=uzawa2(u,v,p,N,a);
end
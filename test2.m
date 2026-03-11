N=256;
e=1e-100;
a=1.5;
cre;
t=0;
[ru,rv,rp] = Residual(u,v,p,f,g,D,N);
re=recal(ru,rv,rp,N);
er=ercal(u,v,uu,vv,N);
r0=re;
ite=0;
while re>1e-8*r0
    tic
    [u,v,p] = uzawa(u,v,p,f,g,N,e,1);
    t=t+toc;
    [ru,rv,rp] = Residual(u,v,p,f,g,D,N);
    re=recal(ru,rv,rp,N);
    er=ercal(u,v,uu,vv,N);
    ite=ite+1;
end
disp('N:');disp(N);
disp('time(s):');disp(t);
disp('iteration:');disp(ite);
disp('residual:');disp(re);
disp('error:');disp(er);
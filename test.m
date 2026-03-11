N=1024;
L1=4;
m1=2;m2=2;
L=N/L1;
cre;
t=0;
[ru,rv,rp] = Residual(u,v,p,f,g,D,N);
re=recal(ru,rv,rp,N);
er=ercal(u,v,uu,vv,N);
r0=re;
ite=0;
while re>1e-8*r0
    tic
    [u,v,p] = MAC(u,v,p,f,g,D,m1,m2,L1,N);
    t=t+toc;
    [ru,rv,rp] = Residual(u,v,p,f,g,D,N);
    re=recal(ru,rv,rp,N);
    er=ercal(u,v,uu,vv,N);
    ite=ite+1;
end
disp('N:');disp(N);
disp('L:');disp(L);
disp('v1:');disp(m1);
disp('v2:');disp(m2);
disp('time(s):');disp(t);
disp('iteration:');disp(ite);
disp('residual:');disp(re);
disp('error:');disp(er);


N=1024;
e=1e-20;
a=1;
m1=2;
m2=2;
L=2;
cre;
t=0;
[ru,rv,rp] = Residual(u,v,p,f,g,D,N);
re=recal(ru,rv,rp,N);
er=ercal(u,v,uu,vv,N);
r0=re;
ite=0;
b2=-1;
while re>1e-8*r0
    tic
    [u,v,p,b2] = inuzawa(u,v,p,f,g,a,m1,m2,L,N,e,b2);
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
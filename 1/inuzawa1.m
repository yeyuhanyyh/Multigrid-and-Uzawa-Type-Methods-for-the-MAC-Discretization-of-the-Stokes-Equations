function [u,v] = inuzawa1(u,v,f,g,m1,m2,L,N,e,b2)
if b2==-1
k=0;
[ru,rv] = Residual3(u,v,f,g,N);
while k<2
    u0=zeros(N+1,N);v0=zeros(N,N+1);
    [zu,zv]=precon(u0,v0,ru,rv,m1,m2,L,N);
    k=k+1;
    if k==1
        pu=zu;pv=zv;rou=product2(ru,rv,zu,zv,N);
    else
        rou1=rou;rou=product2(ru,rv,zu,zv,N);
        b=rou/rou1;
        pu=zu+b*pu;pv=zv+b*pv;
    end
    [wu,wv]=productAp(pu,pv,N);
    a=rou/product2(pu,pv,wu,wv,N);
    u=u+a*pu; v=v+a*pv;
    ru=ru-a*wu; rv=rv-a*wv;
end
disp('ite(PCG):');disp(k);
else
k=0;
[ru,rv] = Residual3(u,v,f,g,N);
rou=b2;
while k<N*N&&rou>e*e*b2
    u0=zeros(N+1,N);v0=zeros(N,N+1);
    [zu,zv]=precon(u0,v0,ru,rv,m1,m2,L,N);
    k=k+1;
    if k==1
        pu=zu;pv=zv;rou=product2(ru,rv,zu,zv,N);
    else
        rou1=rou;rou=product2(ru,rv,zu,zv,N);
        b=rou/rou1;
        pu=zu+b*pu;pv=zv+b*pv;
    end
    [wu,wv]=productAp(pu,pv,N);
    a=rou/product2(pu,pv,wu,wv,N);
    u=u+a*pu; v=v+a*pv;
    ru=ru-a*wu; rv=rv-a*wv;
end
disp('ite(PCG):');disp(k);
end
end
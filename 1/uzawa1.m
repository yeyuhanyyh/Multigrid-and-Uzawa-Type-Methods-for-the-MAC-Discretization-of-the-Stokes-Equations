function [u,v] = uzawa1(u,v,p,f,g,N,e)
k=0;
[ru,rv] = Residual2(u,v,p,f,g,N);
rou=product(ru,rv,N);
u1=zeros(N+1,N);v1=zeros(N,N+1);
[bu,bv] = Residual2(u1,v1,p,f,g,N);
b2=product(bu,bv,N);
while k<N*N&&rou>e*e*b2
    k=k+1;
    if k==1
        pu=ru;pv=rv;
    else
        b=rou/rou1;
        pu=ru+b*pu;pv=rv+b*pv;
    end
    [wu,wv]=productAp(pu,pv,N);
    a=rou/product2(pu,pv,wu,wv,N);
    u=u+a*pu; v=v+a*pv;
    ru=ru-a*wu; rv=rv-a*wv;rou1=rou;rou=product(ru,rv,N); 
end
disp('ite(CG):');disp(k);
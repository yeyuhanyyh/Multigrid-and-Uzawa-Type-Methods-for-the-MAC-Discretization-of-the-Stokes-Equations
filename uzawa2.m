function [p,b2] = uzawa2(u,v,p,N,a)
h=1/N;
b2=0;
for i=1:N
    for j=1:N
        d=(u(i+1,j)-u(i,j)+v(i,j+1)-v(i,j))/h;
        p(i,j)=p(i,j)-d*a;
        b2=b2+d*d;
    end
end
end
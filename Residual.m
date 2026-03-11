function [ru,rv,rp] = Residual(u,v,p,f,g,D,N)
%calculating the residual for Q1
h=1/N;

for j=1:N
    ru(1,j)=0;
end
for i=2:N
    j=1;
    f1=(3*u(i,j)-u(i-1,j)-u(i+1,j)-u(i,j+1)+h*(p(i,j)-p(i-1,j)))/(h*h);
    ru(i,j)=f(i,j)-f1;
  
for j=2:N-1
     f1=(4*u(i,j)-u(i-1,j)-u(i+1,j)-u(i,j-1)-u(i,j+1)+h*(p(i,j)-p(i-1,j)))/(h*h);
     ru(i,j)=f(i,j)-f1;
end

    j=N;
    f1=(3*u(i,j)-u(i-1,j)-u(i+1,j)-u(i,j-1)+h*(p(i,j)-p(i-1,j)))/(h*h);
     ru(i,j)=f(i,j)-f1;
end
for j=1:N
    ru(N+1,j)=0;
end

for i=1:N
    rv(i,1)=0;
end
for j=2:N
    i=1;
    g1=(3*v(i,j)-v(i+1,j)-v(i,j-1)-v(i,j+1)+h*(p(i,j)-p(i,j-1)))/(h*h);
    rv(i,j)=g(i,j)-g1;   
    
    for i=2:N-1
    g1=(4*v(i,j)-v(i-1,j)-v(i+1,j)-v(i,j-1)-v(i,j+1)+h*(p(i,j)-p(i,j-1)))/(h*h);
    rv(i,j)=g(i,j)-g1;
    end

    i=N;
    g1=(3*v(i,j)-v(i-1,j)-v(i,j-1)-v(i,j+1)+h*(p(i,j)-p(i,j-1)))/(h*h);
    rv(i,j)=g(i,j)-g1;
end

for i=1:N
    rv(i,N+1)=0;
end

for i=1:N
    for j=1:N
        rp(i,j)=D(i,j)-(u(i,j)-u(i+1,j)+v(i,j)-v(i,j+1))/h;
    end
end
end


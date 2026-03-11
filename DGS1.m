function [u,v] = DGS1(u,v,p,f,g,n)
%first step of DGS
h=1/n;
for j=1:n
    u(1,j)=0;
end
for i=2:n
    j=1;
    u(i,j)=(h*h*f(i,j)+u(i-1,j)+u(i+1,j)+u(i,j+1)...
    -h*(p(i,j)-p(i-1,j)))/3;
  
for j=2:n-1
     u(i,j)=(h*h*f(i,j)+u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)...
     -h*(p(i,j)-p(i-1,j)))/4;
end

    j=n;
    u(i,j)=(h*h*f(i,j)+u(i-1,j)+u(i+1,j)+u(i,j-1)...
    -h*(p(i,j)-p(i-1,j)))/3;
end
for j=1:n
    u(n+1,j)=0;
end

for i=1:n
    v(i,1)=0;
end
for j=2:n
    i=1;
    v(i,j)=(h*h*g(i,j)+v(i+1,j)+v(i,j-1)+v(i,j+1)...
    -h*(p(i,j)-p(i,j-1)))/3;   

    for i=2:n-1
        v(i,j)=(h*h*g(i,j)+v(i-1,j)+v(i+1,j)+v(i,j-1)+v(i,j+1)...
        -h*(p(i,j)-p(i,j-1)))/4;
    end
    i=n;
    v(i,j)=(h*h*g(i,j)+v(i-1,j)+v(i,j-1)+v(i,j+1)...
    -h*(p(i,j)-p(i,j-1)))/3;
end
for i=1:n
    v(i,n+1)=0;
end
end

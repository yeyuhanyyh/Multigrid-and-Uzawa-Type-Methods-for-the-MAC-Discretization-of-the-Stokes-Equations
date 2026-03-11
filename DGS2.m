function [u,v,p] = DGS2(u,v,p,D,n)
% second step of DGS 
h=1/n;
for i=2:n-1
    for j=2:n-1
        r=(u(i,j)-u(i+1,j)+v(i,j)-v(i,j+1))/h-D(i,j);
        d=r*h/4;
        u(i,j)=u(i,j)-d;
        u(i+1,j)=u(i+1,j)+d;
        v(i,j)=v(i,j)-d;
        v(i,j+1)=v(i,j+1)+d;
        p(i,j)=p(i,j)+r;
        p(i+1,j)=p(i+1,j)-r/4;
        p(i-1,j)=p(i-1,j)-r/4;
        p(i,j+1)=p(i,j+1)-r/4;
        p(i,j-1)=p(i,j-1)-r/4;
    end
end
j=n;
for i=2:n-1
    r=(u(i,j)-u(i+1,j)+v(i,j)-v(i,j+1))/h-D(i,j);
    d=r*h/3;
    u(i,j)=u(i,j)-d;
    u(i+1,j)=u(i+1,j)+d;
    v(i,j)=v(i,j)-d;
    p(i,j)=p(i,j)+r;
    p(i+1,j)=p(i+1,j)-r/3;
    p(i-1,j)=p(i-1,j)-r/3;
    p(i,j-1)=p(i,j-1)-r/3;
end

i=n;
for j=2:n-1
    r=(u(i,j)-u(i+1,j)+v(i,j)-v(i,j+1))/h-D(i,j);
    d=r*h/3;
    u(i,j)=u(i,j)-d;
    v(i,j)=v(i,j)-d;
    v(i,j+1)=v(i,j+1)+d;
    p(i,j)=p(i,j)+r;
    p(i-1,j)=p(i-1,j)-r/3;
    p(i,j+1)=p(i,j+1)-r/3;
    p(i,j-1)=p(i,j-1)-r/3;
end

i=1;
for j=2:n-1
    r=(u(i,j)-u(i+1,j)+v(i,j)-v(i,j+1))/h-D(i,j);
    d=r*h/3;
    u(i+1,j)=u(i+1,j)+d;
    v(i,j)=v(i,j)-d;
    v(i,j+1)=v(i,j+1)+d;
    p(i,j)=p(i,j)+r;
    p(i+1,j)=p(i+1,j)-r/3;
    p(i,j+1)=p(i,j+1)-r/3;
    p(i,j-1)=p(i,j-1)-r/3;
end




j=1;
for i=2:n-1
    r=(u(i,j)-u(i+1,j)+v(i,j)-v(i,j+1))/h-D(i,j);
    d=r*h/3;
    u(i,j)=u(i,j)-d;
    u(i+1,j)=u(i+1,j)+d;
    v(i,j+1)=v(i,j+1)+d;
    p(i,j)=p(i,j)+r;
    p(i+1,j)=p(i+1,j)-r/3;
    p(i-1,j)=p(i-1,j)-r/3;
    p(i,j+1)=p(i,j+1)-r/3;
end

i=1; j=1;
r=(u(i,j)-u(i+1,j)+v(i,j)-v(i,j+1))/h-D(i,j);
d=r*h/2;
u(i+1,j)=u(i+1,j)+d;
v(i,j+1)=v(i,j+1)+d;
p(i,j)=p(i,j)+r;

p(i,j+1)=p(i,j+1)-r/2;
p(i+1,j)=p(i+1,j)-r/2;

i=1; j=n;
r=(u(i,j)-u(i+1,j)+v(i,j)-v(i,j+1))/h-D(i,j);
d=r*h/2;
u(i+1,j)=u(i+1,j)+d;
v(i,j)=v(i,j)-d;
p(i,j)=p(i,j)+r;
p(i+1,j)=p(i+1,j)-r/2;
p(i,j-1)=p(i,j-1)-r/2;

i=n; j=1;
r=(u(i,j)-u(i+1,j)+v(i,j)-v(i,j+1))/h-D(i,j);
d=r*h/2;
u(i,j)=u(i,j)-d;
v(i,j+1)=v(i,j+1)+d;
p(i,j)=p(i,j)+r;
p(i-1,j)=p(i-1,j)-r/2; 
p(i,j+1)=p(i,j+1)-r/2;

i=n; j=n;
r=(u(i,j)-u(i+1,j)+v(i,j)-v(i,j+1))/h-D(i,j);
d=r*h/2;
u(i,j)=u(i,j)-d;
v(i,j)=v(i,j)-d;
p(i,j)=p(i,j)+r;
p(i-1,j)=p(i-1,j)-r/2; 
p(i,j-1)=p(i,j-1)-r/2;


for i=1:n
    v(i,n+1)=0;
end

for j=1:n
    u(n+1,j)=0;
end

end



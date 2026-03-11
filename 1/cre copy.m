%set up basical values
h=1/N;
p=zeros(N,N);
u=zeros(N+1,N);
v=zeros(N,N+1);
D=zeros(N,N);
for i=2:N
    x=(i-1)*h;
    b(i)=-2*pi*(1-cos(2*pi*x));
    t(i)=2*pi*(1-cos(2*pi*x));
end
for j=2:N
    y=(j-1)*h;
    l(j)=2*pi*(1-cos(2*pi*y)); 
    r(j)=-2*pi*(1-cos(2*pi*y));
end

for i=1:N+1
    for j=1:N
        x=(i-1)*h;
        y=(j-0.5)*h;
        f(i,j)=-4*pi*pi*(2*cos(2*pi*x)-1)*sin(2*pi*y)+x*x;
    end
end

for i=1:N
    for j=1:N+1
        x=(i-0.5)*h;
        y=(j-1)*h;
        g(i,j)=4*pi*pi*(2*cos(2*pi*y)-1)*sin(2*pi*x);
    end
end

for i=2:N
    f(i,1)=f(i,1)+b(i)/h;
    f(i,N)=f(i,N)+t(i)/h;
end
for j=2:N
    g(1,j)=g(1,j)+l(j)/h;
    g(N,j)=g(N,j)+r(j)/h;
end

for i=1:N+1
    for j=1:N
        x=(i-1)*h;
        y=(j-0.5)*h;
        uu(i,j)=(1-cos(2*pi*x))*sin(2*pi*y);
    end
end
for i=1:N
    for j=1:N+1
        x=(i-0.5)*h;
        y=(j-1)*h;
        vv(i,j)=-(1-cos(2*pi*y))*sin(2*pi*x);
    end
end
for i=1:N
    for j=1:N
        x=(i-0.5)*h;
        y=(j-0.5)*h;
        pp(i,j)=x*x*x/3-1/12;
    end
end

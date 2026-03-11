function [f1,g1]=setpre(p,f,g,N)
h=1/N;
for j=1:N
    f1(1,j)=0;
end

for i=2:N
    for j=1:N
    f1(i,j)=f(i,j)-(p(i,j)-p(i-1,j))/h;
    end
end

for j=1:N
    f1(N+1,j)=0;
end

for i=1:N
    g1(i,1)=0;
end

for j=2:N
    for i=1:N
    g1(i,j)=g(i,j)-(p(i,j)-p(i,j-1))/h;
    end
end

for i=1:N
    g1(i,N+1)=0;
end

end
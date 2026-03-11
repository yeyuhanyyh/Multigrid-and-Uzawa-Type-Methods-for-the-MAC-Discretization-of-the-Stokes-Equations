function [Au,Av]= productAp(pu,pv,N)
%calculate Pt*A*P
ans=0;
h=1/N;
for j=1:N
    Au(1,j)=0;
end
for i=2:N
    j=1;
    Au(i,j)=(3*pu(i,j)-pu(i-1,j)-pu(i+1,j)-pu(i,j+1))/(h*h);
  
for j=2:N-1
     Au(i,j)=(4*pu(i,j)-pu(i-1,j)-pu(i+1,j)-pu(i,j+1)-pu(i,j-1))/(h*h);
end
    j=N;
    Au(i,j)=(3*pu(i,j)-pu(i-1,j)-pu(i+1,j)-pu(i,j-1))/(h*h);

end
for j=1:N
    Au(N+1,j)=0;
end

for i=1:N
    Av(i,1)=0;
end
for j=2:N
    i=1;
    Av(i,j)=(3*pv(i,j)-pv(i,j-1)-pv(i+1,j)-pv(i,j+1))/(h*h);

    for i=2:N-1
        Av(i,j)=(4*pv(i,j)-pv(i-1,j)-pv(i,j-1)-pv(i+1,j)-pv(i,j+1))/(h*h);
    end
    i=N;
     Av(i,j)=(3*pv(i,j)-pv(i-1,j)-pv(i,j-1)-pv(i,j+1))/(h*h);
end
for i=1:N
    Av(i,N+1)=0;
end
end
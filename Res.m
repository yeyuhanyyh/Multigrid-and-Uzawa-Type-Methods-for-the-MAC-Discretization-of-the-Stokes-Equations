function [ruc,rvc,rpc] = Res(ru,rv,rp,N)
%restrict from N*N to (N/2)*(N/2)
for j=1:N/2
    ruc(1,j)=0;
    ruc(N/2+1,j)=0;
end
for i=2:N/2
    for j=1:N/2
        ruc(i,j)=(ru(2*i,2*j)+ru(2*i,2*j-1)+ru(2*i-2,2*j)+ru(2*i-2,2*j-1))/8 ...
        +(ru(2*i-1,2*j)+ru(2*i-1,2*j-1))/4;
    end
end

for i=1:N/2
    rvc(i,1)=0;
    rvc(i,N/2+1)=0;
end
for j=2:N/2
    for i=1:N/2
        rvc(i,j)=(rv(2*i,2*j)+rv(2*i-1,2*j)+rv(2*i,2*j-2)+rv(2*i-1,2*j-2))/8 ...
        +(rv(2*i,2*j-1)+rv(2*i-1,2*j-1))/4;
    end
end

for i=1:N/2
    for j=1:N/2
        rpc(i,j)=(rp(2*i,2*j)+rp(2*i-1,2*j)+rp(2*i,2*j-1)+rp(2*i-1,2*j-1))/4;
    end
end
end
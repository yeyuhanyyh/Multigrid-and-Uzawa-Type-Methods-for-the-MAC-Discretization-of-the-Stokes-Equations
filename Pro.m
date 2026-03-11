function [eu,ev,ep] = Pro(euc,evc,epc,N)
%prolong from (N/2)*(N/2) to N*N 
for i=1:N/2
    for j=1:N/2
        ep(2*i,2*j)=epc(i,j);
        ep(2*i-1,2*j)=epc(i,j);
        ep(2*i,2*j-1)=epc(i,j);
        ep(2*i-1,2*j-1)=epc(i,j);
    end
end

for j=1:N/2
    for i=1:N/2
        eu(2*i-1,2*j-1)=euc(i,j);
        eu(2*i-1,2*j)=euc(i,j);
        eu(2*i,2*j-1)=(euc(i,j)+euc(i+1,j))/2;
        eu(2*i,2*j)=(euc(i,j)+euc(i+1,j))/2;
    end
end

for j=1:N
    eu(1,j)=0;eu(N+1,j)=0;
end

for i=1:N/2
    for j=1:N/2
        ev(2*i-1,2*j-1)=evc(i,j);
        ev(2*i,2*j-1)=evc(i,j);
        ev(2*i-1,2*j)=(evc(i,j)+evc(i,j+1))/2;
        ev(2*i,2*j)=(evc(i,j)+evc(i,j+1))/2;
    end
end

for i=1:N
    ev(i,1)=0; ev(i,N+1)=0;
end
end


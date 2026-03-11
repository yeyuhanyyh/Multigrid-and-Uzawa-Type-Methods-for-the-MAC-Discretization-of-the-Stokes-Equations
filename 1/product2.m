function ans= product2(ru,rv,ru1,rv2,N)
%calculate Rt*R
ans=0;
for i=1:N+1
    for j=1:N
        ans=ans+ru(i,j)*ru1(i,j)+rv(j,i)*rv2(j,i);
    end
end
end
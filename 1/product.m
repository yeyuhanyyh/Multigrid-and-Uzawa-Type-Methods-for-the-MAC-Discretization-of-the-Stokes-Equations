function ans= product(ru,rv,N)
%calculate Rt*R
ans=0;
for i=1:N
    for j=1:N
        ans=ans+ru(i,j)*ru(i,j)+rv(i,j)*rv(i,j);
    end
end
end
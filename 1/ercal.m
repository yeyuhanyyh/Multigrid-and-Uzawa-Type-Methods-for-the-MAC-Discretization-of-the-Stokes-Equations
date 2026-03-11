function er=ercal(u,v,uu,vv,N)
%calculate the norm of error
    er=0;
    for i=1:N+1
        for j=1:N
        er=er+(u(i,j)-uu(i,j))*(u(i,j)-uu(i,j))+(v(j,i)-vv(j,i))*(v(j,i)-vv(j,i));
        end
    end
    er=sqrt(er)/N;
end


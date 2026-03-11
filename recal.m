function re=recal(ru,rv,rp,N)
%calculating the norm of residual
    re=0;
    for i=1:N+1
        for j=1:N
        re=re+ru(i,j)*ru(i,j)+rv(j,i)*rv(j,i);
        end
    end
    for i=1:N
        for j=1:N
            re=re+rp(i,j)*rp(i,j);
        end
    end
    re=sqrt(re)/N;
end
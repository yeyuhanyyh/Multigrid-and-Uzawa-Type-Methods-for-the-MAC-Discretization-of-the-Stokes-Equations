U=zeros(2*N*(N+1),1);
DA=diag(A);
LA=A-triu(A);
M=inv(DA-LA);
temp=0;
n=N*(N+1)
while temp<100
U=U+M*(F-B*P-A*U);
for i=2:N-1
    for j=2:N-1
        R(i,j)=U((i-1)*N+j)-U(i*N+j)+U((i-1)*(N+1)+j+n)-U((i-1)*(N+1)+j+n+1);
        U((i-1)*N+j)=U((i-1)*N+j)-R(i,j)*h/4;
        U((i-1)*(N+1)+j+n)=U((i-1)*(N+1)+j+n)-R(i,j)*h/4;
        U(i*N+j)=U(i*N+j)+R(i,j)*h/4;
        U((i-1)*(N+1)+j+n+1)=U((i-1)*(N+1)+j+n+1)+R(i,j)*h/4;
        P(N*(i-1)+j)=P(N*(i-1)+j)+R(i,j);
        P(N*(i-2)+j)=P(N*(i-2)+j)-R(i,j)/4;
        P(N*i+j)=P(N*i+j)-R(i,j)/4;
        P(N*(i-1)+j-1)=P(N*(i-1)+j-1)-R(i,j)/4;
        P(N*(i-1)+j+1)=P(N*(i-1)+j+1)-R(i,j)/4;
    end
end

for i=2:N-1
    j=N;
    R(i,j)=U((i-1)*N+j)-U(i*N+j)+U((i-1)*(N+1)+j+n)-U((i-1)*(N+1)+j+n+1);
    U((i-1)*N+j)=U((i-1)*N+j)-R(i,j)*h/3;
    U((i-1)*(N+1)+j+n)=U((i-1)*(N+1)+j+n)-R(i,j)*h/3;
    U(i*N+j)=U(i*N+j)+R(i,j)*h/3;
    P(N*(i-1)+j)=P(N*(i-1)+j)+R(i,j);
    P(N*(i-2)+j)=P(N*(i-2)+j)-R(i,j)/3;
    P(N*i+j)=P(N*i+j)-R(i,j)/3;
    P(N*(i-1)+j-1)=P(N*(i-1)+j-1)-R(i,j)/3;
end

for i=2:N-1
    j=1;
    R(i,j)=U((i-1)*N+j)-U(i*N+j)+U((i-1)*(N+1)+j+n)-U((i-1)*(N+1)+j+n+1);
    U((i-1)*N+j)=U((i-1)*N+j)-R(i,j)*h/3;
    U((i-1)*(N+1)+j+n)=U((i-1)*(N+1)+j+n)-R(i,j)*h/3;
    U(i*N+j)=U(i*N+j)+R(i,j)*h/3;
    U((i-1)*(N+1)+j+n+1)=U((i-1)*(N+1)+j+n+1)+R(i,j)*h/4;
    P(N*(i-1)+j)=P(N*(i-1)+j)+R(i,j);
    P(N*(i-2)+j)=P(N*(i-2)+j)-R(i,j)/3;
    P(N*i+j)=P(N*i+j)-R(i,j)/3;
    P(N*(i-1)+j+1)=P(N*(i-1)+j+1)-R(i,j)/3;
end

for j=2:N-1
    i=N;
    R(i,j)=U((i-1)*N+j)-U(i*N+j)+U((i-1)*(N+1)+j+n)-U((i-1)*(N+1)+j+n+1);
    U((i-1)*N+j)=U((i-1)*N+j)-R(i,j)*h/3;
    U((i-1)*(N+1)+j+n)=U((i-1)*(N+1)+j+n)-R(i,j)*h/3;
    U((i-1)*(N+1)+j+n+1)=U((i-1)*(N+1)+j+n+1)+R(i,j)*h/3;
    P(N*(i-1)+j)=P(N*(i-1)+j)+R(i,j);
    P(N*(i-2)+j)=P(N*(i-2)+j)-R(i,j)/3;
    P(N*(i-1)+j-1)=P(N*(i-1)+j-1)-R(i,j)/3;
    P(N*(i-1)+j+1)=P(N*(i-1)+j+1)-R(i,j)/3;
end

for j=2:N-1
    i=1;
    R(i,j)=U((i-1)*N+j)-U(i*N+j)+U((i-1)*(N+1)+j+n)-U((i-1)*(N+1)+j+n+1);
    U((i-1)*(N+1)+j+n)=U((i-1)*(N+1)+j+n)-R(i,j)*h/3;
    U(i*N+j)=U(i*N+j)+R(i,j)*h/3;
    U((i-1)*(N+1)+j+n+1)=U((i-1)*(N+1)+j+n+1)+R(i,j)*h/3;
    P(N*(i-1)+j)=P(N*(i-1)+j)+R(i,j);
    P(N*i+j)=P(N*i+j)-R(i,j)/3;
    P(N*(i-1)+j-1)=P(N*(i-1)+j-1)-R(i,j)/3;
    P(N*(i-1)+j+1)=P(N*(i-1)+j+1)-R(i,j)/3;
end

i=1;j=1; R(i,j)=U((i-1)*N+j)-U(i*N+j)+U((i-1)*(N+1)+j+n)-U((i-1)*(N+1)+j+n+1);
 U(i*N+j)=U(i*N+j)+R(i,j)*h/2;
 U((i-1)*(N+1)+j+n+1)=U((i-1)*(N+1)+j+n+1)+R(i,j)*h/2;
 P(N*(i-1)+j)=P(N*(i-1)+j)+R(i,j);
 P(N*i+j)=P(N*i+j)-R(i,j)/2;
 P(N*(i-1)+j+1)=P(N*(i-1)+j+1)-R(i,j)/2;

 i=1;j=N;
 R(i,j)=U((i-1)*N+j)-U(i*N+j)+U((i-1)*(N+1)+j+n)-U((i-1)*(N+1)+j+n+1);
 U((i-1)*(N+1)+j+n)=U((i-1)*(N+1)+j+n)-R(i,j)*h/2;
 U(i*N+j)=U(i*N+j)+R(i,j)*h/2;
 P(N*(i-1)+j)=P(N*(i-1)+j)+R(i,j);
 P(N*i+j)=P(N*i+j)-R(i,j)/2;
 P(N*(i-1)+j-1)=P(N*(i-1)+j-1)-R(i,j)/2;

 i=N;j=1;
 R(i,j)=U((i-1)*N+j)-U(i*N+j)+U((i-1)*(N+1)+j+n)-U((i-1)*(N+1)+j+n+1);
 U((i-1)*N+j)=U((i-1)*N+j)-R(i,j)*h/2;
 U((i-1)*(N+1)+j+n+1)=U((i-1)*(N+1)+j+n+1)+R(i,j)*h/2;
 P(N*(i-1)+j)=P(N*(i-1)+j)+R(i,j);
 P(N*(i-2)+j)=P(N*(i-2)+j)-R(i,j)/2;
 P(N*(i-1)+j+1)=P(N*(i-1)+j+1)-R(i,j)/2;

 i=N;j=N;
 R(i,j)=U((i-1)*N+j)-U(i*N+j)+U((i-1)*(N+1)+j+n)-U((i-1)*(N+1)+j+n+1);
 U((i-1)*N+j)=U((i-1)*N+j)-R(i,j)*h/2;
 U((i-1)*(N+1)+j+n)=U((i-1)*(N+1)+j+n)-R(i,j)*h/2;
 P(N*(i-1)+j)=P(N*(i-1)+j)+R(i,j);
 P(N*(i-2)+j)=P(N*(i-2)+j)-R(i,j)/2;
 P(N*(i-1)+j-1)=P(N*(i-1)+j-1)-R(i,j)/2;

 temp=temp+1;
end



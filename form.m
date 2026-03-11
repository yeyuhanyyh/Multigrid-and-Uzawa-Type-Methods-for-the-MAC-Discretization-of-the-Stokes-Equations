function [u,v,p] = form(N)
%forming initial value of u,v,p
p=zeros(N,N);
u=zeros(N+1,N);
v=zeros(N,N+1);
end
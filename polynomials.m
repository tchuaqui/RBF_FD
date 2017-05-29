function [ pol,dpol,d2pol ] = polynomials( x_dados,n,p)
pol=zeros(5,n); 
pol(1,:)=ones(1,n); 
pol(2,:)=x_dados;
pol(3,:)=x_dados.^2;
pol(4,:)=x_dados.^3;
pol(5,:)=x_dados.^4;

dpol=zeros(5,n);
dpol(1,:)=zeros(1,n); 
dpol(2,:)=ones(1,n);
dpol(3,:)=2*x_dados;
dpol(4,:)=3*x_dados.^2;
dpol(5,:)=4*x_dados.^3;

d2pol=zeros(5,n);
d2pol(1,:)=zeros(1,n); 
d2pol(2,:)=zeros(1,n);
d2pol(3,:)=2*ones(1,n);
d2pol(4,:)=6*x_dados;
d2pol(5,:)=12*x_dados.^2;

pol=pol(1:p,:);
dpol=dpol(1:p,:);
d2pol=d2pol(1:p,:);
end



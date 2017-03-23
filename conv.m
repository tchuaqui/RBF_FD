function [ freq,freq_exacta ] = conv( n )
x_inicial=0;x_final=100e-3;
L=x_final-x_inicial;
cfstr='cl';
x_dados=[0:L/(n-1):L];dist=x_dados(3)-x_dados(1);
[xi,xj]=meshgrid(x_dados);
x_central=find(x_dados==0.5);
%%
k=5/6;
propmec=[5e-3 6e10 2.3e10 7500]; %[esp E G rho]
propel=[-16.492145 -16.492145 2.588549e-8]; %[e31 e32 ezz]
esp=[propmec(1) propmec(1)];
rho=propmec(4);

e31=propel(1);
ezz=propel(3);

Q11=[propmec(2) propmec(2)];
Q55=[propmec(3) propmec(3)];

% I=1*h^3/12;
% E=1/I;rho=1/A;
% G=E/(2*(1+0.25));
%%
ncamadas=2;
h=0;
for i=1:ncamadas
h=h+esp(i);    
end

z=zeros(ncamadas+1,1);
z(1)=-h/2;
zm=zeros(ncamadas,1);
for i=2:ncamadas+1
z(i)=z(i-1)+esp(i-1);
zm(i-1)=(z(i-1)+z(i))/2;
end

I0=zeros(ncamadas,1); I1=zeros(ncamadas,1); I2=zeros(ncamadas,1); J0=0; J1=0; J2=0;
B11=0; B55=0; C11=0; D11=0; 
for i=1:ncamadas    
I0(i)=z(i+1)-z(i);
I1(i)=(z(i+1)^2-z(i)^2)/2;
I2(i)=(z(i+1)^3-z(i)^3)/3;

J0=J0+rho*I0(i);
J1=J1+rho*I1(i);
J2=J2+rho*I2(i);

B11=B11+Q11(i)*I0(i);
C11=C11+Q11(i)*I1(i);
B55=B55+k*Q55(i)*I0(i);
% Atheta=Atheta+Q11(i)*I1(i);
% Btheta=Btheta-k*Q55(i)*I0(i);
D11=D11+Q11(i)*I2(i);   %c1theta
% C2theta=C2theta-k*Q55(i)*I0(i); 
% Mu=Mu+Q11(i)*I1(i);
% Mtheta=Mtheta+Q11(i)*I2(i);
end

F1=e31*I0(1)/esp(1);
F2=e31*I0(2)/esp(2);
H1=e31*I1(1)/esp(1);
H2=e31*I1(2)/esp(2);
I_1=ezz*I0(1)/(esp(1)^2);
I_2=ezz*I0(2)/(esp(2)^2);
%%
rhs_1_u=zeros(3,numel(x_dados));  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rhs_1_w=zeros(3,numel(x_dados));
rhs_1_theta=zeros(3,numel(x_dados));

rhs_2_u=zeros(3,numel(x_dados));    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rhs_2_w=zeros(3,numel(x_dados));
rhs_2_theta=zeros(3,numel(x_dados));

rhs_3_u=zeros(3,numel(x_dados));    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rhs_3_w=zeros(3,numel(x_dados));     %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rhs_3_theta=zeros(3,numel(x_dados));  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

arhs_1_u=zeros(3,numel(x_dados));         %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
arhs_1_w=zeros(3,numel(x_dados));
arhs_1_theta=zeros(3,numel(x_dados));

arhs_2_u=zeros(3,numel(x_dados));              %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
arhs_2_w=zeros(3,numel(x_dados));
arhs_2_theta=zeros(3,numel(x_dados));

arhs_3_u=zeros(3,numel(x_dados));              %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
arhs_3_w=zeros(3,numel(x_dados));                 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
arhs_3_theta=zeros(3,numel(x_dados));            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c=0.1;  %2*dist/sqrt(sqrt(3))   
for i=1:numel(x_dados)
if i==1
        sub_dominio=[x_dados(i), x_dados(i+1), x_dados(i+2)];
        switch cfstr
            case {'ss'}
        rhs_1_u(:,i)=g(c,x_dados(i),sub_dominio(:));        
        rhs_1_w(:,i)=0;
        rhs_1_theta(:,i)=0;
                
        rhs_2_u(:,i)=0;        
        rhs_2_w(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(:,i)=0;
        
        rhs_3_u(:,i)=C11*dgdx(c,x_dados(i),sub_dominio(:));
        rhs_3_w(:,i)=0;
        rhs_3_theta(:,i)=D11*dgdx(c,x_dados(i),sub_dominio(:));
        
        arhs_1_u(:,i)=0;
        arhs_1_w(:,i)=0;
        arhs_1_theta(:,i)=0;
        
        arhs_2_u(:,i)=0;
        arhs_2_w(:,i)=0;
        arhs_2_theta(:,i)=0;
        
        arhs_3_u(:,i)=0;
        arhs_3_w(:,i)=0;
        arhs_3_theta(:,i)=0;
            case {'cc'}
        rhs_1_u(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_1_w(:,i)=0;
        rhs_1_theta(:,i)=0;        
                
        rhs_2_u(:,i)=0;
        rhs_2_w(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(:,i)=0;
        
        rhs_3_u(:,i)=0;
        rhs_3_w(:,i)=0;
        rhs_3_theta(:,i)=g(c,x_dados(i),sub_dominio(:));
        
        arhs_1_u(:,i)=0;
        arhs_1_w(:,i)=0;
        arhs_1_theta(:,i)=0;
        
        arhs_2_u(:,i)=0;
        arhs_2_w(:,i)=0;
        arhs_2_theta(:,i)=0;
        
        arhs_3_u(:,i)=0;
        arhs_3_w(:,i)=0;
        arhs_3_theta(:,i)=0;
            case {'cl'}
        rhs_1_u(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_1_w(:,i)=0;
        rhs_1_theta(:,i)=0;        
                
        rhs_2_u(:,i)=0;
        rhs_2_w(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(:,i)=0;
        
        rhs_3_u(:,i)=0;
        rhs_3_w(:,i)=0;
        rhs_3_theta(:,i)=g(c,x_dados(i),sub_dominio(:));
        
        arhs_1_u(:,i)=0;
        arhs_1_w(:,i)=0;
        arhs_1_theta(:,i)=0;
        
        arhs_2_u(:,i)=0;
        arhs_2_w(:,i)=0;
        arhs_2_theta(:,i)=0;
        
        arhs_3_u(:,i)=0;
        arhs_3_w(:,i)=0;
        arhs_3_theta(:,i)=0;
        end
[Axi,Axj]=meshgrid(sub_dominio);
matriz_pesos=g(c,Axi,Axj);

pesos_1_u(i,:)=matriz_pesos\rhs_1_u(:,i);
pesos_1_w(i,:)=matriz_pesos\rhs_1_w(:,i);
pesos_1_theta(i,:)=matriz_pesos\rhs_1_theta(:,i);
pesos_2_u(i,:)=matriz_pesos\rhs_2_u(:,i);
pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);
pesos_3_u(i,:)=matriz_pesos\rhs_3_u(:,i);
pesos_3_w(i,:)=matriz_pesos\rhs_3_w(:,i);
pesos_3_theta(i,:)=matriz_pesos\rhs_3_theta(:,i);

apesos_1_u(i,:)=matriz_pesos\arhs_1_u(:,i);
apesos_1_w(i,:)=matriz_pesos\arhs_1_w(:,i);
apesos_1_theta(i,:)=matriz_pesos\arhs_1_theta(:,i);
apesos_2_u(i,:)=matriz_pesos\arhs_2_u(:,i);
apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);
apesos_3_u(i,:)=matriz_pesos\arhs_3_u(:,i);
apesos_3_w(i,:)=matriz_pesos\arhs_3_w(:,i);
apesos_3_theta(i,:)=matriz_pesos\arhs_3_theta(:,i);
        
    elseif  i==numel(x_dados)
        sub_dominio=[x_dados(i), x_dados(i-1), x_dados(i-2)];
         switch cfstr
            case {'ss'}
        rhs_1_u(:,i)=g(c,x_dados(i),sub_dominio(:));        
        rhs_1_w(:,i)=0;
        rhs_1_theta(:,i)=0;
                
        rhs_2_u(:,i)=0;        
        rhs_2_w(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(:,i)=0;
        
        rhs_3_u(:,i)=C11*dgdx(c,x_dados(i),sub_dominio(:));
        rhs_3_w(:,i)=0;
        rhs_3_theta(:,i)=D11*dgdx(c,x_dados(i),sub_dominio(:));
        
        arhs_1_u(:,i)=0;
        arhs_1_w(:,i)=0;
        arhs_1_theta(:,i)=0;
        
        arhs_2_u(:,i)=0;
        arhs_2_w(:,i)=0;
        arhs_2_theta(:,i)=0;
        
        arhs_3_u(:,i)=0;
        arhs_3_w(:,i)=0;
        arhs_3_theta(:,i)=0;
            case {'cc'}
        rhs_1_u(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_1_w(:,i)=0;
        rhs_1_theta(:,i)=0;        
                
        rhs_2_u(:,i)=0;
        rhs_2_w(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(:,i)=0;
        
        rhs_3_u(:,i)=0;
        rhs_3_w(:,i)=0;
        rhs_3_theta(:,i)=g(c,x_dados(i),sub_dominio(:));
        
        arhs_1_u(:,i)=0;
        arhs_1_w(:,i)=0;
        arhs_1_theta(:,i)=0;
        
        arhs_2_u(:,i)=0;
        arhs_2_w(:,i)=0;
        arhs_2_theta(:,i)=0;
        
        arhs_3_u(:,i)=0;
        arhs_3_w(:,i)=0;
        arhs_3_theta(:,i)=0;
        case {'cl'}
        rhs_1_u(:,i)=B11*dgdx(c,x_dados(i),sub_dominio(:));  
        rhs_1_w(:,i)=0;
        rhs_1_theta(:,i)=C11*dgdx(c,x_dados(i),sub_dominio(:));
        
        rhs_2_u(:,i)=0;
        rhs_2_w(:,i)=B55*dgdx(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(:,i)=B55*g(c,x_dados(i),sub_dominio(:));
        
        rhs_3_u(:,i)=C11*dgdx(c,x_dados(i),sub_dominio(:));
        rhs_3_w(:,i)=0;
        rhs_3_theta(:,i)=D11*dgdx(c,x_dados(i),sub_dominio(:));
        
        arhs_1_u(:,i)=0;
        arhs_1_w(:,i)=0;
        arhs_1_theta(:,i)=0;
        
        arhs_2_u(:,i)=0;
        arhs_2_w(:,i)=0;
        arhs_2_theta(:,i)=0;
        
        arhs_3_u(:,i)=0;
        arhs_3_w(:,i)=0;
        arhs_3_theta(:,i)=0;
         end
[Axi,Axj]=meshgrid(sub_dominio);
matriz_pesos=g(c,Axi,Axj);

pesos_1_u(i,:)=matriz_pesos\rhs_1_u(:,i);
pesos_1_w(i,:)=matriz_pesos\rhs_1_w(:,i);
pesos_1_theta(i,:)=matriz_pesos\rhs_1_theta(:,i);
pesos_2_u(i,:)=matriz_pesos\rhs_2_u(:,i);
pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);
pesos_3_u(i,:)=matriz_pesos\rhs_3_u(:,i);
pesos_3_w(i,:)=matriz_pesos\rhs_3_w(:,i);
pesos_3_theta(i,:)=matriz_pesos\rhs_3_theta(:,i);

apesos_1_u(i,:)=matriz_pesos\arhs_1_u(:,i);
apesos_1_w(i,:)=matriz_pesos\arhs_1_w(:,i);
apesos_1_theta(i,:)=matriz_pesos\arhs_1_theta(:,i);
apesos_2_u(i,:)=matriz_pesos\arhs_2_u(:,i);
apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);
apesos_3_u(i,:)=matriz_pesos\arhs_3_u(:,i);
apesos_3_w(i,:)=matriz_pesos\arhs_3_w(:,i);
apesos_3_theta(i,:)=matriz_pesos\arhs_3_theta(:,i);       
    else
        sub_dominio=[x_dados(i), x_dados(i-1), x_dados(i+1)];
        rhs_1_u(:,i)=B11*d2gdx2(c,x_dados(i),sub_dominio(:));
        rhs_1_w(:,i)=0;
        rhs_1_theta(:,i)=C11*d2gdx2(c,x_dados(i),sub_dominio(:));
        
        rhs_2_u(:,i)=0;
        rhs_2_w(:,i)=B55*d2gdx2(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(:,i)=B55*dgdx(c,x_dados(i),sub_dominio(:));
        
        rhs_3_u(:,i)=C11*d2gdx2(c,x_dados(i),sub_dominio(:));
        rhs_3_w(:,i)=-B55*dgdx(c,x_dados(i),sub_dominio(:));
        rhs_3_theta(:,i)=D11*d2gdx2(c,x_dados(i),sub_dominio(:))-B55*g(c,x_dados(i),sub_dominio(:));
        
        arhs_1_u(:,i)=-J0*g(c,x_dados(i),sub_dominio(:));
        arhs_1_w(:,i)=0;
        arhs_1_theta(:,i)=-J1*g(c,x_dados(i),sub_dominio(:));
        
        arhs_2_u(:,i)=0;
        arhs_2_w(:,i)=-J0*g(c,x_dados(i),sub_dominio(:));
        arhs_2_theta(:,i)=0;
        
        arhs_3_u(:,i)=-J1*g(c,x_dados(i),sub_dominio(:));
        arhs_3_w(:,i)=0;
        arhs_3_theta(:,i)=-J2*g(c,x_dados(i),sub_dominio(:));
        
[Axi,Axj]=meshgrid(sub_dominio);
matriz_pesos=g(c,Axi,Axj);

pesos_1_u(i,:)=matriz_pesos\rhs_1_u(:,i);
pesos_1_w(i,:)=matriz_pesos\rhs_1_w(:,i);
pesos_1_theta(i,:)=matriz_pesos\rhs_1_theta(:,i);
pesos_2_u(i,:)=matriz_pesos\rhs_2_u(:,i);
pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);
pesos_3_u(i,:)=matriz_pesos\rhs_3_u(:,i);
pesos_3_w(i,:)=matriz_pesos\rhs_3_w(:,i);
pesos_3_theta(i,:)=matriz_pesos\rhs_3_theta(:,i);

apesos_1_u(i,:)=matriz_pesos\arhs_1_u(:,i);
apesos_1_w(i,:)=matriz_pesos\arhs_1_w(:,i);
apesos_1_theta(i,:)=matriz_pesos\arhs_1_theta(:,i);
apesos_2_u(i,:)=matriz_pesos\arhs_2_u(:,i);
apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);
apesos_3_u(i,:)=matriz_pesos\arhs_3_u(:,i);
apesos_3_w(i,:)=matriz_pesos\arhs_3_w(:,i);
apesos_3_theta(i,:)=matriz_pesos\arhs_3_theta(:,i);        
        
end
end

%ASSEMBLAGEM DOS PESOS NAS MATRIZES DE INERCIA E DE RIGIDEZ
L_total=zeros(3*n,3*n);
A_total=zeros(3*n,3*n);
for i=2:n-1
    for j=i                
L_total(i,j-1)=pesos_1_u(i,2);        
L_total(i,j)=pesos_1_u(i,1);
L_total(i,j+1)=pesos_1_u(i,3);

L_total(i,n+j-1)=pesos_1_w(i,2);
L_total(i,n+j)=pesos_1_w(i,1);
L_total(i,n+j+1)=pesos_1_w(i,3);

L_total(i,2*n+j-1)=pesos_1_theta(i,2);
L_total(i,2*n+j)=pesos_1_theta(i,1);
L_total(i,2*n+j+1)=pesos_1_theta(i,3);

L_total(n+i,j-1)=pesos_2_u(i,2);
L_total(n+i,j)=pesos_2_u(i,1);
L_total(n+i,j+1)=pesos_2_u(i,3);

L_total(n+i,n+j-1)=pesos_2_w(i,2);
L_total(n+i,n+j)=pesos_2_w(i,1);
L_total(n+i,n+j+1)=pesos_2_w(i,3);

L_total(n+i,2*n+j-1)=pesos_2_theta(i,2);
L_total(n+i,2*n+j)=pesos_2_theta(i,1);
L_total(n+i,2*n+j+1)=pesos_2_theta(i,3);

L_total(2*n+i,j-1)=pesos_3_u(i,2);
L_total(2*n+i,j)=pesos_3_u(i,1);
L_total(2*n+i,j+1)=pesos_3_u(i,3);

L_total(2*n+i,n+j-1)=pesos_3_w(i,2);
L_total(2*n+i,n+j)=pesos_3_w(i,1);
L_total(2*n+i,n+j+1)=pesos_3_w(i,3);

L_total(2*n+i,2*n+j-1)=pesos_3_theta(i,2);
L_total(2*n+i,2*n+j)=pesos_3_theta(i,1);
L_total(2*n+i,2*n+j+1)=pesos_3_theta(i,3);

A_total(i,j-1)=apesos_1_u(i,2);        
A_total(i,j)=apesos_1_u(i,1);
A_total(i,j+1)=apesos_1_u(i,3);

A_total(i,n+j-1)=apesos_1_w(i,2);
A_total(i,n+j)=apesos_1_w(i,1);
A_total(i,n+j+1)=apesos_1_w(i,3);

A_total(i,2*n+j-1)=apesos_1_theta(i,2);
A_total(i,2*n+j)=apesos_1_theta(i,1);
A_total(i,2*n+j+1)=apesos_1_theta(i,3);

A_total(n+i,j-1)=apesos_2_u(i,2);
A_total(n+i,j)=apesos_2_u(i,1);
A_total(n+i,j+1)=apesos_2_u(i,3);

A_total(n+i,n+j-1)=apesos_2_w(i,2);
A_total(n+i,n+j)=apesos_2_w(i,1);
A_total(n+i,n+j+1)=apesos_2_w(i,3);

A_total(n+i,2*n+j-1)=apesos_2_theta(i,2);
A_total(n+i,2*n+j)=apesos_2_theta(i,1);
A_total(n+i,2*n+j+1)=apesos_2_theta(i,3);

A_total(2*n+i,j-1)=apesos_3_u(i,2);
A_total(2*n+i,j)=apesos_3_u(i,1);
A_total(2*n+i,j+1)=apesos_3_u(i,3);

A_total(2*n+i,n+j-1)=apesos_3_w(i,2);
A_total(2*n+i,n+j)=apesos_3_w(i,1);
A_total(2*n+i,n+j+1)=apesos_3_w(i,3);

A_total(2*n+i,2*n+j-1)=apesos_3_theta(i,2);
A_total(2*n+i,2*n+j)=apesos_3_theta(i,1);
A_total(2*n+i,2*n+j+1)=apesos_3_theta(i,3);
    end
end
L_total(1,1:3)=pesos_1_u(1,1:3);
L_total(n,n-2:n)=pesos_1_u(n,3:-1:1);

L_total(1,n+1:n+3)=pesos_1_w(1,1:3);
L_total(n,2*n-2:2*n)=pesos_1_w(n,3:-1:1);

L_total(1,2*n+1:2*n+3)=pesos_1_theta(1,1:3);
L_total(n,end-2:end)=pesos_1_theta(n,3:-1:1);

L_total(n+1,1:3)=pesos_2_u(1,1:3);
L_total(2*n,n-2:n)=pesos_2_u(n,3:-1:1);

L_total(n+1,n+1:n+3)=pesos_2_w(1,1:3);
L_total(2*n,2*n-2:2*n)=pesos_2_w(n,3:-1:1);

L_total(n+1,2*n+1:2*n+3)=pesos_2_theta(1,1:3);
L_total(2*n,end-2:end)=pesos_2_theta(n,3:-1:1);

L_total(2*n+1,1:3)=pesos_3_u(1,1:3);
L_total(end,n-2:n)=pesos_3_u(n,3:-1:1);

L_total(2*n+1,n+1:n+3)=pesos_3_w(1,1:3);
L_total(end,2*n-2:2*n)=pesos_3_w(n,3:-1:1);

L_total(2*n+1,2*n+1:2*n+3)=pesos_3_theta(1,1:3);
L_total(end,end-2:end)=pesos_3_theta(n,3:-1:1);


A_total(1,1:3)=apesos_1_u(1,1:3);
A_total(n,n-2:n)=apesos_1_u(n,3:-1:1);

A_total(1,n+1:n+3)=apesos_1_w(1,1:3);
A_total(n,2*n-2:2*n)=apesos_1_w(n,3:-1:1);

A_total(1,2*n+1:2*n+3)=apesos_1_theta(1,1:3);
A_total(n,end-2:end)=apesos_1_theta(n,3:-1:1);

A_total(n+1,1:3)=apesos_2_u(1,1:3);
A_total(2*n,n-2:n)=apesos_2_u(n,3:-1:1);

A_total(n+1,n+1:n+3)=apesos_2_w(1,1:3);
A_total(2*n,2*n-2:2*n)=apesos_2_w(n,3:-1:1);

A_total(n+1,2*n+1:2*n+3)=apesos_2_theta(1,1:3);
A_total(2*n,end-2:end)=apesos_2_theta(n,3:-1:1);

A_total(2*n+1,1:3)=apesos_3_u(1,1:3);
A_total(end,n-2:n)=apesos_3_u(n,3:-1:1);

A_total(2*n+1,n+1:n+3)=apesos_3_w(1,1:3);
A_total(end,2*n-2:2*n)=apesos_3_w(n,3:-1:1);

A_total(2*n+1,2*n+1:2*n+3)=apesos_3_theta(1,1:3);
A_total(end,end-2:end)=apesos_3_theta(n,3:-1:1);

%% EIGENVALUE PROBLEM
[lambda_vec,lambda]=eig(L_total,A_total);
[V,D]=eig(L_total,A_total,'qz');SS5=L_total*V-A_total*V*D;

T=L_total*lambda_vec-lambda*A_total*lambda_vec;

lambda=diag(lambda,0); 
[lambda,indice]=sort(lambda);
lambda_vec=lambda_vec(:,(indice(:)));

eigval = lambda(1,1);eigvec = lambda_vec(:,1);


m=1; E=6e10; I=I2(1)+I2(2); A=h; G=2.3e10;
sol_exacta=(m*pi/L)^2*sqrt((E*I)/(rho*A))*sqrt(1-(((m*pi/L)^2*E*I)/(k*G*A+(m*pi/L)^2*E*I)));
sol_exacta_norm=sol_exacta*L^2*sqrt(rho*A/(E*I));

p=1;
lambda_mode_w(1:n,p)=lambda_vec(n+1:2*n,p);
lambda_mode_phi_x(1:n,p)=lambda_vec(2*n+1:end,p);
lambda_mode=[lambda_mode_w;lambda_mode_phi_x];


freq_exacta=sol_exacta/(2*pi);
freq=sqrt(lambda(p))/(2*pi);
end


clear all
x_inicial=0;x_final=1;
contador_n=0;
for n=11:2:91;
    contador_n=contador_n+1;
x_dados=[0:1/(n-1):1];dist=x_dados(3)-x_dados(1);
[xi,xj]=meshgrid(x_dados);
x_central=find(x_dados==0.5);

%%
lsobreh=10;
h=1/lsobreh;
carga=1;A=1*h;k=5/6;
I=1*h^3/12;
E=1/I;rho=1/A;
G=E/(2*(1+0.25));

% rho=1;
% carga=1;A=1;k=5/6;
% E=1;I=1;G=E/(1+0.3);
L=x_final-x_inicial;


% solucao_exacta(1:n)=(-L.^3.*carga(:).*x_dados(:)+2.*carga(:).*x_dados(:).^3.*L-carga(:).*x_dados(:).^4)./(24.*E.*I)+(carga(:).*x_dados(:).^2-L.*carga(:).*x_dados(:))./(2.*G.*A.*k);
% solucao_exacta(n+1:2*n)=(-carga(:).*L.^3+3.*carga(:).*2.*x_dados(:).^2.*L-4.*carga(:).*x_dados(:).^3)./(24.*E.*I);

%viga simplesmente apoiada
rhs_total=zeros(2*n,1);
rhs_total(1:n,1)=carga;
rhs_total(1,1)=0;
rhs_total(n,1)=0;

rhs_1_w=zeros(3,numel(x_dados));
rhs_1_theta=zeros(3,numel(x_dados));

rhs_2_w=zeros(3,numel(x_dados));
rhs_2_theta=zeros(3,numel(x_dados));

arhs_1_w=zeros(3,numel(x_dados));
arhs_1_theta=zeros(3,numel(x_dados));

arhs_2_w=zeros(3,numel(x_dados));
arhs_2_theta=zeros(3,numel(x_dados));

contador_c=0;
for c=[1]%2*dist/sqrt(sqrt(3))
contador_c=contador_c+1;    
for i=1:numel(x_dados)
    if i==1
        sub_dominio=[x_dados(i), x_dados(i+1), x_dados(i+2)];
        rhs_1_w(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_1_theta(:,i)=0;
        
        rhs_2_w(:,i)=0;
        rhs_2_theta(:,i)=dgdx(c,x_dados(i),sub_dominio(:));
        
        arhs_1_w(:,i)=0;
        arhs_1_theta(:,i)=0;
        
        arhs_2_w(:,i)=0;
        arhs_2_theta(:,i)=0;
    
[Axi,Axj]=meshgrid(sub_dominio);
matriz_pesos=g(c,Axi,Axj);

pesos_1_w(i,:)=matriz_pesos\rhs_1_w(:,i);
pesos_1_theta(i,:)=matriz_pesos\rhs_1_theta(:,i);
pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);

apesos_1_w(i,:)=matriz_pesos\arhs_1_w(:,i);
apesos_1_theta(i,:)=matriz_pesos\arhs_1_theta(:,i);
apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);
        
    elseif  i==numel(x_dados)
        sub_dominio=[x_dados(i), x_dados(i-1), x_dados(i-2)];
        rhs_1_w(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_1_theta(:,i)=0;
        
        rhs_2_w(:,i)=0;
        rhs_2_theta(:,i)=dgdx(c,x_dados(i),sub_dominio(:));
        
        arhs_1_w(:,i)=0;
        arhs_1_theta(:,i)=0;
        
        arhs_2_w(:,i)=0;
        arhs_2_theta(:,i)=0;

[Axi,Axj]=meshgrid(sub_dominio);
matriz_pesos=g(c,Axi,Axj);

pesos_1_w(i,:)=matriz_pesos\rhs_1_w(:,i);
pesos_1_theta(i,:)=matriz_pesos\rhs_1_theta(:,i);
pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);

apesos_1_w(i,:)=matriz_pesos\arhs_1_w(:,i);
apesos_1_theta(i,:)=matriz_pesos\arhs_1_theta(:,i);
apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);        
    else
        
        sub_dominio=[x_dados(i), x_dados(i-1), x_dados(i+1)];
        rhs_1_w(:,i)=(G*A*k)*d2gdx2(c,x_dados(i),sub_dominio(:));
        rhs_1_theta(:,i)=(G*A*k)*dgdx(c,x_dados(i),sub_dominio(:));
        
        rhs_2_w(:,i)=-(G*A*k)*dgdx(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(:,i)=E*I*d2gdx2(c,x_dados(i),sub_dominio(:))-(G*A*k)*g(c,x_dados(i),sub_dominio(:));
        
        arhs_1_w(:,i)=-rho*A*g(c,x_dados(i),sub_dominio(:));
        arhs_1_theta(:,i)=0;
        
        arhs_2_w(:,i)=0;
        arhs_2_theta(:,i)=-rho*I*g(c,x_dados(i),sub_dominio(:));
        
[Axi,Axj]=meshgrid(sub_dominio);
matriz_pesos=g(c,Axi,Axj);

pesos_1_w(i,:)=matriz_pesos\rhs_1_w(:,i);
pesos_1_theta(i,:)=matriz_pesos\rhs_1_theta(:,i);
pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);

apesos_1_w(i,:)=matriz_pesos\arhs_1_w(:,i);
apesos_1_theta(i,:)=matriz_pesos\arhs_1_theta(:,i);
apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);       
        
    end
        


% [Axi,Axj]=meshgrid(sub_dominio);
% matriz_pesos=g(c,Axi,Axj);
% 
% pesos_1_w(i,:)=matriz_pesos\rhs_1_w(:,i);
% pesos_1_theta(i,:)=matriz_pesos\rhs_1_theta(:,i);
% pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
% pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);
% 
% apesos_1_w(i,:)=matriz_pesos\arhs_1_w(:,i);
% apesos_1_theta(i,:)=matriz_pesos\arhs_1_theta(:,i);
% apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
% apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);
end

L_total=zeros(2*n,2*n);
A_total=zeros(2*n,2*n);
for i=2:n-1
    for j=i
L_total(i,j-1)=pesos_1_w(i,2);        
L_total(i,j)=pesos_1_w(i,1);
L_total(i,j+1)=pesos_1_w(i,3);

L_total(i,n+j-1)=pesos_1_theta(i,2);
L_total(i,n+j)=pesos_1_theta(i,1);
L_total(i,n+j+1)=pesos_1_theta(i,3);

L_total(n+i,j-1)=pesos_2_w(i,2);
L_total(n+i,j)=pesos_2_w(i,1);
L_total(n+i,j+1)=pesos_2_w(i,3);

L_total(n+i,n+j-1)=pesos_2_theta(i,2);
L_total(n+i,n+j)=pesos_2_theta(i,1);
L_total(n+i,n+j+1)=pesos_2_theta(i,3);

A_total(i,j-1)=apesos_1_w(i,2);        
A_total(i,j)=apesos_1_w(i,1);
A_total(i,j+1)=apesos_1_w(i,3);

A_total(i,n+j-1)=apesos_1_theta(i,2);
A_total(i,n+j)=apesos_1_theta(i,1);
A_total(i,n+j+1)=apesos_1_theta(i,3);

A_total(n+i,j-1)=apesos_2_w(i,2);
A_total(n+i,j)=apesos_2_w(i,1);
A_total(n+i,j+1)=apesos_2_w(i,3);

A_total(n+i,n+j-1)=apesos_2_theta(i,2);
A_total(n+i,n+j)=apesos_2_theta(i,1);
A_total(n+i,n+j+1)=apesos_2_theta(i,3);

    end
end
L_total(1,1:3)=pesos_1_w(1,1:3);
L_total(n,n-2:n)=pesos_1_w(n,3:-1:1);

L_total(1,n+1:n+3)=pesos_1_theta(1,1:3);
L_total(n,end-2:end)=pesos_1_theta(n,3:-1:1);

L_total(n+1,1:3)=pesos_2_w(1,1:3);
L_total(end,n-2:n)=pesos_2_w(n,3:-1:1);

L_total(n+1,n+1:n+3)=pesos_2_theta(1,1:3);
L_total(end,end-2:end)=pesos_2_theta(n,3:-1:1);

A_total(1,1:3)=apesos_1_w(1,1:3);
A_total(n,n-2:n)=apesos_1_w(n,3:-1:1);

A_total(1,n+1:n+3)=apesos_1_theta(1,1:3);
A_total(n,end-2:end)=apesos_1_theta(n,3:-1:1);

A_total(n+1,1:3)=apesos_2_w(1,1:3);
A_total(end,n-2:n)=apesos_2_w(n,3:-1:1);

A_total(n+1,n+1:n+3)=apesos_2_theta(1,1:3);
A_total(end,end-2:end)=apesos_2_theta(n,3:-1:1);

% solucao_final=L_total\rhs_total;



% figure(3)
% subplot(1,2,1);plot(x_dados,solucao_exacta(1:n));hold on
% subplot(1,2,1);plot(x_dados,solucao_final(1:n)','r.');hold on
% subplot(1,2,2);plot(x_dados,solucao_exacta(n+1:end));hold on
% subplot(1,2,2);plot(x_dados,solucao_final(n+1:end)','r.');hold on

dt=10^-5;
clear gg
gg=x_central;

%_______________condiçoes iniciais______________

forca=zeros(2*n,1);
deslocamentos_inicial(1:2*n,1)=0;
velocidade_inicial(1:2*n,1)=0;
aceleracao_inicial=((-L_total*deslocamentos_inicial+rhs_total)\A_total)';
deslocamentos_inicial(1:2*n,1)=0;

%____________parametros Newmark____reddy____________
alfa=3/2;gama=8/5;
a1=(1-alfa)*dt;a2=alfa*dt;a3=2/(gama*dt^2);a4=a3*dt;a5=(1-gama)/gama;

 aceleracao=aceleracao_inicial;
 velocidade=velocidade_inicial;
 deslocamentos=deslocamentos_inicial;
 
dmatriz_chapeu=L_total+a3*A_total;
cond_number(1,contador_c)=rcond(dmatriz_chapeu);

for passo=1:70000
t=passo*dt;

forca=rhs_total';

forca_chapeu=rhs_total+A_total*(a3*deslocamentos+a4*velocidade+a5*aceleracao);

deslocamentos_new=dmatriz_chapeu\forca_chapeu;

aceleracao_new=a3*(deslocamentos_new-deslocamentos)-a4*velocidade-a5*aceleracao;
velocidade_new=velocidade+a1*aceleracao+a2*aceleracao_new;

%-----------------------------------
aceleracao=aceleracao_new;
velocidade=velocidade_new;
deslocamentos=deslocamentos_new;

 w_final(passo+1,1,contador_c)=t;
 w_final(passo+1,2,contador_c)=deslocamentos(gg);
 w_final(passo+1,3,contador_c)=c;
 
 
end

w_exact=-0.0268;t_exact=32242; %navier
erro1(contador_c)=abs((w_exact-w_final(t_exact,2,contador_c))/w_exact)*100;
end
N(contador_n)=n;
erro_n(contador_n)=abs((w_exact-w_final(t_exact,2,contador_c))/w_exact)*100;
end
figure(3);
AA=w_final(32242,3,:);
plot(AA(:),erro1(:),'r.');hold on

% figure(4);plot(w_final(1:500:end,1),w_final(1:500:end,2),'yellow.');hold on
% figure(4)
% plot(w_final(:,1),w_final(:,2),'black.');hold on


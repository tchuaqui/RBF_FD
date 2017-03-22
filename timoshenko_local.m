clear all
x_inicial=0;x_final=1;
%  x_dados = x_inicial + (x_final-(x_inicial)) * rand(1,18);
%  x_dados=unique([0,x_dados,0.5,1]);
contador=1;
% n_total=numel(x_dados);
% x_dados_reg=[0:1/(n_total-1):1];
for n=21
x_dados=[0:1/(n-1):1];
[xi,xj]=meshgrid(x_dados);
x_central=find(x_dados==0.5);
%%
for lsobreh=1:10:50;
h=1/lsobreh;
carga=1;A=1*h;k=5/6;
I=1*h^3/12;
E=1;rho=1/A;
G=E/(2*(1+0.25));
% carga=1;A=1;k=5/6;
% E=1;I=1;G=E/(2*(1+0.3));

L=x_final-x_inicial;
solucao_exacta(1:n)=(-L.^3.*carga.*x_dados(:)+2.*carga.*x_dados(:).^3.*L-carga.*x_dados(:).^4)./(24.*E.*I)+(carga.*x_dados(:).^2-L.*carga.*x_dados(:))./(2.*G.*A.*k);
solucao_exacta(n+1:2*n)=(-carga.*L.^3+3.*carga.*2.*x_dados(:).^2.*L-4.*carga.*x_dados(:).^3)./(24.*E.*I);

% solucao_exacta(1:n)=(-L.^3.*carga.*x_dados_reg(:)+2.*carga.*x_dados_reg(:).^3.*L-carga.*x_dados_reg(:).^4)./(24.*E.*I)+(carga.*x_dados_reg(:).^2-L.*carga.*x_dados_reg(:))./(2.*G.*A.*k);
% solucao_exacta(n+1:2*n)=(-carga.*L.^3+3.*carga.*2.*x_dados_reg(:).^2.*L-4.*carga.*x_dados_reg(:).^3)./(24.*E.*I);


%viga simplesmente apoiada
rhs_total=zeros(2*n,1);
rhs_total(1:n,1)=carga;
rhs_total(1,1)=0;
rhs_total(n,1)=0;

rhs_1_w=zeros(3,numel(x_dados));
rhs_1_theta=zeros(3,numel(x_dados));

rhs_2_w=zeros(3,numel(x_dados));
rhs_2_theta=zeros(3,numel(x_dados));

for c=2/sqrt(3)
    
for i=1:numel(x_dados)
    if i==1
        sub_dominio=[x_dados(i), x_dados(i+1), x_dados(i+2)];
        rhs_1_w(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_1_theta(:,i)=0;
        
        rhs_2_w(:,i)=0;
        rhs_2_theta(:,i)=dgdx(c,x_dados(i),sub_dominio(:));
        
    elseif  i==numel(x_dados)
        sub_dominio=[x_dados(i), x_dados(i-1), x_dados(i-2)];
        rhs_1_w(:,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_1_theta(:,i)=0;
        
        rhs_2_w(:,i)=0;
        rhs_2_theta(:,i)=dgdx(c,x_dados(i),sub_dominio(:));
        
    else
        
        sub_dominio=[x_dados(i), x_dados(i-1), x_dados(i+1)];
        rhs_1_w(:,i)=(G*A*k)*d2gdx2(c,x_dados(i),sub_dominio(:));
        rhs_1_theta(:,i)=(G*A*k)*dgdx(c,x_dados(i),sub_dominio(:));
        
        rhs_2_w(:,i)=-(G*A*k)*dgdx(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(:,i)=E*I*d2gdx2(c,x_dados(i),sub_dominio(:))-(G*A*k)*g(c,x_dados(i),sub_dominio(:));
    end
    
[Axi,Axj]=meshgrid(sub_dominio);
matriz_pesos=g(c,Axi,Axj);
pesos_1_w(i,:)=matriz_pesos\rhs_1_w(:,i);
pesos_1_theta(i,:)=matriz_pesos\rhs_1_theta(:,i);
pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);
end

L_total=zeros(2*n,2*n);

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

solucao_final=L_total\rhs_total;


num_pontos(contador)=numel(x_dados);
erro_rel_deformada(contador)=100*abs(solucao_exacta(x_central)-solucao_final(x_central))/abs(solucao_exacta(x_central));
shape(contador)=c;
num_cond(contador)=rcond(L_total);
max_disp(contador)=solucao_final(x_central);
asobreh(contador)=lsobreh;
contador=contador+1; 

end
end
end

%  figure(4)
%  plot(num_pontos,erro_rel_deformada,'rx');hold on

 figure(44)
 plot(asobreh,max_disp,'r.');hold on

figure(1)
plot(asobreh,num_cond,'b.');hold on
% 
% figure(3)
% plot(shape,erro_rel_deformada,'rx');hold on

% figure(2)
% %subplot(1,2,1);plot(x_dados,solucao_exacta(1:n));hold on
% subplot(1,2,1);plot(x_dados,solucao_final(1:n)','r.');hold on
% %subplot(1,2,2);plot(x_dados,solucao_exacta(n+1:end));hold on
% subplot(1,2,2);plot(x_dados,solucao_final(n+1:end)','r.');hold on

function [ freq ] = conv( n,p,c,mode,st) %cf,ncamadas, distribuicao quad ou lin do potencial
x_inicial=0;x_final=100e-3;
L=x_final-x_inicial;
cfstr='cl';
ncamadas=2;    %n camadas
esp_al=8e-3; esp_piezo=1e-3;
x_dados=[0:L/(n-1):L];dist=x_dados(3)-x_dados(1);
[xi,xj]=meshgrid(x_dados);
x_central=find(x_dados==0.5);
[ pol,dpol,d2pol ] = polynomials( x_dados,n,p);%%
%%
k=5/6;
switch ncamadas
    case 2
propmec=[5e-3 6e10 2.3e10 7500]; %[esp E G rho]
propel=[-16.492145 -16.492145 2.588549e-8]; %[e31 e32 ezz]
esp=[propmec(1) propmec(1)];
rho=[propmec(4) propmec(4)];

e31=propel(1);
ezz=propel(3);

Q11=[propmec(2) propmec(2)];
Q55=[propmec(3) propmec(3)];
    case 3
propmec=[esp_piezo 6e10 2.3e10 7500]; %[esp E G rho] %propriedades do piezo
propel=[-16.492145 -16.492145 2.588549e-8]; %[e31 e32 ezz]
propal=[esp_al 70.3e9 0.345 2710]; %[esp E mu rho]
esp=[propmec(1) propal(1) propmec(1)];
rho=[propmec(4) propal(4) propmec(4)];
e31=propel(1);
ezz=propel(3);
Q11=[propmec(2) propal(2) propmec(2)];
Q55=[propmec(3) propal(2)/(2*(1+propal(3))) propmec(3)];
end
%%
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

J0=J0+rho(i)*I0(i);
J1=J1+rho(i)*I1(i);
J2=J2+rho(i)*I2(i);

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
% G=0; E=0;                                                       %distribuiçao linear do potencial
G=(e31^2)*(I2(1)-I0(1)*zm(1)^2)/ezz+(e31^2)*(I2(end)-I0(end)*zm(end)^2)/ezz; 
E=(e31^2)*(I1(1)-I0(1)*zm(1))/ezz+(e31^2)*(I1(end)-I0(end)*zm(end))/ezz;  %distribuiçao quadratica do potencial


F1=e31*I0(1)/esp(1);
F2=e31*I0(end)/esp(end);
H1=e31*I1(1)/esp(1);
H2=e31*I1(end)/esp(end);
I_1=ezz*I0(1)/(esp(1)^2);
I_2=ezz*I0(end)/(esp(end)^2);

%% INICIACAO DAS MATRIZES (PARA DETERMINACAO DOS PESOS) RIGIDEZ E INERCIA

rhs_2_w=zeros(st+p,numel(x_dados));
rhs_2_theta=zeros(st+p,numel(x_dados));
rhs_2_phis=zeros(st+p,numel(x_dados));   
rhs_2_phia=zeros(st+p,numel(x_dados));   
rhs_3_w=zeros(st+p,numel(x_dados));     
rhs_3_theta=zeros(st+p,numel(x_dados));  
rhs_3_phis=zeros(st+p,numel(x_dados));       
rhs_3_phia=zeros(st+p,numel(x_dados));      

rhs_phiphis=zeros(st+p,numel(x_dados));  
rhs_phiphia=zeros(st+p,numel(x_dados));   

arhs_2_w=zeros(st+p,numel(x_dados));
arhs_2_theta=zeros(st+p,numel(x_dados));
arhs_3_w=zeros(st+p,numel(x_dados));                 
arhs_3_theta=zeros(st+p,numel(x_dados));   


%% CONDICOES DE FRONTEIRA E DOMINIO 
for i=1:numel(x_dados)
if i==1
        sub_dominio=x_dados(1:st);
        switch cfstr
            case {'ss'}     
        rhs_2_w(1:st,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(1:end,i)=0;
        rhs_3_w(1:end,i)=0;
        rhs_3_theta(1:st,i)=(D11+G-(1/B11)*(C11+E)^2)*dgdx(c,x_dados(i),sub_dominio(:));

        arhs_2_w(1:end,i)=0;
        arhs_2_theta(1:end,i)=0;
        arhs_3_w(1:end,i)=0;
        arhs_3_theta(1:end,i)=0;
        
        rhs_3_phis(1:st,i)=H1*g(c,x_dados(i),sub_dominio(:));   
        rhs_3_phia(1:st,i)=H2*g(c,x_dados(i),sub_dominio(:));   
        rhs_phiphis(1:st,i)=I_1*g(c,x_dados(i),sub_dominio(:));  
        rhs_phiphia(1:st,i)=I_2*g(c,x_dados(i),sub_dominio(:));  
        
        if p~=0
        rhs_2_w(st+1:end,i)=pol(1:p,i);
        rhs_3_theta(st+1:end,i)=(D11+G-(1/B11)*(C11+E)^2)*dpol(1:p,i);
        rhs_3_phis(st+1:end,i)=H1*pol(1:p,i);
        rhs_3_phia(st+1:end,i)=H2*pol(1:p,i);
        rhs_phiphis(st+1:end,i)=I_1*pol(1:p,i);
        rhs_phiphia(st+1:end,i)=I_2*pol(1:p,i);
        end
            case {'cc'}
        rhs_2_w(1:st,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(1:end,i)=0;   
        rhs_3_w(1:end,i)=0;
        rhs_3_theta(1:st,i)=g(c,x_dados(i),sub_dominio(:));
        
        arhs_2_w(1:end,i)=0;
        arhs_2_theta(1:end,i)=0;
        arhs_3_w(1:end,i)=0;
        arhs_3_theta(1:end,i)=0;
  
        rhs_3_phis(1:end,i)=0;  
        rhs_3_phia(1:end,i)=0;  
        rhs_phiphis(1:st,i)=I_1*g(c,x_dados(i),sub_dominio(:));  
        rhs_phiphia(1:st,i)=I_2*g(c,x_dados(i),sub_dominio(:));  
        
        if p~=0
        rhs_2_w(st+1:end,i)=pol(1:p,i);
        rhs_3_theta(st+1:end,i)=pol(1:p,i);
        rhs_phiphis(st+1:end,i)=I_1*pol(1:p,i);
        rhs_phiphia(st+1:end,i)=I_2*pol(1:p,i);
        end
            case {'cl'}
        rhs_2_w(1:st,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(1:end,i)=0;
        rhs_3_w(1:end,i)=0;
        rhs_3_theta(1:st,i)=g(c,x_dados(i),sub_dominio(:));
        
        arhs_2_w(1:end,i)=0;
        arhs_2_theta(1:end,i)=0;
        arhs_3_w(1:end,i)=0;
        arhs_3_theta(1:end,i)=0;

        rhs_3_phis(1:end,i)=0;   
        rhs_3_phia(1:end,i)=0;   
        rhs_phiphis(1:st,i)=I_1*g(c,x_dados(i),sub_dominio(:));  
        rhs_phiphia(1:st,i)=I_2*g(c,x_dados(i),sub_dominio(:));  
        
        if p~=0
        rhs_2_w(st+1:end,i)=pol(1:p,i);
        rhs_3_theta(st+1:end,i)=pol(1:p,i);
        rhs_phiphis(st+1:end,i)=I_1*pol(1:p,i);
        rhs_phiphia(st+1:end,i)=I_2*pol(1:p,i);
        end
        end
[Axi,Axj]=meshgrid(sub_dominio);
m_aux1=pol(1:p,1:st); m_aux2=zeros(p,p);
matriz_pesos=[g(c,Axi,Axj) m_aux1';    
              m_aux1 m_aux2]; 

pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);
pesos_3_w(i,:)=matriz_pesos\rhs_3_w(:,i);
pesos_3_theta(i,:)=matriz_pesos\rhs_3_theta(:,i);

apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);
apesos_3_w(i,:)=matriz_pesos\arhs_3_w(:,i);
apesos_3_theta(i,:)=matriz_pesos\arhs_3_theta(:,i);
  
pesos_3_phis(i,:)=matriz_pesos\rhs_3_phis(:,i);   
pesos_3_phia(i,:)=matriz_pesos\rhs_3_phia(:,i);   
pesos_phiphis(i,:)=matriz_pesos\rhs_phiphis(:,i); 
pesos_phiphia(i,:)=matriz_pesos\rhs_phiphia(:,i);        
elseif  i==numel(x_dados)
        sub_dominio=x_dados(n-st+1:n);
        switch cfstr
            case {'ss'}        
        rhs_2_w(1:st,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(1:end,i)=0;
        rhs_3_w(1:end,i)=0;
        rhs_3_theta(1:st,i)=(D11+G-(1/B11)*(C11+E)^2)*dgdx(c,x_dados(i),sub_dominio(:));
        
        arhs_2_w(1:end,i)=0;
        arhs_2_theta(1:end,i)=0;
        arhs_3_w(1:end,i)=0;
        arhs_3_theta(1:end,i)=0;

        rhs_3_phis(1:st,i)=H1*g(c,x_dados(i),sub_dominio(:));   
        rhs_3_phia(1:st,i)=H2*g(c,x_dados(i),sub_dominio(:));   
        rhs_phiphis(1:st,i)=I_1*g(c,x_dados(i),sub_dominio(:));  
        rhs_phiphia(1:st,i)=I_2*g(c,x_dados(i),sub_dominio(:));  
        
        if p~=0
        rhs_2_w(st+1:end,i)=pol(1:p,i);
        rhs_3_theta(st+1:end,i)=(D11+G-(1/B11)*(C11+E)^2)*dpol(1:p,i);
        rhs_3_phis(st+1:end,i)=H1*pol(1:p,i);
        rhs_3_phia(st+1:end,i)=H2*pol(1:p,i);
        rhs_phiphis(st+1:end,i)=I_1*pol(1:p,i);
        rhs_phiphia(st+1:end,i)=I_2*pol(1:p,i);
        end
            case {'cc'}
        rhs_2_w(1:st,i)=g(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(1:end,i)=0;
        rhs_3_w(1:end,i)=0;
        rhs_3_theta(1:st,i)=g(c,x_dados(i),sub_dominio(:));

        arhs_2_w(1:end,i)=0;
        arhs_2_theta(1:end,i)=0;
        arhs_3_w(1:end,i)=0;
        arhs_3_theta(1:end,i)=0;
        
        rhs_3_phis(1:end,i)=0;   
        rhs_3_phia(1:end,i)=0;   
        rhs_phiphis(1:st,i)=I_1*g(c,x_dados(i),sub_dominio(:));  
        rhs_phiphia(1:st,i)=I_2*g(c,x_dados(i),sub_dominio(:));  
        
        if p~=0
        rhs_2_w(st+1:end,i)=pol(1:p,i);
        rhs_3_theta(st+1:end,i)=pol(1:p,i);
        rhs_phiphis(st+1:end,i)=I_1*pol(1:p,i);
        rhs_phiphia(st+1:end,i)=I_2*pol(1:p,i);
        end
        case {'cl'}
        rhs_2_w(1:st,i)=B55*dgdx(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(1:st,i)=B55*g(c,x_dados(i),sub_dominio(:));
        rhs_3_w(1:end,i)=0;
        rhs_3_theta(1:st,i)=(D11+G-(1/B11)*(C11+E)^2)*dgdx(c,x_dados(i),sub_dominio(:));
        
        arhs_2_w(1:end,i)=0;
        arhs_2_theta(1:end,i)=0;
        arhs_3_w(1:end,i)=0;
        arhs_3_theta(1:end,i)=0;
          
        rhs_3_phis(1:st,i)=H1*g(c,x_dados(i),sub_dominio(:));   
        rhs_3_phia(1:st,i)=H2*g(c,x_dados(i),sub_dominio(:));   
        rhs_phiphis(1:st,i)=I_1*g(c,x_dados(i),sub_dominio(:));  
        rhs_phiphia(1:st,i)=I_2*g(c,x_dados(i),sub_dominio(:));  
        
        if p~=0
        rhs_2_w(st+1:end,i)=B55*dpol(1:p,i);
        rhs_2_theta(st+1:end,i)=B55*pol(1:p,i);
        rhs_3_theta(st+1:end,i)=(D11+G-(1/B11)*(C11+E)^2)*dpol(1:p,i);
        rhs_3_phis(st+1:end,i)=H1*pol(1:p,i);
        rhs_3_phia(st+1:end,i)=H2*pol(1:p,i);
        rhs_phiphis(st+1:end,i)=I_1*pol(1:p,i);
        rhs_phiphia(st+1:end,i)=I_2*pol(1:p,i);
        end
        end
[Axi,Axj]=meshgrid(sub_dominio);
m_aux1=pol(1:p,i-st+1:i); m_aux2=zeros(p,p);
matriz_pesos=[g(c,Axi,Axj) m_aux1';    
              m_aux1 m_aux2]; 

pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);
pesos_3_w(i,:)=matriz_pesos\rhs_3_w(:,i);
pesos_3_theta(i,:)=matriz_pesos\rhs_3_theta(:,i);

apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);
apesos_3_w(i,:)=matriz_pesos\arhs_3_w(:,i);
apesos_3_theta(i,:)=matriz_pesos\arhs_3_theta(:,i);    
 
pesos_3_phis(i,:)=matriz_pesos\rhs_3_phis(:,i);   
pesos_3_phia(i,:)=matriz_pesos\rhs_3_phia(:,i);   
pesos_phiphis(i,:)=matriz_pesos\rhs_phiphis(:,i);  
pesos_phiphia(i,:)=matriz_pesos\rhs_phiphia(:,i);         

else
    if i<=ceil(st/2)
        sub_dominio=x_dados(1:st);
        m_aux1=pol(1:p,1:st); 
    elseif i>ceil(st/2) && i<n-floor(st/2)
        sub_dominio=x_dados(i-floor(st/2):i+floor(st/2));
        m_aux1=pol(1:p,i-floor(st/2):i+floor(st/2)); 
    else
        sub_dominio=x_dados(n-st+1:n);
        m_aux1=pol(1:p,n-st+1:n); 
    end
        rhs_2_w(1:st,i)=B55*d2gdx2(c,x_dados(i),sub_dominio(:));
        rhs_2_theta(1:st,i)=B55*dgdx(c,x_dados(i),sub_dominio(:));
        rhs_3_w(1:st,i)=-B55*dgdx(c,x_dados(i),sub_dominio(:));
        rhs_3_theta(1:st,i)=(D11+G-(1/B11)*(C11+E)^2)*d2gdx2(c,x_dados(i),sub_dominio(:))-B55*g(c,x_dados(i),sub_dominio(:));
        
        arhs_2_w(1:st,i)=-J0*g(c,x_dados(i),sub_dominio(:));
        arhs_2_theta(1:end,i)=0;
        arhs_3_w(1:end,i)=0;
        arhs_3_theta(1:st,i)=-J2*g(c,x_dados(i),sub_dominio(:));
          
        rhs_3_phis(1:st,i)=H1*dgdx(c,x_dados(i),sub_dominio(:));   
        rhs_3_phia(1:st,i)=H2*dgdx(c,x_dados(i),sub_dominio(:));  
        rhs_phiphis(1:st,i)=I_1*g(c,x_dados(i),sub_dominio(:));  
        rhs_phiphia(1:st,i)=I_2*g(c,x_dados(i),sub_dominio(:)); 
        
        if p~=0
        rhs_2_w(st+1:end,i)=B55*d2pol(1:p,i);
        rhs_2_theta(st+1:end,i)=B55*dpol(1:p,i);
        rhs_3_w(st+1:end,i)=-B55*dpol(1:p,i);
        rhs_3_theta(st+1:end,i)=(D11+G-(1/B11)*(C11+E)^2)*d2pol(1:p,i)-B55*pol(1:p,i);
        rhs_3_phis(st+1:end,i)=H1*dpol(1:p,i);
        rhs_3_phia(st+1:end,i)=H2*dpol(1:p,i);
        rhs_phiphis(st+1:end,i)=I_1*pol(1:p,i);
        rhs_phiphia(st+1:end,i)=I_2*pol(1:p,i);
        arhs_2_w(st+1:end,i)=-J0*pol(1:p,i);
        arhs_3_theta(st+1:end,i)=-J2*pol(1:p,i);
        end
[Axi,Axj]=meshgrid(sub_dominio);
m_aux2=zeros(p,p);
matriz_pesos=[g(c,Axi,Axj) m_aux1';    
              m_aux1 m_aux2]; 

pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);
pesos_3_w(i,:)=matriz_pesos\rhs_3_w(:,i);
pesos_3_theta(i,:)=matriz_pesos\rhs_3_theta(:,i);

apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);
apesos_3_w(i,:)=matriz_pesos\arhs_3_w(:,i);
apesos_3_theta(i,:)=matriz_pesos\arhs_3_theta(:,i);   

pesos_3_phis(i,:)=matriz_pesos\rhs_3_phis(:,i);  
pesos_3_phia(i,:)=matriz_pesos\rhs_3_phia(:,i);  
pesos_phiphis(i,:)=matriz_pesos\rhs_phiphis(:,i);  
pesos_phiphia(i,:)=matriz_pesos\rhs_phiphia(:,i);       
        
end
end

%ASSEMBLAGEM DOS PESOS NAS MATRIZES DE INERCIA E DE RIGIDEZ
L_total=zeros(2*n,2*n);
A_total=zeros(2*n,2*n);

for i=1:n
 if i<=ceil(st/2)
  L_total(i,1:st)=pesos_2_w(i,1:st);
  L_total(i,n+1:n+st)=pesos_2_theta(i,1:st);
  L_total(n+i,1:st)=pesos_3_w(i,1:st);
  L_total(n+i,n+1:n+st)=pesos_3_theta(i,1:st);
  
  A_total(i,1:st)=apesos_2_w(i,1:st);
  A_total(i,n+1:n+st)=apesos_2_theta(i,1:st);
  A_total(n+i,1:st)=apesos_3_w(i,1:st);
  A_total(n+i,n+1:n+st)=apesos_3_theta(i,1:st);
 elseif i>ceil(st/2) && i<n-floor(st/2)
  L_total(i,i-floor(st/2):i+floor(st/2))=pesos_2_w(i,1:st);
  L_total(i,n+i-floor(st/2):n+i+floor(st/2))=pesos_2_theta(i,1:st);
  L_total(n+i,i-floor(st/2):i+floor(st/2))=pesos_3_w(i,1:st);
  L_total(n+i,n+i-floor(st/2):n+i+floor(st/2))=pesos_3_theta(i,1:st);
  
  A_total(i,i-floor(st/2):i+floor(st/2))=apesos_2_w(i,1:st);
  A_total(i,n+i-floor(st/2):n+i+floor(st/2))=apesos_2_theta(i,1:st);
  A_total(n+i,i-floor(st/2):i+floor(st/2))=apesos_3_w(i,1:st);
  A_total(n+i,n+i-floor(st/2):n+i+floor(st/2))=apesos_3_theta(i,1:st);
 else
  L_total(i,n-st+1:n)=pesos_2_w(i,1:st);
  L_total(i,end-st+1:end)=pesos_2_theta(i,1:st);
  L_total(n+i,n-st+1:n)=pesos_3_w(i,1:st);
  L_total(n+i,end-st+1:end)=pesos_3_theta(i,1:st);
  
  A_total(i,n-st+1:n)=apesos_2_w(i,1:st);
  A_total(i,end-st+1:end)=apesos_2_theta(i,1:st);
  A_total(n+i,n-st+1:n)=apesos_3_w(i,1:st);
  A_total(n+i,end-st+1:end)=apesos_3_theta(i,1:st);
 end
end

%% EIGENVALUE PROBLEM
[lambda_vec,lambda]=eig(L_total,A_total);
[V,D]=eig(L_total,A_total,'qz');SS5=L_total*V-A_total*V*D;

T=L_total*lambda_vec-lambda*A_total*lambda_vec;

lambda=diag(lambda,0); 
[lambda,indice]=sort(lambda);
lambda_vec=lambda_vec(:,(indice(:)));

eigval = lambda(1,1);eigvec = lambda_vec(:,1);

% 
% m=1; E=6e10; I=I2(1)+I2(2); A=h; G=2.3e10;
% sol_exacta=(m*pi/L)^2*sqrt((E*I)/(rho*A))*sqrt(1-(((m*pi/L)^2*E*I)/(k*G*A+(m*pi/L)^2*E*I)));
% sol_exacta_norm=sol_exacta*L^2*sqrt(rho*A/(E*I));

pp=mode;
lambda_mode_w(1:n,pp)=lambda_vec(1:n,pp);
lambda_mode_phi_x(1:n,pp)=lambda_vec(n+1:2*n,pp);
lambda_mode=[lambda_mode_w;lambda_mode_phi_x];


% freq_exacta=sol_exacta/(2*pi);
freq=sqrt(lambda(pp))/(2*pi);
end


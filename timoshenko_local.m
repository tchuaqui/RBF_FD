clear all
x_inicial=0;x_final=100e-3;
L=x_final-x_inicial;
cfstr='cl';
n=5;
x_dados=[0:L/(n-1):L];dist=x_dados(3)-x_dados(1);
[xi,xj]=meshgrid(x_dados);
x_central=find(x_dados==0.5);
carga=1000;
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
% G=0;                                                        %distribuiçao linear do potencial
G=(e31^2)*(I2(1)-I0(1)*zm(1)^2)/ezz+(e31^2)*(I2(2)-I0(2)*zm(2)^2)/ezz;   %distribuiçao quadratica do potencial

F1=e31*I0(1)/esp(1);
F2=e31*I0(2)/esp(2);
H1=e31*I1(1)/esp(1);
H2=e31*I1(2)/esp(2);
I_1=ezz*I0(1)/(esp(1)^2);
I_2=ezz*I0(2)/(esp(2)^2);




%% INICIACAO DAS MATRIZES (PARA DETERMINACAO DOS PESOS) RIGIDEZ E INERCIA
rhs_1_w=[zeros(4,numel(x_dados))];                  
rhs_1_theta=[zeros(4,numel(x_dados))];
rhs_2_w=[zeros(4,numel(x_dados))];
rhs_2_theta=[zeros(4,numel(x_dados))];
rhs_2_phis=[zeros(4,numel(x_dados))];
rhs_2_phia=[zeros(4,numel(x_dados))];
arhs_1_w=[zeros(4,numel(x_dados))];
arhs_1_theta=[zeros(4,numel(x_dados))];
arhs_2_w=[zeros(4,numel(x_dados))];
arhs_2_theta=[zeros(4,numel(x_dados))];
rhs_phiphis=[zeros(4,numel(x_dados))];
rhs_phiphia=[zeros(4,numel(x_dados))];


% 
% rhs_1_w=[zeros(5,numel(x_dados))];                  
% rhs_1_theta=[zeros(5,numel(x_dados))];
% rhs_2_w=[zeros(5,numel(x_dados))];
% rhs_2_theta=[zeros(5,numel(x_dados))];
% rhs_2_phis=[zeros(5,numel(x_dados))];
% rhs_2_phia=[zeros(5,numel(x_dados))];
% arhs_1_w=[zeros(5,numel(x_dados))];
% arhs_1_theta=[zeros(5,numel(x_dados))];
% arhs_2_w=[zeros(5,numel(x_dados))];
% arhs_2_theta=[zeros(5,numel(x_dados))];
% rhs_phiphis=[zeros(5,numel(x_dados))];
% rhs_phiphia=[zeros(5,numel(x_dados))];

c=2;  %2*dist/sqrt(sqrt(3))        %SHAPE PARAMETER   

%% CONDICOES DE FRONTEIRA E DOMINIO
for i=1:numel(x_dados)
if i==1
        sub_dominio=[x_dados(i), x_dados(i+1), x_dados(i+2)];
        switch cfstr
            case {'ss'}
        rhs_1_w(1:3,i)=g(c,x_dados(i),sub_dominio(:));           
        rhs_1_w(4,i)=1;
%         rhs_1_w(5,i)=x_dados(i);
        rhs_1_theta(1:3,i)=0;                                    
        rhs_1_theta(4,i)=0;
%         rhs_1_theta(5,i)=0;
        rhs_2_w(1:3,i)=0;                                        
        rhs_2_w(4,i)=0;
%         rhs_2_w(5,i)=0;
        rhs_2_theta(1:3,i)=dgdx(c,x_dados(i),sub_dominio(:));    
        rhs_2_theta(4,i)=0;
%         rhs_2_theta(5,i)=1;
        
        arhs_1_w(1:3,i)=0;                                      
        arhs_1_w(4,i)=0;  
%         arhs_1_w(5,i)=0;
        arhs_1_theta(1:3,i)=0;                                   
        arhs_1_theta(4,i)=0;
%         arhs_1_theta(5,i)=0;
        arhs_2_w(1:3,i)=0;                                       
        arhs_2_w(4,i)=0;
%         arhs_2_w(5,i)=0;
        arhs_2_theta(1:3,i)=0;                                   
        arhs_2_theta(4,i)=0;
%         arhs_2_theta(5,i)=0;
        
        rhs_2_phis(1:3,i)=H1*dgdx(c,x_dados(i),sub_dominio(:));     
        rhs_2_phis(4,i)=0;   %%%%%%   
%         rhs_2_phis(5,i)=H1;
        rhs_2_phia(1:3,i)=H2*dgdx(c,x_dados(i),sub_dominio(:));    
        rhs_2_phia(4,i)=0;   %%%%%%
%         rhs_2_phia(5,i)=H2;
        rhs_phiphis(1:3,i)=I_1*g(c,x_dados(i),sub_dominio(:));      
        rhs_phiphis(4,i)=I_1;   %%%%%%%%%%
%         rhs_phiphis(5,i)=I_1*x_dados(i);
        rhs_phiphia(1:3,i)=I_2*g(c,x_dados(i),sub_dominio(:));      
        rhs_phiphia(4,i)=I_2;   %%%%%%%%%%%
%         rhs_phiphia(5,i)=I_2*x_dados(i);
           
            case {'cc'}
        rhs_1_w(1:3,i)=g(c,x_dados(i),sub_dominio(:));                
        rhs_1_w(4,i)=1;    
%         rhs_1_w(5,i)=x_dados(i);  
        rhs_1_theta(1:3,i)=0;                                         
        rhs_1_theta(4,i)=0;
%         rhs_1_theta(5,i)=0;
        rhs_2_w(1:3,i)=0;                                              
        rhs_2_w(4,i)=0;
%         rhs_2_w(5,i)=0;
        rhs_2_theta(1:3,i)=g(c,x_dados(i),sub_dominio(:));              
        rhs_2_theta(4,i)=1;
%         rhs_2_theta(5,i)=x_dados(i);
        
        arhs_1_w(1:3,i)=0;                                             
        arhs_1_w(4,i)=0;
%         arhs_1_w(5,i)=0;
        arhs_1_theta(1:3,i)=0;                                         
        arhs_1_theta(4,i)=0;
%         arhs_1_theta(5,i)=0;
        arhs_2_w(1:3,i)=0;                                             
        arhs_2_w(4,i)=0;
%         arhs_2_w(5,i)=0;
        arhs_2_theta(1:3,i)=0;                                         
        arhs_2_theta(4,i)=0;
%         arhs_2_theta(5,i)=0;
        
        rhs_2_phis(1:3,i)=H1*dgdx(c,x_dados(i),sub_dominio(:));      
        rhs_2_phis(4,i)=0; %%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phis(5,i)=H1;
        rhs_2_phia(1:3,i)=H2*dgdx(c,x_dados(i),sub_dominio(:));     
        rhs_2_phia(4,i)=0;%%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phia(5,i)=H2;
        rhs_phiphis(1:3,i)=I_1*g(c,x_dados(i),sub_dominio(:));       
        rhs_phiphis(4,i)=I_1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphis(5,i)=I_1*x_dados(i);
        rhs_phiphia(1:3,i)=I_2*g(c,x_dados(i),sub_dominio(:));      
        rhs_phiphia(4,i)=I_2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphia(5,i)=I_2*x_dados(i);
            
            case {'cl'}
        rhs_1_w(1:3,i)=g(c,x_dados(i),sub_dominio(:));            
        rhs_1_w(4,i)=1;
%         rhs_1_w(5,i)=x_dados(i);
        rhs_1_theta(1:3,i)=0;                                     
        rhs_1_theta(4,i)=0;
%         rhs_1_theta(5,i)=0;
        rhs_2_w(1:3,i)=0;                                         
        rhs_2_w(4,i)=0;
%         rhs_2_w(5,i)=0;
        rhs_2_theta(1:3,i)=g(c,x_dados(i),sub_dominio(:));        
        rhs_2_theta(4,i)=1;
%         rhs_2_theta(5,i)=x_dados(i);
        
        arhs_1_w(1:3,i)=0;                                        
        arhs_1_w(4,i)=0;
%         arhs_1_w(5,i)=0;
        arhs_1_theta(1:3,i)=0;                                    
        arhs_1_theta(4,i)=0;
%         arhs_1_theta(5,i)=0;
        arhs_2_w(1:3,i)=0;                                        
        arhs_2_w(4,i)=0;
%         arhs_2_w(5,i)=0;
        arhs_2_theta(1:3,i)=0;                                    
        arhs_2_theta(4,i)=0;
%         arhs_2_theta(5,i)=0;
        
        rhs_2_phis(1:3,i)=H1*dgdx(c,x_dados(i),sub_dominio(:));   
        rhs_2_phis(4,i)=0;   %%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phis(5,i)=H1;
        rhs_2_phia(1:3,i)=H2*dgdx(c,x_dados(i),sub_dominio(:));   
        rhs_2_phia(4,i)=0;    %%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phia(5,i)=H2;
        rhs_phiphis(1:3,i)=I_1*g(c,x_dados(i),sub_dominio(:));    
        rhs_phiphis(4,i)=I_1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphis(5,i)=I_1*x_dados(i);
        rhs_phiphia(1:3,i)=I_2*g(c,x_dados(i),sub_dominio(:));    
        rhs_phiphia(4,i)=I_2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphia(5,i)=I_2*x_dados(i);
        end

[Axi,Axj]=meshgrid(sub_dominio);
matriz_pesos=[g(c,Axi,Axj) ones(3,1);    %!!!!!!!!!!!!!!!!!!!!
              ones(1,3) 0];
       
% matriz_pesos=[g(c,Axi,Axj) ones(3,1) [x_dados(i);x_dados(i+1);x_dados(i+2)];    %!!!!!!!!!!!!!!!!!!!!
%               ones(1,3) 0 0;
%               [x_dados(i) x_dados(i+1) x_dados(i+2)] 0 0]; 

pesos_1_w(i,:)=matriz_pesos\rhs_1_w(:,i);
pesos_1_theta(i,:)=matriz_pesos\rhs_1_theta(:,i);
pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);

apesos_1_w(i,:)=matriz_pesos\arhs_1_w(:,i);
apesos_1_theta(i,:)=matriz_pesos\arhs_1_theta(:,i);
apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);

pesos_2_phis(i,:)=matriz_pesos\rhs_2_phis(:,i);   %%%%%%%%%%%%%%%%%%%%%%%%
pesos_2_phia(i,:)=matriz_pesos\rhs_2_phia(:,i);   %%%%%%%%%%%%%%%%%%%%%%%%%%%
pesos_phiphis(i,:)=matriz_pesos\rhs_phiphis(:,i);  %%%%%%%%%%%%%%%%%%%%
pesos_phiphia(i,:)=matriz_pesos\rhs_phiphia(:,i);         %%%%%%%%%%%%%%%%%%%%%%%%%
        
    elseif  i==numel(x_dados)
        sub_dominio=[x_dados(i), x_dados(i-1), x_dados(i-2)];
         switch cfstr
            case {'ss'}
        rhs_1_w(1:3,i)=g(c,x_dados(i),sub_dominio(:));            
        rhs_1_w(4,i)=1;
%         rhs_1_w(5,i)=x_dados(i);
        rhs_1_theta(1:3,i)=0;                                     
        rhs_1_theta(4,i)=0;
%         rhs_1_theta(5,i)=0;
        rhs_2_w(1:3,i)=0;                                         
        rhs_2_w(4,i)=0;
%         rhs_2_w(5,i)=0;
        rhs_2_theta(1:3,i)=dgdx(c,x_dados(i),sub_dominio(:));     
        rhs_2_theta(4,i)=0;
%         rhs_2_theta(5,i)=1;
        
        arhs_1_w(1:3,i)=0;                                        
        arhs_1_w(4,i)=0;
%         arhs_1_w(5,i)=0;
        arhs_1_theta(1:3,i)=0;                                    
        arhs_1_theta(4,i)=0;
%         arhs_1_theta(5,i)=0;
        arhs_2_w(1:3,i)=0;                                        
        arhs_2_w(4,i)=0;
%         arhs_2_w(5,i)=0;
        arhs_2_theta(1:3,i)=0;                                    
        arhs_2_theta(4,i)=0;
%         arhs_2_theta(5,i)=0;
        
        rhs_2_phis(1:3,i)=H1*dgdx(c,x_dados(i),sub_dominio(:));    
        rhs_2_phis(4,i)=0;   %%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phis(5,i)=H1;
        rhs_2_phia(1:3,i)=H2*dgdx(c,x_dados(i),sub_dominio(:));    
        rhs_2_phia(4,i)=0;     %%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phia(5,i)=H2;
        rhs_phiphis(1:3,i)=I_1*g(c,x_dados(i),sub_dominio(:));     
        rhs_phiphis(4,i)=I_1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphis(5,i)=I_1*x_dados(i);
        rhs_phiphia(1:3,i)=I_2*g(c,x_dados(i),sub_dominio(:));     
        rhs_phiphia(4,i)=I_2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphia(5,i)=I_2*x_dados(i);
        
             case {'cc'}
        rhs_1_w(1:3,i)=g(c,x_dados(i),sub_dominio(:));            
        rhs_1_w(4,i)=1;
%         rhs_1_w(5,i)=x_dados(i);
        rhs_1_theta(1:3,i)=0;                                      
        rhs_1_theta(4,i)=0;
%         rhs_1_theta(5,i)=0;
        rhs_2_w(1:3,i)=0;                                         
        rhs_2_w(4,i)=0;
%         rhs_2_w(5,i)=0;
        rhs_2_theta(1:3,i)=g(c,x_dados(i),sub_dominio(:));        
        rhs_2_theta(4,i)=1;
%         rhs_2_theta(5,i)=x_dados(i);
        
        arhs_1_w(1:3,i)=0;                                      
        arhs_1_w(4,i)=0;
%         arhs_1_w(5,i)=0;
        arhs_1_theta(1:3,i)=0;                                  
        arhs_1_theta(4,i)=0;
%         arhs_1_theta(5,i)=0;
        arhs_2_w(1:3,i)=0;                                      
        arhs_2_w(4,i)=0;
%         arhs_2_w(5,i)=0;
        arhs_2_theta(1:3,i)=0;                                  
        arhs_2_theta(4,i)=0;
%         arhs_2_theta(5,i)=0;
        
        rhs_2_phis(1:3,i)=H1*dgdx(c,x_dados(i),sub_dominio(:));   
        rhs_2_phis(4,i)=0;   %%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phis(5,i)=H1;
        rhs_2_phia(1:3,i)=H2*dgdx(c,x_dados(i),sub_dominio(:));  
        rhs_2_phia(4,i)=0;  %%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phia(5,i)=H2;
        rhs_phiphis(1:3,i)=I_1*g(c,x_dados(i),sub_dominio(:));    
        rhs_phiphis(4,i)=I_1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphis(5,i)=I_1*x_dados(i);
        rhs_phiphia(1:3,i)=I_2*g(c,x_dados(i),sub_dominio(:));   
        rhs_phiphia(4,i)=I_2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphia(5,i)=I_2*x_dados(i);
       
             case {'cl'}
        rhs_1_w(1:3,i)=dgdx(c,x_dados(i),sub_dominio(:));        
        rhs_1_w(4,i)=0;
%         rhs_1_w(5,i)=1;
        rhs_1_theta(1:3,i)=g(c,x_dados(i),sub_dominio(:));      
        rhs_1_theta(4,i)=1;
%         rhs_1_theta(5,i)=x_dados(i);
        rhs_2_w(1:3,i)=0;                                        
        rhs_2_w(4,i)=0;
%         rhs_2_w(5,i)=0;
        rhs_2_theta(1:3,i)=dgdx(c,x_dados(i),sub_dominio(:));    
        rhs_2_theta(4,i)=0;
%         rhs_2_theta(5,i)=1; 
        
        arhs_1_w(1:3,i)=0;                       
        arhs_1_w(4,i)=0;
%         arhs_1_w(5,i)=0;
        arhs_1_theta(1:3,i)=0;                    
        arhs_1_theta(4,i)=0;
%         arhs_1_theta(5,i)=0;
        arhs_2_w(1:3,i)=0;                        
        arhs_2_w(4,i)=0;
%         arhs_2_w(5,i)=0;
        arhs_2_theta(1:3,i)=0;                   
        arhs_2_theta(4,i)=0;
%         arhs_2_theta(5,i)=0;
        
        rhs_2_phis(1:3,i)=H1*dgdx(c,x_dados(i),sub_dominio(:));  
        rhs_2_phis(4,i)=0;   %%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phis(5,i)=H1;
        rhs_2_phia(1:3,i)=H2*dgdx(c,x_dados(i),sub_dominio(:));  
        rhs_2_phia(4,i)=0;  %%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phia(5,i)=H2;
        rhs_phiphis(1:3,i)=I_1*g(c,x_dados(i),sub_dominio(:));   
        rhs_phiphis(4,i)=I_1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphis(5,i)=I_1*x_dados(i);
        rhs_phiphia(1:3,i)=I_2*g(c,x_dados(i),sub_dominio(:));   
        rhs_phiphia(4,i)=I_2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphia(5,i)=I_2*x_dados(i);
         end
[Axi,Axj]=meshgrid(sub_dominio);
matriz_pesos=[g(c,Axi,Axj) ones(3,1);    %!!!!!!!!!!!!!!!!!!!!
              ones(1,3) 0];
% matriz_pesos=[g(c,Axi,Axj) ones(3,1) [x_dados(i);x_dados(i-1);x_dados(i-2)];    %!!!!!!!!!!!!!!!!!!!!
%               ones(1,3) 0 0;
%               [x_dados(i) x_dados(i-1) x_dados(i-2)] 0 0]; 
          
pesos_1_w(i,:)=matriz_pesos\rhs_1_w(:,i);
pesos_1_theta(i,:)=matriz_pesos\rhs_1_theta(:,i);
pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);

apesos_1_w(i,:)=matriz_pesos\arhs_1_w(:,i);
apesos_1_theta(i,:)=matriz_pesos\arhs_1_theta(:,i);
apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i);   

pesos_2_phis(i,:)=matriz_pesos\rhs_2_phis(:,i);   %%%%%%%%%%%%%%%%%%%%%%%%
pesos_2_phia(i,:)=matriz_pesos\rhs_2_phia(:,i);   %%%%%%%%%%%%%%%%%%%%%%%%%%%
pesos_phiphis(i,:)=matriz_pesos\rhs_phiphis(:,i);  %%%%%%%%%%%%%%%%%%%%
pesos_phiphia(i,:)=matriz_pesos\rhs_phiphia(:,i);         %%%%%%%%%%%%%%%%%%%%%%%%%
else
         
        sub_dominio=[x_dados(i), x_dados(i-1), x_dados(i+1)];
        
        rhs_1_w(1:3,i)=B55*d2gdx2(c,x_dados(i),sub_dominio(:));                                             
        rhs_1_w(4,i)=0;
%         rhs_1_w(5,i)=0;
        rhs_1_theta(1:3,i)=B55*dgdx(c,x_dados(i),sub_dominio(:));                                          
        rhs_1_theta(4,i)=0;
%         rhs_1_theta(5,i)=B55;
        rhs_2_w(1:3,i)=-B55*dgdx(c,x_dados(i),sub_dominio(:));                                             
        rhs_2_w(4,i)=0;
%         rhs_2_w(5,i)=-B55;
        rhs_2_theta(1:3,i)=(D11+G)*d2gdx2(c,x_dados(i),sub_dominio(:))-B55*g(c,x_dados(i),sub_dominio(:)); 
        rhs_2_theta(4,i)=-B55;
%         rhs_2_theta(5,i)=-B55*x_dados(i);
        
        arhs_1_w(1:3,i)=-J0*g(c,x_dados(i),sub_dominio(:));      
        arhs_1_w(4,i)=-J0;
%         arhs_1_w(5,i)=-J0*x_dados(i);
        arhs_1_theta(1:3,i)=0;                                   
        arhs_1_theta(4,i)=0;
%         arhs_1_theta(5,i)=0;
        arhs_2_w(1:3,i)=0;                                        
        arhs_2_w(4,i)=0;
%         arhs_2_w(5,i)=0;
        arhs_2_theta(1:3,i)=-J2*g(c,x_dados(i),sub_dominio(:));  
        arhs_2_theta(4,i)=-J2;
%         arhs_2_theta(5,i)=-J2*x_dados(i);
        
        rhs_2_phis(1:3,i)=H1*dgdx(c,x_dados(i),sub_dominio(:));   
        rhs_2_phis(4,i)=0;  %%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phis(5,i)=H1;
        rhs_2_phia(1:3,i)=H2*dgdx(c,x_dados(i),sub_dominio(:));   
        rhs_2_phia(4,i)=0;  %%%%%%%%%%%%%%%%%%%%%%
%         rhs_2_phia(5,i)=H2;
        rhs_phiphis(1:3,i)=I_1*g(c,x_dados(i),sub_dominio(:));    
        rhs_phiphis(4,i)=I_1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphis(5,i)=I_1*x_dados(i);
        rhs_phiphia(1:3,i)=I_2*g(c,x_dados(i),sub_dominio(:));    
        rhs_phiphia(4,i)=I_2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rhs_phiphia(5,i)=I_2*x_dados(i);
        
[Axi,Axj]=meshgrid(sub_dominio);
matriz_pesos=[g(c,Axi,Axj) ones(3,1);    %!!!!!!!!!!!!!!!!!!!!
              ones(1,3) 0];
% matriz_pesos=[g(c,Axi,Axj) ones(3,1) [x_dados(i);x_dados(i-1);x_dados(i+1)];    %!!!!!!!!!!!!!!!!!!!!
%               ones(1,3) 0 0;
%               [x_dados(i) x_dados(i-1) x_dados(i+1)] 0 0]; 

pesos_1_w(i,:)=matriz_pesos\rhs_1_w(:,i);
pesos_1_theta(i,:)=matriz_pesos\rhs_1_theta(:,i);
pesos_2_w(i,:)=matriz_pesos\rhs_2_w(:,i);
pesos_2_theta(i,:)=matriz_pesos\rhs_2_theta(:,i);

apesos_1_w(i,:)=matriz_pesos\arhs_1_w(:,i);
apesos_1_theta(i,:)=matriz_pesos\arhs_1_theta(:,i);
apesos_2_w(i,:)=matriz_pesos\arhs_2_w(:,i);
apesos_2_theta(i,:)=matriz_pesos\arhs_2_theta(:,i); 

pesos_2_phis(i,:)=matriz_pesos\rhs_2_phis(:,i);   %%%%%%%%%%%%%%%%%%%%%%%%
pesos_2_phia(i,:)=matriz_pesos\rhs_2_phia(:,i);   %%%%%%%%%%%%%%%%%%%%%%%%%%%
pesos_phiphis(i,:)=matriz_pesos\rhs_phiphis(:,i);  %%%%%%%%%%%%%%%%%%%%
pesos_phiphia(i,:)=matriz_pesos\rhs_phiphia(:,i);         %%%%%%%%%%%%%%%%%%%%%%%%%
end
end

%%
%   ASSEMBLAGEM DOS PESOS NAS MATRIZES DE INERCIA E DE RIGIDEZ

L_total=zeros(2*n,2*n);
A_total=zeros(2*n,2*n);

K_tphis=zeros(n,n);   %%%%%%%%%%%%%%%%%%%%%%
K_tphia=zeros(n,n);   %%%%%%%%%%%%%%%%%%%%%%
K_phiphis=zeros(n,n);  %%%%%%%%%%%%%%%%%%%%%%%%
K_phiphia=zeros(n,n);  %%%%%%%%%%%%%%%%%%%%%%%%%%%

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

K_tphis(i,j-1)=pesos_2_phis(i,2);  %%%%%%%%%%%%%%%%%%%%
K_tphis(i,j)=pesos_2_phis(i,1);    %%%%%%%%%%%%%%%%%%%
K_tphis(i,j+1)=pesos_2_phis(i,3);  %%%%%%%%%%%%%%%%%%%%

K_tphia(i,j-1)=pesos_2_phia(i,2);   %%%%%%%%%%%%%%%%%%%%%%%%
K_tphia(i,j)=pesos_2_phia(i,1);    %%%%%%%%%%%%%%%%%%%%%%%%%
K_tphia(i,j+1)=pesos_2_phia(i,3);  %%%%%%%%%%%%%%%%%%%%%%%%%%

K_phiphis(i,j-1)=pesos_phiphis(i,2);   %%%%%%%%%%%%%%%%%%%%%%%%
K_phiphis(i,j)=pesos_phiphis(i,1);    %%%%%%%%%%%%%%%%%%%%%%%%%
K_phiphis(i,j+1)=pesos_phiphis(i,3);  %%%%%%%%%%%%%%%%%%%%%%%%%%

K_phiphia(i,j-1)=pesos_phiphia(i,2);   %%%%%%%%%%%%%%%%%%%%%%%%
K_phiphia(i,j)=pesos_phiphia(i,1);    %%%%%%%%%%%%%%%%%%%%%%%%%
K_phiphia(i,j+1)=pesos_phiphia(i,3);  %%%%%%%%%%%%%%%%%%%%%%%%%%

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
K_tphis(1,1:3)=pesos_2_phis(1,1:3);    %%%%%%%%%%%%%%%%%%%%%%%%
K_tphis(n,n-2:n)=pesos_2_phis(n,3:-1:1); %%%%%%%%%%%%%%%%%%%%%

K_tphia(1,1:3)=pesos_2_phia(1,1:3);   %%%%%%%%%%%%%%%%%%%%%%%%%%
K_tphia(n,n-2:n)=pesos_2_phia(n,3:-1:1);  %%%%%%%%%%%%%%%%%%%%%%%%%

K_phiphis(1,1:3)=pesos_phiphis(1,1:3);   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_phiphis(n,n-2:n)=pesos_phiphis(n,3:-1:1);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
K_phiphia(1,1:3)=pesos_phiphia(1,1:3);   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_phiphia(n,n-2:n)=pesos_phiphia(n,3:-1:1);   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% L_TOTAL FINAL (PODE SER SMART, CLOSE CIRCUIT)
L_total(n+1:2*n,n+1:2*n)=L_total(n+1:2*n,n+1:2*n)+K_tphis*(K_phiphis^-1)*K_tphis;   %SMART (1 SENSOR E 1 ATUADOR)
% L_total(n+1:2*n,n+1:2*n)=L_total(n+1:2*n,n+1:2*n);                                  %CLOSED-CIRCUIT 0V


% TERMOS DE FRONTEIRA
L_total(1,1:end)=0; L_total(n,1:end)=0; L_total(n+1,1:end)=0; L_total(end,1:end)=0;

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
lambda_mode_w(1:n,p)=lambda_vec(1:n,p);
lambda_mode_phi_x(1:n,p)=lambda_vec(n+1:end,p);
lambda_mode=[lambda_mode_w;lambda_mode_phi_x];

figure(1)
subplot(1,3,1);plot(x_dados, lambda_mode_w(:,p));hold on;title(['w(' num2str(m) ')_{exact} = ' num2str(sol_exacta/(2*pi),'%6.4f')]);legend(['w(' num2str(m) ')=' num2str(sqrt(lambda(p))/(2*pi),'%6.4f')])
subplot(1,3,2);plot(x_dados, lambda_mode_phi_x(:,p));hold on;title(['w(' num2str(m) ')_{exact} = ' num2str(sol_exacta/(2*pi),'%6.4f')]);legend(['w(' num2str(m) ')=' num2str(sqrt(lambda(p))/(2*pi),'%6.4f')])

%%
%----NEWMARK  (TEM DE SER COM SMART)
% Gv=0.0001;
Gv=0;

C_tt=K_tphia*(K_phiphis^-1)*K_tphis;
C_total=zeros(2*n,2*n);
C_total(n+1:2*n,n+1:2*n)=C_tt;
C_total=-Gv*C_total;
C_total(n+1,n+1:end)=0;
C_total(end,n+1:end)=0;

%cond. iniciais 
vetor_carga=zeros(2*n,1);
vetor_carga(2:n-1)=carga;
solucao_estatica=L_total\vetor_carga;
x_0=solucao_estatica; v_0=zeros(2*n,1);
vetor_f=zeros(2*n,1);   %caso se queira impôr força inicial ao inves de deslocamento inicial
a_0=pinv(A_total)*(vetor_f-C_total*v_0-L_total*x_0);
delta=1/2; alpha=1/4;
% delta_t=1/(freq*200);   %delta t
delta_t=1/100000;
a0=1/(alpha*delta_t^2); a1=delta/(alpha*delta_t); a2=1/(alpha*delta_t); a3=1/(2*alpha)-1;
a4=delta/alpha-1; a5=(delta_t/2)*(delta/alpha-2); a6=delta_t*(1-delta); a7=delta*delta_t;

%rigidez efetiva
K_efe=L_total+a0*A_total+a1*C_total;

% t_final=10/freq;   % t final
t_final=0.02;
n_t=int64(t_final/delta_t+1);
t=zeros(n_t,1);
x_t=zeros(2*n,n_t); x_t(:,1)=x_0;
v_t=zeros(2*n,n_t); v_t(:,1)=v_0;
a_t=zeros(2*n,n_t); a_t(:,1)=a_0;
for i=2:n_t
  t(i)=t(i-1)+delta_t;
  F_efe=A_total*(a0*x_t(:,i-1)+a2*v_t(:,i-1)+a3*a_t(:,i-1))+C_total*(a1*x_t(:,i-1)+a4*v_t(:,i-1)+a5*a_t(:,i-1));
  x_t(:,i)=K_efe\F_efe;
  a_t(:,i)=a0*(x_t(:,i)-x_t(:,i-1))-a2*v_t(:,i-1)-a3*a_t(:,i-1);
  v_t(:,i)=v_t(:,i-1)+a6*a_t(:,i-1)+a7*a_t(:,i);
end
switch cfstr
    case {'cc'}
      x_max=x_t(ceil(n/2),:);  
    case {'ss'}
      x_max=x_t(ceil(n/2),:);  
    case{'cl'}  
      x_max=x_t(n,:);
end

freq=sqrt(lambda(p))/(2*pi);
yharm=x_0(n)*cos(2*pi*freq*t);      %COM A AMPLITUDE IGUAL AO CASO CL
figure(2)
plot(t,yharm);
hold on
plot(t,x_max);      
hold on
%%
X=fft(x_max);
X_mag=abs(X(1:ceil(n_t/2)));
[pk_vals, pk_locs]=findpeaks(X_mag);
%remove peaks below threshold
% inds=find(X_mag(pk_locs)<1);
% pk_locs(inds)=[];

%determine frequencies
pk_freqs=zeros(length(pk_locs),1);
for i=1:length(pk_locs)
pk_freqs(i)=(pk_locs(i)-1)/t_final;
end

figure(3)
plot(X_mag);
hold on
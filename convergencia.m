clear all
nmax=100;
p=3;   %n� de polinomios
mode=2; %modo natural
c=1;  %shape parameter c=0.5 for p=0, (c=1 for p=1 mau), (c=4 p=2 mt mau), c=1:20 for p=3
y=zeros(1,nmax-9); 
% y_exacta=zeros(1,nmax-9);
freq_ansys=[491.3 2929 7658];  %frequencias para bimorph (2 camadas)
for i=10:nmax
    c=2/sqrt(i);
  [freq]=conv(i,p,c,mode); 
  y(i-9)=freq;
%   y_exacta(i-9)=freq_exacta;
  disp(i/nmax*100);
end
y_ex_cl=freq_ansys(mode)*ones(1,nmax-9);
plot(10:nmax,y);
hold on
% plot(10:nmax,y_exacta); %para SS sem piezoelectricidade
plot(10:nmax,y_ex_cl); %para CL open-circuit ansys
hold on
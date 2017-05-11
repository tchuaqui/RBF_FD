clear all
nmax=150;
p=0;   %nº de polinomios
c=5;  %shape parameter
y=zeros(1,nmax-9); y_exacta=zeros(1,nmax-9);
for i=5:nmax
  [freq,freq_exacta]=conv(i,p,c); 
  y(i-9)=freq;
  y_exacta(i-9)=freq_exacta;
  disp(i/nmax*100);
end
y_ex_cl=464.3*ones(1,nmax-9);
plot(10:nmax,y);
hold on
% plot(10:nmax,y_exacta); %para SS sem piezoelectricidade
plot(10:nmax,y_ex_cl); %para CL open-circuit ansys
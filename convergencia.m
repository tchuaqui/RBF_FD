clear all
nmax=100;
p=3;   %nº de polinomios
c=2;  %shape parameter c=0.5 for p=0, (c=1 for p=1 mau), (c=4 p=2 mt mau), c=1:20 for p=3
y=zeros(1,nmax-9); y_exacta=zeros(1,nmax-9);
for i=10:nmax
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
hold on
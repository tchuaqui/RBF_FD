clear all
nmax=100;
p=3;   %n� de polinomios
c=1;  %shape parameter c=0.5 for p=0, (c=1 for p=1 mau), (c=4 p=2 mt mau), c=1:20 for p=3
mode=3; %1st, 2nd or 3rd mode
y=zeros(1,nmax-9); y_exacta=zeros(1,nmax-9);
for i=10:nmax
  [freq,freq_exacta]=conv(i,p,c,mode); 
  y(i-9)=freq;
  y_exacta(i-9)=freq_exacta;
  disp(i/nmax*100);
end
f_mode=[464.3 2783 7093];
y_ex_cl=f_mode(mode)*ones(1,nmax-9);
plot(10:nmax,y);
hold on
% plot(10:nmax,y_exacta); %para SS sem piezoelectricidade
plot(10:nmax,y_ex_cl); %para CL open-circuit ansys
hold on
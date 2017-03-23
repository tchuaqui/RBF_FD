nmax=150;
y=zeros(1,nmax-9); y_exacta=zeros(1,nmax-9);
for i=10:nmax
  [freq,freq_exacta]=conv(i); 
  y(i-9)=freq;
  y_exacta(i-9)=freq_exacta;
  disp(i/nmax*100);
end
y_ex_cl=453*ones(1,nmax-9);
plot(10:nmax,y);
hold on
% % plot(10:nmax,y_exacta); %para SS
plot(10:nmax,y_ex_cl); %para CL
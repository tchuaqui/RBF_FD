
clear all

syms c xi xj  
% 
% rbf=sqrt((xi-xj)^2+c^2); %Multiquadrica   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rbf=sqrt((1/c^2)*(xi-xj)^2+1);
pol=sym('1');
% rbf=((xi-xj)^2+c^2)^(-1/2); %Multiquadrica Inversa
% rbf=sqrt((xi-xj)^2+c^2)^2*log(sqrt((xi-xj)^2+c^2));%TPS 
% rbf=(sqrt((xi-xj)^2+c^2))^(7);%piece-wise polynomial, c impar
% rbf=exp(-3^2*sqrt(((xi-xj)^2+c^2))^2);%gauss

g=vectorize(char(rbf));
dgdx_s=diff(rbf,xj);dgdx=vectorize(char(dgdx_s));
d2gdx2_s=diff(diff(rbf,xj),xj);d2gdx2=vectorize(char(d2gdx2_s));


p=vectorize(char(pol));
dpdx_s=diff(pol,xj);dpdx=vectorize(char(dpdx_s));
d2pdx2_s=diff(diff(pol,xj),xj);d2pdx2=vectorize(char(d2pdx2_s));


%----create g.m--------
f=fopen('g.m','w');
fprintf(f,'function y=g(c,xi,xj)\n');
fprintf(f,'y=');fprintf(f,g);fprintf(f,';');
status=fclose(f);
%----create dgdx.m--------
f=fopen('dgdx.m','w');
fprintf(f,'function y=dgdx(c,xi,xj)\n');
fprintf(f,'y=');fprintf(f,dgdx);fprintf(f,';');
status=fclose(f);
%----create d2gdx2.m--------
f=fopen('d2gdx2.m','w');
fprintf(f,'function y=d2gdx2(c,xi,xj)\n');
fprintf(f,'y=');fprintf(f,d2gdx2);fprintf(f,';');
status=fclose(f);
%%
%----create p.m--------
f=fopen('pol.m','w');
fprintf(f,'function y=p(c,xi,xj)\n');
fprintf(f,'y=');fprintf(f,p);fprintf(f,';');
status=fclose(f);
%----create dpdx.m--------
f=fopen('dpoldx.m','w');
fprintf(f,'function y=dpdx(c,xi,xj)\n');
fprintf(f,'y=');fprintf(f,dpdx);fprintf(f,';');
status=fclose(f);
%----create d2pdx2.m--------
f=fopen('d2poldx2.m','w');
fprintf(f,'function y=d2pdx2(c,xi,xj)\n');
fprintf(f,'y=');fprintf(f,d2pdx2);fprintf(f,';');
status=fclose(f);


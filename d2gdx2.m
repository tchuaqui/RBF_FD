function y=d2gdx2(c,xi,xj)
y=1./((xi - xj).^2 + c.^2).^(1./2) - (2.*xi - 2.*xj).^2./(4.*((xi - xj).^2 + c.^2).^(3./2));
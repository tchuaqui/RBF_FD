function y=dgdx(c,xi,xj)
y=-(2.*xi - 2.*xj)./(2.*c.^2.*((xi - xj).^2./c.^2 + 1).^(1./2));
function temp=f(x)
  %if(x<1)
  %temp=-2*pi*x.*cos(pi*x)+(pi^2*(1.^x+x.^2)).*sin(pi*x);
  %else
  temp=-2*pi*x.*cos(pi*x)+(pi^2*(1.^x+x.^2)).*sin(pi*x)+q(x)*sin(pi*x);
  %endif
endfunction

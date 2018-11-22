function x=mod_thomas(a,d,r,n)
  % symetry
  a=[a 0];
  temp=a(1);
  a(1)=a(1)/d(1);
  r(1)=r(1)/d(1);
  for i=2:n
    denom=d(i)-temp*a(i-1);
    r(i)=(r(i)-temp*r(i-1))/denom;
    temp=a(i);
    a(i)=a(i)/denom;    
  endfor
  x(n)=r(n);
  for i=n-1:-1:1
    x(i)=r(i)-a(i)*x(i+1);    
  endfor
endfunction

function value_at_point=valueAtPoint(c,h,x,X,a,n)
  if(x<X(2))
  value_at_point=1/h*c(1)*(x-X(1));
  elseif (x>X(n+1))
    value_at_point=1/h*c(n)*(X(n+2)-x);
    else
      temp=floor((x-a)/h);
      value_at_point=(c(temp+1)*(x-X(temp+1))+c(temp)*(X(temp+2)-x))/h;
    endif
endfunction

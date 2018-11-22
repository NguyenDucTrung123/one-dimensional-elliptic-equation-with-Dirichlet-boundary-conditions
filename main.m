function main(n,a,b,x)
  %gio tam coi co bien 0,1
  %n la so diem chia khac a,b( a=x_0;b=x_{n+1})
  %=> co (n+1) doan do dai h
  
  tic
  format long
  F=zeros(n,1);
  h=(b-a)/(n+1);
  X=[a:h:b];
  
  %evaluate approximate itegral by 2-point Gauss Quadrature on each interval [x_i,x_(i+1)]
  _X=(X(1:n+1)+X(2:n+2))./2; %size: 1*(n+1)
  _Xp=_X+h/(2*sqrt(3))*ones(1,n+1);
  _Xm=_X-h/(2*sqrt(3))*ones(1,n+1);
  
  %function p is applied to Gauss points 
  %{
  p_p=feval('p',_Xp);
  p_m=feval('p',_Xm);
  %}
  for i=1:n+1
    p_p(i)=p(_Xp(i));
    p_m(i)=p(_Xm(i));
  endfor
  
  p_pm=p_p+p_m;
 
 %function q is applied to Gauss points
  %q_p=feval('q',_Xp);
  %q_m=feval('q',_Xm);
  for i=1:n+1
    q_p(i)=q(_Xp(i));
    q_m(i)=q(_Xm(i));
  endfor
  
  q_pm0=q_m+q_p;
  q_pm1=(4/3-2/sqrt(3))*q_m+(4/3+2/sqrt(3))*q_p;
  q_pm2=(4/3+2/sqrt(3))*q_m+(4/3-2/sqrt(3))*q_p;
  
  %diag(B):
  diag_B=1/(2*h)*(p_pm(1:n)+p_pm(2:n+1))+h/8*(q_pm1(1:n)+q_pm2(2:n+1));
  diag_B_1=h/12*q_pm0(2:n)-1/(2*h)*p_pm(2:n);
  
  %function f is applied to Gauss points
  %f_p=feval('f',_Xp);
  %f_m=feval('f',_Xm);
  for i=1:n+1
    f_p(i)=f(_Xp(i));
    f_m(i)=f(_Xm(i));
  endfor
  
  f_pm1=(1-1/sqrt(3))*f_m+(1+1/sqrt(3))*f_p;
  f_pm2=(1+1/sqrt(3))*f_m+(1-1/sqrt(3))*f_p;
  
  %F:
  F=h/4*(f_pm1(1:n)+f_pm2(2:n+1));
  %Find c by three different methods:
  %tic;
  %c=thomas(diag_B_1,diag_B,diag_B_1,F,n);
  %t=toc;
  %fprintf("2)time for thomas to solve the system of linear equations: %f \n",t);
  %tic;
  c=mod_thomas(diag_B_1,diag_B,F,n);
  %(in case use mod_thomas to solve the system of linear equations)
  %t=toc;
  %fprintf("2)time for mod_thomas to solve the system of linear equations: %f \n",t);
  %fprintf("2)time for tridiag to solve the system of linear equations: %f \n",t);
  %approximate function is applied to a particular x in [a,b]:
  t1=toc;
  if(x<X(2))
  value_at_point=1/h*c(1)*(x-X(1));
  elseif (x>X(n+1))
    value_at_point=1/h*c(n)*(X(n+2)-x);
    else
      temp=floor((x-a)/h);
      value_at_point=(c(temp+1)*(x-X(temp+1))+c(temp)*(X(temp+2)-x))/h;
    endif
    fprintf("1)The number of finite elements: %d\n",n+1);
    fprintf("2)Time : %f \n",t1)
    fprintf("3)Approximate value of function when applied to x=%f is: %f\n ",x,value_at_point);
    %point error: abs(exSol(x)-value_at_point);
    % Now, we evaluate error between approximate and exact solution:
    % approximate solution is applied to Gauss points:
    xw=MonoGaussPoints(64);
    A=(b-a)/2;
    B=(b+a)/2;
    %abscissa = A*xw(:,2)+(b+a)/2*ones(64,1);
    s=0;
    for i=1:64
      temp=A*xw(i,2)+B;
      s=s+xw(i,1)*(exSol(temp)-valueAtPoint(c,h,temp,X,a,n))^2;
    endfor
    fprintf("4)L2-error between exact solution and approximation is: %.8g\n",sqrt(s));     
    
    %{
    c_p=1/2*(1+1/sqrt(3))*c;
    c_m=1/2*(1-1/sqrt(3))*c;
    U_h_x=zeros(1,n+1);
    U_h_xx=zeros(1,n+1);
    U_h_x(1)=c_m(1);
    U_h_x(n+1)=c_p(n);
    U_h_x(2:n)=c_p(1:n-1)+c_m(2:n);
    U_h_xx(1)=c_p(1);
    U_h_xx(n+1)=c_m(n);
    U_h_xx(2:n)=c_m(1:n-1)+c_p(2:n);
    % exact solution is applied to Gauss points:
    U_x=exSol(_Xm);
    U_xx=exSol(_Xp);
    % approximate L2-norm of (exSol-approximateSol):
    E=h/2*(sum((U_x-U_h_x).^2)+sum(U_xx-U_h_xx).^2);
    fprintf("4)L2-error between exact solution and approximation is: %.8g\n",sqrt(E));
%}
endfunction

function [residual, g1, g2] = nobond_noadj_linear_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                    columns: equations in order of declaration
%                                                    rows: variables in declaration order
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: equations in order of declaration
%                                                       rows: variables in declaration order
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 11, 1);

%
% Model equations
%

T49 = params(2)*params(17)/(1-params(14))^2;
T103 = (-(params(17)^(-2)))*(1-params(4)+params(13)*params(3)*(1-params(16))/params(15));
T114 = params(13)*params(3)*(1-params(16))/(params(17)*params(15)^2);
lhs =params(14)*y(1)+params(18)*y(3)+y(4);
rhs =y(2);
residual(1)= lhs-rhs;
lhs =y(1);
rhs =(-(1-params(3)))*params(13)/params(14)*y(6)+(1-params(3))*(1-params(16))/params(14)*y(7)-y(3)*params(13)*(1-params(3))*(1-params(16))*params(14)^(-2);
residual(2)= lhs-rhs;
lhs =y(1);
rhs =y(2)*params(2)/(1-params(14))+y(3)*T49;
residual(3)= lhs-rhs;
lhs =y(2)*(-(params(17)^(-2)))*(1-params(16)*params(12))-y(6)*params(12)/params(17)-params(16)*params(12)/params(17)*y(11);
rhs =params(1)*y(9);
residual(4)= lhs-rhs;
residual(5) = y(7)+(1-params(4))*y(5)-y(5)-y(2);
lhs =y(11)*params(12)*params(15)+params(12)*y(5);
rhs =y(7);
residual(6)= lhs-rhs;
lhs =y(7);
rhs =params(13)*y(10)+y(5)*params(3)*params(13)/params(15)+y(3)*(1-params(3))*params(13)/params(14);
residual(7)= lhs-rhs;
lhs =y(8);
rhs =y(5)-(1-params(4))*y(5);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =y(2)*T103-y(6)*params(3)*params(13)/(params(17)*params(15))+y(7)*params(3)*(1-params(16))/(params(17)*params(15))-y(5)*T114;
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(10)*params(5)+y(11)*params(7)+params(9)*x(1);
residual(10)= lhs-rhs;
lhs =y(11);
rhs =y(10)*params(8)+y(11)*params(6)+params(10)*x(2);
residual(11)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(11, 11);

  %
  % Jacobian matrix
  %

  g1(1,1)=params(14);
  g1(1,2)=(-1);
  g1(1,3)=params(18);
  g1(1,4)=1;
  g1(2,1)=1;
  g1(2,3)=params(13)*(1-params(3))*(1-params(16))*params(14)^(-2);
  g1(2,6)=(-((-(1-params(3)))*params(13)/params(14)));
  g1(2,7)=(-((1-params(3))*(1-params(16))/params(14)));
  g1(3,1)=1;
  g1(3,2)=(-(params(2)/(1-params(14))));
  g1(3,3)=(-T49);
  g1(4,2)=(-(params(17)^(-2)))*(1-params(16)*params(12));
  g1(4,6)=(-(params(12)/params(17)));
  g1(4,9)=(-params(1));
  g1(4,11)=(-(params(16)*params(12)/params(17)));
  g1(5,2)=(-1);
  g1(5,5)=1-params(4)-1;
  g1(5,7)=1;
  g1(6,5)=params(12);
  g1(6,7)=(-1);
  g1(6,11)=params(12)*params(15);
  g1(7,3)=(-((1-params(3))*params(13)/params(14)));
  g1(7,5)=(-(params(3)*params(13)/params(15)));
  g1(7,7)=1;
  g1(7,10)=(-params(13));
  g1(8,5)=(-(1-(1-params(4))));
  g1(8,8)=1;
  g1(9,2)=(-T103);
  g1(9,5)=T114;
  g1(9,6)=params(3)*params(13)/(params(17)*params(15));
  g1(9,7)=(-(params(3)*(1-params(16))/(params(17)*params(15))));
  g1(9,9)=1;
  g1(10,10)=1-params(5);
  g1(10,11)=(-params(7));
  g1(11,10)=(-params(8));
  g1(11,11)=1-params(6);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],11,121);
end
end

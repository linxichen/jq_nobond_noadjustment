function [residual, g1, g2] = nobond_noadj_static(y, x, params)
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

T28 = (1-y(6))*(1-params(3))*params(11)*exp(y(10))*y(5)^params(3);
T30 = y(3)^(-params(3));
T44 = y(2)^(-1);
T55 = params(3)*exp(y(10))*(1-y(6))*params(11)*y(5)^(params(3)-1);
T56 = y(3)^(1-params(3));
T58 = 1-params(4)+T55*T56;
T102 = getPowerDeriv(y(2),(-1),1);
T128 = T56*params(3)*exp(y(10))*(1-y(6))*params(11)*getPowerDeriv(y(5),params(3)-1,1);
T150 = T56*y(5)^(params(3)-1)*params(3)*exp(y(10))*(-params(11));
residual(1) = y(1)*y(3)+y(4)-y(2);
lhs =T28*T30;
rhs =y(1);
residual(2)= lhs-rhs;
lhs =y(1);
rhs =y(2)*params(2)/(1-y(3));
residual(3)= lhs-rhs;
lhs =(1-y(6)*params(12)*exp(y(11)))*T44;
rhs =T44*params(1)*T58;
residual(4)= lhs-rhs;
residual(5) = y(5)*(1-params(4))+y(7)-y(5)-y(2);
lhs =y(5)*params(12)*exp(y(11));
rhs =y(7);
residual(6)= lhs-rhs;
lhs =y(7);
rhs =T56*y(5)^params(3)*params(11)*exp(y(10));
residual(7)= lhs-rhs;
lhs =y(8);
rhs =y(5)-y(5)*(1-params(4));
residual(8)= lhs-rhs;
lhs =y(9);
rhs =T44*T58;
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(10)*params(5)+y(11)*params(7)+params(9)*x(1);
residual(10)= lhs-rhs;
lhs =y(11);
rhs =y(11)*params(6)+y(10)*params(8)+params(10)*x(2);
residual(11)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(11, 11);

  %
  % Jacobian matrix
  %

  g1(1,1)=y(3);
  g1(1,2)=(-1);
  g1(1,3)=y(1);
  g1(1,4)=1;
  g1(2,1)=(-1);
  g1(2,3)=T28*getPowerDeriv(y(3),(-params(3)),1);
  g1(2,5)=T30*(1-y(6))*(1-params(3))*params(11)*exp(y(10))*getPowerDeriv(y(5),params(3),1);
  g1(2,6)=T30*y(5)^params(3)*exp(y(10))*params(11)*(-(1-params(3)));
  g1(2,10)=T28*T30;
  g1(3,1)=1;
  g1(3,2)=(-(params(2)/(1-y(3))));
  g1(3,3)=(-(y(2)*params(2)/((1-y(3))*(1-y(3)))));
  g1(4,2)=(1-y(6)*params(12)*exp(y(11)))*T102-T58*params(1)*T102;
  g1(4,3)=(-(T44*params(1)*T55*getPowerDeriv(y(3),1-params(3),1)));
  g1(4,5)=(-(T44*params(1)*T128));
  g1(4,6)=T44*(-(params(12)*exp(y(11))))-T44*params(1)*T150;
  g1(4,10)=(-(T44*params(1)*T55*T56));
  g1(4,11)=T44*(-(y(6)*params(12)*exp(y(11))));
  g1(5,2)=(-1);
  g1(5,5)=1-params(4)-1;
  g1(5,7)=1;
  g1(6,5)=params(12)*exp(y(11));
  g1(6,7)=(-1);
  g1(6,11)=y(5)*params(12)*exp(y(11));
  g1(7,3)=(-(y(5)^params(3)*params(11)*exp(y(10))*getPowerDeriv(y(3),1-params(3),1)));
  g1(7,5)=(-(T56*params(11)*exp(y(10))*getPowerDeriv(y(5),params(3),1)));
  g1(7,7)=1;
  g1(7,10)=(-(T56*y(5)^params(3)*params(11)*exp(y(10))));
  g1(8,5)=(-(1-(1-params(4))));
  g1(8,8)=1;
  g1(9,2)=(-(T58*T102));
  g1(9,3)=(-(T44*T55*getPowerDeriv(y(3),1-params(3),1)));
  g1(9,5)=(-(T44*T128));
  g1(9,6)=(-(T44*T150));
  g1(9,9)=1;
  g1(9,10)=(-(T44*T55*T56));
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

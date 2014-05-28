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

residual = zeros( 12, 1);

%
% Model equations
%

T24 = y(5)^params(4);
T27 = y(3)^(-params(4));
T41 = y(2)^(-params(9));
T50 = y(5)^(params(4)-1);
T52 = y(3)^(1-params(4));
T54 = 1-params(5)+params(4)*(1-y(6))*y(11)*T50*T52;
T75 = y(11)/params(16);
T79 = y(12)/params(17);
T81 = T79^params(12);
T108 = getPowerDeriv(y(2),(-params(9)),1);
residual(1) = y(1)*y(3)+y(4)-y(2);
lhs =(1-y(6))*(1-params(4))*y(11)*T24*T27;
rhs =y(1);
residual(2)= lhs-rhs;
lhs =y(1);
rhs =params(3)*y(2)^params(9)/(1-y(3));
residual(3)= lhs-rhs;
lhs =(1-y(6)*y(12))*T41;
rhs =T41*params(1)*T54;
residual(4)= lhs-rhs;
residual(5) = y(5)*(1-params(5))+y(7)-y(5)-y(2);
lhs =y(5)*y(12);
rhs =y(7);
residual(6)= lhs-rhs;
lhs =y(7);
rhs =T52*y(11)*T24;
residual(7)= lhs-rhs;
lhs =y(8);
rhs =y(5)-y(5)*(1-params(5));
residual(8)= lhs-rhs;
lhs =y(9);
rhs =T41*T54;
residual(9)= lhs-rhs;
lhs =T75;
rhs =T75^params(10)*T81*exp(params(14)*x(1));
residual(10)= lhs-rhs;
lhs =T79;
rhs =T79^params(11)*T75^params(13)*exp((-params(15))*x(2));
residual(11)= lhs-rhs;
lhs =y(10);
rhs =y(7)/y(3);
residual(12)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(12, 12);

  %
  % Jacobian matrix
  %

  g1(1,1)=y(3);
  g1(1,2)=(-1);
  g1(1,3)=y(1);
  g1(1,4)=1;
  g1(2,1)=(-1);
  g1(2,3)=(1-y(6))*(1-params(4))*y(11)*T24*getPowerDeriv(y(3),(-params(4)),1);
  g1(2,5)=T27*(1-y(6))*(1-params(4))*y(11)*getPowerDeriv(y(5),params(4),1);
  g1(2,6)=T27*T24*y(11)*(-(1-params(4)));
  g1(2,11)=T27*(1-y(6))*(1-params(4))*T24;
  g1(3,1)=1;
  g1(3,2)=(-(params(3)*getPowerDeriv(y(2),params(9),1)/(1-y(3))));
  g1(3,3)=(-(params(3)*y(2)^params(9)/((1-y(3))*(1-y(3)))));
  g1(4,2)=(1-y(6)*y(12))*T108-T54*params(1)*T108;
  g1(4,3)=(-(T41*params(1)*params(4)*(1-y(6))*y(11)*T50*getPowerDeriv(y(3),1-params(4),1)));
  g1(4,5)=(-(T41*params(1)*T52*params(4)*(1-y(6))*y(11)*getPowerDeriv(y(5),params(4)-1,1)));
  g1(4,6)=T41*(-y(12))-T41*params(1)*T52*T50*params(4)*(-y(11));
  g1(4,11)=(-(T41*params(1)*T52*T50*(1-y(6))*params(4)));
  g1(4,12)=T41*(-y(6));
  g1(5,2)=(-1);
  g1(5,5)=1-params(5)-1;
  g1(5,7)=1;
  g1(6,5)=y(12);
  g1(6,7)=(-1);
  g1(6,12)=y(5);
  g1(7,3)=(-(y(11)*T24*getPowerDeriv(y(3),1-params(4),1)));
  g1(7,5)=(-(T52*y(11)*getPowerDeriv(y(5),params(4),1)));
  g1(7,7)=1;
  g1(7,11)=(-(T24*T52));
  g1(8,5)=(-(1-(1-params(5))));
  g1(8,8)=1;
  g1(9,2)=(-(T54*T108));
  g1(9,3)=(-(T41*params(4)*(1-y(6))*y(11)*T50*getPowerDeriv(y(3),1-params(4),1)));
  g1(9,5)=(-(T41*T52*params(4)*(1-y(6))*y(11)*getPowerDeriv(y(5),params(4)-1,1)));
  g1(9,6)=(-(T41*T52*T50*params(4)*(-y(11))));
  g1(9,9)=1;
  g1(9,11)=(-(T41*T52*T50*(1-y(6))*params(4)));
  g1(10,11)=1/params(16)-exp(params(14)*x(1))*T81*1/params(16)*getPowerDeriv(T75,params(10),1);
  g1(10,12)=(-(exp(params(14)*x(1))*T75^params(10)*1/params(17)*getPowerDeriv(T79,params(12),1)));
  g1(11,11)=(-(exp((-params(15))*x(2))*T79^params(11)*1/params(16)*getPowerDeriv(T75,params(13),1)));
  g1(11,12)=1/params(17)-exp((-params(15))*x(2))*T75^params(13)*1/params(17)*getPowerDeriv(T79,params(11),1);
  g1(12,3)=(-((-y(7))/(y(3)*y(3))));
  g1(12,7)=(-(1/y(3)));
  g1(12,10)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],12,144);
end
end

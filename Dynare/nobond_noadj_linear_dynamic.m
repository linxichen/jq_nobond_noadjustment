function [residual, g1, g2, g3] = nobond_noadj_linear_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           columns: equations in order of declaration
%                                                           rows: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(11, 1);
T49 = params(2)*params(17)/(1-params(14))^2;
T105 = (-(params(17)^(-2)))*(1-params(4)+params(13)*params(3)*(1-params(16))/params(15));
T116 = params(13)*params(3)*(1-params(16))/(params(17)*params(15)^2);
lhs =params(14)*y(4)+params(18)*y(6)+y(7);
rhs =y(5);
residual(1)= lhs-rhs;
lhs =y(4);
rhs =(-(1-params(3)))*params(13)/params(14)*y(9)+(1-params(3))*(1-params(16))/params(14)*y(10)-y(6)*params(13)*(1-params(3))*(1-params(16))*params(14)^(-2);
residual(2)= lhs-rhs;
lhs =y(4);
rhs =y(5)*params(2)/(1-params(14))+y(6)*T49;
residual(3)= lhs-rhs;
lhs =y(5)*(-(params(17)^(-2)))*(1-params(16)*params(12))-y(9)*params(12)/params(17)-params(16)*params(12)/params(17)*y(14);
rhs =params(1)*y(15);
residual(4)= lhs-rhs;
residual(5) = y(10)+(1-params(4))*y(1)-y(8)-y(5);
lhs =y(14)*params(12)*params(15)+params(12)*y(8);
rhs =y(10);
residual(6)= lhs-rhs;
lhs =y(10);
rhs =params(13)*y(13)+y(1)*params(3)*params(13)/params(15)+y(6)*(1-params(3))*params(13)/params(14);
residual(7)= lhs-rhs;
lhs =y(11);
rhs =y(8)-(1-params(4))*y(1);
residual(8)= lhs-rhs;
lhs =y(12);
rhs =y(5)*T105-y(9)*params(3)*params(13)/(params(17)*params(15))+y(10)*params(3)*(1-params(16))/(params(17)*params(15))-y(1)*T116;
residual(9)= lhs-rhs;
lhs =y(13);
rhs =params(5)*y(2)+params(7)*y(3)+params(9)*x(it_, 1);
residual(10)= lhs-rhs;
lhs =y(14);
rhs =y(2)*params(8)+y(3)*params(6)+params(10)*x(it_, 2);
residual(11)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(11, 17);

  %
  % Jacobian matrix
  %

  g1(1,4)=params(14);
  g1(1,5)=(-1);
  g1(1,6)=params(18);
  g1(1,7)=1;
  g1(2,4)=1;
  g1(2,6)=params(13)*(1-params(3))*(1-params(16))*params(14)^(-2);
  g1(2,9)=(-((-(1-params(3)))*params(13)/params(14)));
  g1(2,10)=(-((1-params(3))*(1-params(16))/params(14)));
  g1(3,4)=1;
  g1(3,5)=(-(params(2)/(1-params(14))));
  g1(3,6)=(-T49);
  g1(4,5)=(-(params(17)^(-2)))*(1-params(16)*params(12));
  g1(4,9)=(-(params(12)/params(17)));
  g1(4,15)=(-params(1));
  g1(4,14)=(-(params(16)*params(12)/params(17)));
  g1(5,5)=(-1);
  g1(5,1)=1-params(4);
  g1(5,8)=(-1);
  g1(5,10)=1;
  g1(6,8)=params(12);
  g1(6,10)=(-1);
  g1(6,14)=params(12)*params(15);
  g1(7,6)=(-((1-params(3))*params(13)/params(14)));
  g1(7,1)=(-(params(3)*params(13)/params(15)));
  g1(7,10)=1;
  g1(7,13)=(-params(13));
  g1(8,1)=1-params(4);
  g1(8,8)=(-1);
  g1(8,11)=1;
  g1(9,5)=(-T105);
  g1(9,1)=T116;
  g1(9,9)=params(3)*params(13)/(params(17)*params(15));
  g1(9,10)=(-(params(3)*(1-params(16))/(params(17)*params(15))));
  g1(9,12)=1;
  g1(10,2)=(-params(5));
  g1(10,13)=1;
  g1(10,3)=(-params(7));
  g1(10,16)=(-params(9));
  g1(11,2)=(-params(8));
  g1(11,3)=(-params(6));
  g1(11,14)=1;
  g1(11,17)=(-params(10));
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],11,289);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],11,4913);
end
end

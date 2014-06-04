function [residual, g1, g2, g3] = nobond_noadj_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(12, 1);
T24 = y(1)^params(3);
T27 = y(6)^(-params(3));
T38 = y(5)^(-1);
T43 = params(1)*y(16)^(-1);
T53 = y(8)^(params(3)-1);
T56 = y(17)^(1-params(3));
T58 = 1-params(4)+params(3)*(1-y(18))*y(19)*T53*T56;
T70 = y(6)^(1-params(3));
T79 = y(1)^(params(3)-1);
T82 = 1-params(4)+T70*params(3)*(1-y(9))*y(14)*T79;
T88 = y(2)/params(11);
T90 = T88^params(5);
T93 = y(3)/params(12);
T95 = T93^params(7);
T105 = T93^params(6);
residual(1) = y(4)*y(6)+y(7)-y(5);
lhs =(1-y(9))*(1-params(3))*y(14)*T24*T27;
rhs =y(4);
residual(2)= lhs-rhs;
lhs =y(4);
rhs =y(5)*params(2)/(1-y(6));
residual(3)= lhs-rhs;
lhs =(1-y(9)*y(15))*T38;
rhs =T43*T58;
residual(4)= lhs-rhs;
residual(5) = y(1)*(1-params(4))+y(10)-y(8)-y(5);
lhs =y(15)*y(8);
rhs =y(10);
residual(6)= lhs-rhs;
lhs =y(10);
rhs =y(14)*T24*T70;
residual(7)= lhs-rhs;
lhs =y(11);
rhs =y(8)-y(1)*(1-params(4));
residual(8)= lhs-rhs;
lhs =y(12);
rhs =T38*T82;
residual(9)= lhs-rhs;
lhs =y(14)/params(11);
rhs =T90*T95*exp(params(9)*x(it_, 1));
residual(10)= lhs-rhs;
lhs =y(15)/params(12);
rhs =T105*T88^params(8)*exp(params(10)*x(it_, 2));
residual(11)= lhs-rhs;
lhs =y(13);
rhs =y(10)/y(6);
residual(12)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(12, 21);

  %
  % Jacobian matrix
  %

  g1(1,4)=y(6);
  g1(1,5)=(-1);
  g1(1,6)=y(4);
  g1(1,7)=1;
  g1(2,4)=(-1);
  g1(2,6)=(1-y(9))*(1-params(3))*y(14)*T24*getPowerDeriv(y(6),(-params(3)),1);
  g1(2,1)=T27*(1-y(9))*(1-params(3))*y(14)*getPowerDeriv(y(1),params(3),1);
  g1(2,9)=T27*T24*y(14)*(-(1-params(3)));
  g1(2,14)=T27*(1-y(9))*(1-params(3))*T24;
  g1(3,4)=1;
  g1(3,5)=(-(params(2)/(1-y(6))));
  g1(3,6)=(-(y(5)*params(2)/((1-y(6))*(1-y(6)))));
  g1(4,5)=(1-y(9)*y(15))*getPowerDeriv(y(5),(-1),1);
  g1(4,16)=(-(T58*params(1)*getPowerDeriv(y(16),(-1),1)));
  g1(4,17)=(-(T43*params(3)*(1-y(18))*y(19)*T53*getPowerDeriv(y(17),1-params(3),1)));
  g1(4,8)=(-(T43*T56*params(3)*(1-y(18))*y(19)*getPowerDeriv(y(8),params(3)-1,1)));
  g1(4,9)=T38*(-y(15));
  g1(4,18)=(-(T43*T56*T53*params(3)*(-y(19))));
  g1(4,19)=(-(T43*T56*T53*params(3)*(1-y(18))));
  g1(4,15)=T38*(-y(9));
  g1(5,5)=(-1);
  g1(5,1)=1-params(4);
  g1(5,8)=(-1);
  g1(5,10)=1;
  g1(6,8)=y(15);
  g1(6,10)=(-1);
  g1(6,15)=y(8);
  g1(7,6)=(-(y(14)*T24*getPowerDeriv(y(6),1-params(3),1)));
  g1(7,1)=(-(T70*y(14)*getPowerDeriv(y(1),params(3),1)));
  g1(7,10)=1;
  g1(7,14)=(-(T24*T70));
  g1(8,1)=1-params(4);
  g1(8,8)=(-1);
  g1(8,11)=1;
  g1(9,5)=(-(T82*getPowerDeriv(y(5),(-1),1)));
  g1(9,6)=(-(T38*params(3)*(1-y(9))*y(14)*T79*getPowerDeriv(y(6),1-params(3),1)));
  g1(9,1)=(-(T38*T70*params(3)*(1-y(9))*y(14)*getPowerDeriv(y(1),params(3)-1,1)));
  g1(9,9)=(-(T38*T70*T79*params(3)*(-y(14))));
  g1(9,12)=1;
  g1(9,14)=(-(T38*T70*T79*(1-y(9))*params(3)));
  g1(10,2)=(-(exp(params(9)*x(it_, 1))*T95*1/params(11)*getPowerDeriv(T88,params(5),1)));
  g1(10,14)=1/params(11);
  g1(10,3)=(-(exp(params(9)*x(it_, 1))*T90*1/params(12)*getPowerDeriv(T93,params(7),1)));
  g1(10,20)=(-(T90*T95*params(9)*exp(params(9)*x(it_, 1))));
  g1(11,2)=(-(exp(params(10)*x(it_, 2))*T105*1/params(11)*getPowerDeriv(T88,params(8),1)));
  g1(11,3)=(-(exp(params(10)*x(it_, 2))*T88^params(8)*1/params(12)*getPowerDeriv(T93,params(6),1)));
  g1(11,15)=1/params(12);
  g1(11,21)=(-(T105*T88^params(8)*params(10)*exp(params(10)*x(it_, 2))));
  g1(12,6)=(-((-y(10))/(y(6)*y(6))));
  g1(12,10)=(-(1/y(6)));
  g1(12,13)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],12,441);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],12,9261);
end
end

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
T24 = y(1)^params(4);
T27 = y(6)^(-params(4));
T41 = y(5)^(-params(9));
T46 = params(1)*y(16)^(-params(9));
T56 = y(8)^(params(4)-1);
T59 = y(17)^(1-params(4));
T61 = 1-params(5)+params(4)*(1-y(18))*y(19)*T56*T59;
T73 = y(6)^(1-params(4));
T82 = y(1)^(params(4)-1);
T85 = 1-params(5)+T73*params(4)*(1-y(9))*y(14)*T82;
T91 = y(2)/params(16);
T93 = T91^params(10);
T96 = y(3)/params(17);
T98 = T96^params(12);
T108 = T96^params(11);
residual(1) = y(4)*y(6)+y(7)-y(5);
lhs =(1-y(9))*(1-params(4))*y(14)*T24*T27;
rhs =y(4);
residual(2)= lhs-rhs;
lhs =y(4);
rhs =params(3)*y(5)^params(9)/(1-y(6));
residual(3)= lhs-rhs;
lhs =(1-y(9)*y(15))*T41;
rhs =T46*T61;
residual(4)= lhs-rhs;
residual(5) = y(1)*(1-params(5))+y(10)-y(8)-y(5);
lhs =y(15)*y(8);
rhs =y(10);
residual(6)= lhs-rhs;
lhs =y(10);
rhs =y(14)*T24*T73;
residual(7)= lhs-rhs;
lhs =y(11);
rhs =y(8)-y(1)*(1-params(5));
residual(8)= lhs-rhs;
lhs =y(12);
rhs =T41*T85;
residual(9)= lhs-rhs;
lhs =y(14)/params(16);
rhs =T93*T98*exp(params(14)*x(it_, 1));
residual(10)= lhs-rhs;
lhs =y(15)/params(17);
rhs =T108*T91^params(13)*exp((-params(15))*x(it_, 2));
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
  g1(2,6)=(1-y(9))*(1-params(4))*y(14)*T24*getPowerDeriv(y(6),(-params(4)),1);
  g1(2,1)=T27*(1-y(9))*(1-params(4))*y(14)*getPowerDeriv(y(1),params(4),1);
  g1(2,9)=T27*T24*y(14)*(-(1-params(4)));
  g1(2,14)=T27*(1-y(9))*(1-params(4))*T24;
  g1(3,4)=1;
  g1(3,5)=(-(params(3)*getPowerDeriv(y(5),params(9),1)/(1-y(6))));
  g1(3,6)=(-(params(3)*y(5)^params(9)/((1-y(6))*(1-y(6)))));
  g1(4,5)=(1-y(9)*y(15))*getPowerDeriv(y(5),(-params(9)),1);
  g1(4,16)=(-(T61*params(1)*getPowerDeriv(y(16),(-params(9)),1)));
  g1(4,17)=(-(T46*params(4)*(1-y(18))*y(19)*T56*getPowerDeriv(y(17),1-params(4),1)));
  g1(4,8)=(-(T46*T59*params(4)*(1-y(18))*y(19)*getPowerDeriv(y(8),params(4)-1,1)));
  g1(4,9)=T41*(-y(15));
  g1(4,18)=(-(T46*T59*T56*params(4)*(-y(19))));
  g1(4,19)=(-(T46*T59*T56*params(4)*(1-y(18))));
  g1(4,15)=T41*(-y(9));
  g1(5,5)=(-1);
  g1(5,1)=1-params(5);
  g1(5,8)=(-1);
  g1(5,10)=1;
  g1(6,8)=y(15);
  g1(6,10)=(-1);
  g1(6,15)=y(8);
  g1(7,6)=(-(y(14)*T24*getPowerDeriv(y(6),1-params(4),1)));
  g1(7,1)=(-(T73*y(14)*getPowerDeriv(y(1),params(4),1)));
  g1(7,10)=1;
  g1(7,14)=(-(T24*T73));
  g1(8,1)=1-params(5);
  g1(8,8)=(-1);
  g1(8,11)=1;
  g1(9,5)=(-(T85*getPowerDeriv(y(5),(-params(9)),1)));
  g1(9,6)=(-(T41*params(4)*(1-y(9))*y(14)*T82*getPowerDeriv(y(6),1-params(4),1)));
  g1(9,1)=(-(T41*T73*params(4)*(1-y(9))*y(14)*getPowerDeriv(y(1),params(4)-1,1)));
  g1(9,9)=(-(T41*T73*T82*params(4)*(-y(14))));
  g1(9,12)=1;
  g1(9,14)=(-(T41*T73*T82*(1-y(9))*params(4)));
  g1(10,2)=(-(exp(params(14)*x(it_, 1))*T98*1/params(16)*getPowerDeriv(T91,params(10),1)));
  g1(10,14)=1/params(16);
  g1(10,3)=(-(exp(params(14)*x(it_, 1))*T93*1/params(17)*getPowerDeriv(T96,params(12),1)));
  g1(10,20)=(-(T93*T98*params(14)*exp(params(14)*x(it_, 1))));
  g1(11,2)=(-(exp((-params(15))*x(it_, 2))*T108*1/params(16)*getPowerDeriv(T91,params(13),1)));
  g1(11,3)=(-(exp((-params(15))*x(it_, 2))*T91^params(13)*1/params(17)*getPowerDeriv(T96,params(11),1)));
  g1(11,15)=1/params(17);
  g1(11,21)=(-(T108*T91^params(13)*(-params(15))*exp((-params(15))*x(it_, 2))));
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

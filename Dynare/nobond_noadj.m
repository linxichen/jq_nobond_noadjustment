%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'nobond_noadj';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('nobond_noadj.log');
M_.exo_names = 'epsz';
M_.exo_names_tex = 'epsz';
M_.exo_names_long = 'epsz';
M_.exo_names = char(M_.exo_names, 'epsxxi');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsxxi');
M_.exo_names_long = char(M_.exo_names_long, 'epsxxi');
M_.endo_names = 'w';
M_.endo_names_tex = 'w';
M_.endo_names_long = 'w';
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names_long = char(M_.endo_names_long, 'n');
M_.endo_names = char(M_.endo_names, 'd');
M_.endo_names_tex = char(M_.endo_names_tex, 'd');
M_.endo_names_long = char(M_.endo_names_long, 'd');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'mmu');
M_.endo_names_tex = char(M_.endo_names_tex, 'mmu');
M_.endo_names_long = char(M_.endo_names_long, 'mmu');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'inv');
M_.endo_names_tex = char(M_.endo_names_tex, 'inv');
M_.endo_names_long = char(M_.endo_names_long, 'inv');
M_.endo_names = char(M_.endo_names, 'mk');
M_.endo_names_tex = char(M_.endo_names_tex, 'mk');
M_.endo_names_long = char(M_.endo_names_long, 'mk');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.endo_names = char(M_.endo_names, 'xxi');
M_.endo_names_tex = char(M_.endo_names_tex, 'xxi');
M_.endo_names_long = char(M_.endo_names_long, 'xxi');
M_.param_names = 'bbeta';
M_.param_names_tex = 'bbeta';
M_.param_names_long = 'bbeta';
M_.param_names = char(M_.param_names, 'aalpha');
M_.param_names_tex = char(M_.param_names_tex, 'aalpha');
M_.param_names_long = char(M_.param_names_long, 'aalpha');
M_.param_names = char(M_.param_names, 'ttheta');
M_.param_names_tex = char(M_.param_names_tex, 'ttheta');
M_.param_names_long = char(M_.param_names_long, 'ttheta');
M_.param_names = char(M_.param_names, 'ddelta');
M_.param_names_tex = char(M_.param_names_tex, 'ddelta');
M_.param_names_long = char(M_.param_names_long, 'ddelta');
M_.param_names = char(M_.param_names, 'rrhozz');
M_.param_names_tex = char(M_.param_names_tex, 'rrhozz');
M_.param_names_long = char(M_.param_names_long, 'rrhozz');
M_.param_names = char(M_.param_names, 'rrhoxxixxi');
M_.param_names_tex = char(M_.param_names_tex, 'rrhoxxixxi');
M_.param_names_long = char(M_.param_names_long, 'rrhoxxixxi');
M_.param_names = char(M_.param_names, 'rrhozxxi');
M_.param_names_tex = char(M_.param_names_tex, 'rrhozxxi');
M_.param_names_long = char(M_.param_names_long, 'rrhozxxi');
M_.param_names = char(M_.param_names, 'rrhoxxiz');
M_.param_names_tex = char(M_.param_names_tex, 'rrhoxxiz');
M_.param_names_long = char(M_.param_names_long, 'rrhoxxiz');
M_.param_names = char(M_.param_names, 'ssigmaepsz');
M_.param_names_tex = char(M_.param_names_tex, 'ssigmaepsz');
M_.param_names_long = char(M_.param_names_long, 'ssigmaepsz');
M_.param_names = char(M_.param_names, 'ssigmaepsxxi');
M_.param_names_tex = char(M_.param_names_tex, 'ssigmaepsxxi');
M_.param_names_long = char(M_.param_names_long, 'ssigmaepsxxi');
M_.param_names = char(M_.param_names, 'zbar');
M_.param_names_tex = char(M_.param_names_tex, 'zbar');
M_.param_names_long = char(M_.param_names_long, 'zbar');
M_.param_names = char(M_.param_names, 'xxibar');
M_.param_names_tex = char(M_.param_names_tex, 'xxibar');
M_.param_names_long = char(M_.param_names_long, 'xxibar');
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 11;
M_.param_nbr = 12;
M_.orig_endo_nbr = 11;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('nobond_noadj_static');
erase_compiled_function('nobond_noadj_dynamic');
M_.lead_lag_incidence = [
 0 4 0;
 0 5 15;
 0 6 16;
 0 7 0;
 1 8 0;
 0 9 17;
 0 10 0;
 0 11 0;
 0 12 0;
 2 13 18;
 3 14 0;]';
M_.nstatic = 5;
M_.nfwrd   = 3;
M_.npred   = 2;
M_.nboth   = 1;
M_.nsfwrd   = 4;
M_.nspred   = 3;
M_.ndynamic   = 6;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(11, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(12, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 48;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
cd ../MATLAB;
mypara;
cd ../Dynare
for i=1:length(M_.params)
deep_parameter_name = M_.param_names(i,:);
eval(['M_.params(i)  = ' deep_parameter_name ' ;'])
end
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
M_.Sigma_e(2, 2) = (1)^2;
M_.sigma_e_is_diagonal = 1;
options_.drop = 2000;
options_.hp_filter = 1600;
options_.irf = 40;
options_.order = 1;
options_.periods = 100000;
var_list_=[];
var_list_ = 'y';
var_list_ = char(var_list_, 'c');
var_list_ = char(var_list_, 'inv');
var_list_ = char(var_list_, 'n');
var_list_ = char(var_list_, 'w');
var_list_ = char(var_list_, 'mmu');
var_list_ = char(var_list_, 'z');
var_list_ = char(var_list_, 'd');
info = stoch_simul(var_list_);
save('nobond_noadj_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('nobond_noadj_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('nobond_noadj_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('nobond_noadj_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('nobond_noadj_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off

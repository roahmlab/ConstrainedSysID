%% SysID Main Script

% Author: Daniel Haugk, 2024, University of Michigan

%% set algorithm options

% error filter value for weights
k = 0.3;

% stopping value for while loops
tol = 0.0001;

% number of starting points for first multi search optimization
MS1 = 3;

% number of starting points for second (friction) multi search optimization
MS2 = 3;

% include regrouping or not
regroup = false;

% search for all parameters (inertial, friction, motor, offset) or not
SearchAll = true;

% search for friction motor and offset parameters or not
SearchFM = true;

% optimize over friction exponents or not
SearchAlpha = false;

% include offset parameter or not
includeOffset = false;

% include friction and motor dynamics or not (only needed when SearchFM is false)
includeFMDynamics = true;

% system includes constraints or not
includeConstraints = true;

% include normalization and weighting
includeWeighting = true;

% choose between projection matrix for general constraints (1) or fully actuated systems (2)
% projection matrix for fully actuated systems might be less computationally expensive!
constraintVariant = 1;

% fmincon options
options = optimoptions('fmincon','MaxFunctionEvaluations',30000,'Display','iter','Algorithm','sqp','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'UseParallel',true);

% collect algorithm options in a struct
AlgOptions = {k,tol,MS1,MS2,regroup,SearchAll,SearchFM,SearchAlpha,includeOffset,includeFMDynamics,includeConstraints,constraintVariant,includeWeighting};


%% set initial conditions or true values

% initial condition for inertial parameters (if they are known use the true values)
% the structure is : XX, XY, XZ, YY, YZ, ZZ, hx, hy, hz, m and repeats itself for each link
X0_ip = [];

% true inertial parameters, if known
X_ip = [];

% initial condition for motor,friction and offset parameters (if they are known use the true values)
% the structure is : Im_1,...,Im_n , Fc_1, Fv_1, beta_1,...., Fc_n, Fv_n, beta_n
X0_fm = [];

% true motor, friction and offset parameters if known
X_fm = [];

% full initial condition vector
X0_1 = [X0_ip;X0_fm];

%%  set lower and upper bounds

% lower bound for inertial parameters (same structure as initial condition)
lb_ip = []';

% lower bound for motor, friction and offset parameters (same structure as initial condition)
lb_fm = []';

% upper bound for inertial parameters (same structure as initial condition)
ub_ip = []';

% upper bound for motor, friction and offset parameters (same structure as initial condition)
ub_fm = []';

% lower bound for friction exponent
% structure: alpha_1,...,alpha_n
lba = []';

% upper bound for friction exponent
% structure: alpha_1,...,alpha_n
uba = []';

% collect bounds in one vector
lb = [lb_ip;lb_fm];
ub = [ub_ip;ub_fm];

%% get data

% this script can be run by itself after setting the algorithm options
run('SysIDData.m');

%% run the system identification

[X, Wfull, alphanew, Ginv] = SysIDAlg(AlgOptions,na_idx,nu_idx,J,lb,ub,lba,uba,T,data,dataFull,W_ip,X0_1,X_ip,X_fm,options);

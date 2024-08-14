%% SysID Algorithm Script

% Script that runs an algorithm specified through the algorithm options

% Author: Daniel Haugk, 2024, University of Michigan

function [X, Wfull, alphanew, Ginv] = SysIDAlg(AlgOptions,na_idx,nu_idx,J,lb,ub,lba,uba,T,data,dataFull,W_ip,X0_1,X_ip,X_fm,options)

% initialize variables
n = length(data(:,1));
m = length(data(1,:))/3;

% initalize algorithm options
regroup = AlgOptions{5};
SearchAll = AlgOptions{6};
SearchFM = AlgOptions{7};
includeWeighting = AlgOptions{13};

% check if regrouping is wanted
if regroup

    % get neceserray things for the regrouping trough QR decomposition
    [Aid, Ad, Kd, ~, Ginv, b, d] = QRDecomposition(W_ip,ones(10*n,1));

else

    % if no regrouping is wanted set variables to empty
    Ginv = [];
    Aid = [];
    Ad = [];
    Kd = [];
    b = [];
    d = [];

end

% assign variables in struct
variables = {n,m,b,d};

%% optimization

% choose between different cases, depending on the wanted algorithm options

% check if all parameters (inertial, motor, friction, offset) need to be found
if SearchAll
    
    % check if normalization and weighting should be included
    if includeWeighting

        [X,Wfull,alphanew] = SysIDAlgFull(variables,AlgOptions,na_idx,nu_idx,J,lb,ub,lba,uba,T,data,dataFull,W_ip,Ginv,Aid,Ad,Kd,X0_1,options);

    else

        [X,Wfull,alphanew] = SysIDAlgFullNoWeights(variables,AlgOptions,na_idx,nu_idx,J,lb,ub,lba,uba,T,data,dataFull,W_ip,Ginv,Aid,Ad,Kd,X0_1,options);

    end

% check if motor, friction and offset parameters need to be found
elseif SearchFM 

    % check if normalization and weighting should be included
    if includeWeighting

        [X,Wfull,alphanew] = SysIDAlgFMParamOnly(variables,AlgOptions,na_idx,nu_idx,J,lb,ub,lba,uba,X_ip,T,data,dataFull,W_ip,X0_1,options);

    else

        [X,Wfull,alphanew] = SysIDAlgFMParamOnlyNoWeights(variables,AlgOptions,na_idx,nu_idx,J,lb,ub,lba,uba,X_ip,T,data,dataFull,W_ip,X0_1,options);

    end

% % check if inertial parameters need to be found
else

    % check if normalization and weighting should be included
    if includeWeighting

        [X,Wfull,alphanew] = SysIDAlgInParamOnly(variables,AlgOptions,na_idx,nu_idx,J,lb,ub,X_fm,T,data,dataFull,W_ip,Ginv,Aid,Ad,Kd,X0_1,options);

    else

        [X,Wfull,alphanew] = SysIDAlgInParamOnlyNoWeights(variables,AlgOptions,na_idx,nu_idx,J,lb,ub,X_fm,T,data,dataFull,W_ip,Ginv,Aid,Ad,Kd,X0_1,options);

    end

end


end
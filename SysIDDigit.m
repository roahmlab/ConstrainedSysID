%% Digit specific data and functions for sysID

% This script is for Digits Limbs only.
% The script is split into 2 parts. 
% In part 1 or step 1 the raw data is extracted from some data file.
% In part 2 or step 2 the processed raw data is used to compute the observation matrix, torque vector and to construct the data matrices

% A similar script has to be created for each new robot!

% A similar script must have the outputs:

% t_ = raw time vector
% pos = raw position vector
% vel = raw velocity vector
% J = constraint Jacobian function
% na_idx = indices of actuated joints of corresponding limb
% nu_idx = indices of unactuated joints of corresponding limb
% W_ip = observation matrix for inertial parameters
% T = torque vector
% data = data matrix containing processed position, velocity and acceleration data of the limb that has to be identified
% dataFull = data matrix containing the data of all joints (if the model was specified as a full model, instead of only a model of the corresponding limb)

% Author: Daniel Haugk, 2024, University of Michigan

function [t_,pos,vel,torq_,J,na_idx,nu_idx,W_ip,T,data,dataFull] = SysIDDigit(datapaths,bodypart,step,q,qdot,qddot,torq,i,constraintVariant)

    % check which bodypart of digit is being identified
    if bodypart == "LeftLeg"
        
        % create model of digit
        model = createFloatingBaseFixedSpringDigit;

        % data rows that correspond to the left leg
        s = 1;
        f = 6;
        
       % create function for constraint jacobian 
       J = @(qv) model.getJLeft(model,qv);

       % specify the index of actuated/unactuated joints
       na_idx = [1 2 3 4 15 16];
       nu_idx = [5 6 7 19 20 21 22 27 28];

       % specify the index of the inertial parameters of the left leg
       p_idx = [1:10 11:20 21:30 31:40 41:50 51:60 61:70 141:150 151:160 181:190 191:200 201:210 211:220 261:270 271:280];

    elseif bodypart == "RightLeg"

        % create model of digit
        model = createFloatingBaseFixedSpringDigit;

        % data rows that correspond to the right leg
        s = 7;
        f = 12;

        % create function for constraint jacobian 
        J = @(qv) model.getJRight(model,qv);

        % specify the index of actuated/unactuated joints
        na_idx = [8 9 10 11 17 18];
        nu_idx = [12 13 14 23 24 25 26 29 30];

        % specify the index of the inertial parameters of the right leg
        p_idx = [71:80 81:90 91:100 101:110 111:120 121:130 131:140 161:170 171:180 221:230 231:240 241:250 251:260 281:290 291:300];

    elseif bodypart == "LeftArm"

        % create model of digits left arm via urdf
        model = create_model_from_urdf('digit-v3-leftarmonly.urdf');

        % data rows that correspond to the left arm
        s = 13;
        f = 16;

        % no constraints for the arms
        J = [];

        % all joints are actuated
        na_idx = [];
        nu_idx = [];

    elseif bodypart == "RightArm"

        % create model of digits right arm via urdf
        model = create_model_from_urdf('digit-v3-rightarmonly.urdf');

        % data rows that correspond to the right arm
        s = 17;
        f = 20;

        % no constraints for the arms
        J = [];

        % all joints are actuated
        na_idx = [];
        nu_idx = [];

    else 
    
        disp("choose a bodypart")

    end

    % first step before data processing, to get the raw data
    if step == 1

        % set variables not used in the first step to empty
        W_ip = [];
        T = [];
        dataFull = [];
    
        % read the raw data from the specified datapaths
        data = readmatrix(datapaths);
    
        % extract position (q), velocity (qd), observed torque (tau_obs) and time (t) from the data matrix
        q = data(:,1:20)';
        qd = data(:,21:40)';
        tau_obs = data(:,41:60)';
        t = data(:,end)';
    
        % use only unique values of t
        [t_, ia, ~] = unique(t);
        t_ = t_(~isnan(t_));
    
        % only use values corresponding to the unique values of t
        m = 1;
        for i = s:f
    
            q_ = q(i,:);
            qd_ = qd(i,:);
            torq__ = tau_obs(i,:);
    
            pos(m,:) = q_(ia);
            vel(m,:) = qd_(ia);
            torq_(m,:) = torq__(ia);
    
            m = m+1;
        end
    
        % remove the last value
        pos(:,end) = [];
        vel(:,end) = [];
        torq_(:,end) = [];
    
    % second step after data processing, to get observation matrix, torque vector and specific data matrices
    elseif step == 2
        
        % set variables not used in this step to empty
        t_ = [];
        torq_ = [];
        
        % check if a leg or arm is being identified
        if bodypart == "LeftLeg" || "RightLeg"
            
            % create zero vector of the full state vectors
            q0 = zeros(30,length(q));
            qd0 = zeros(30,length(q));
            qdd0 = zeros(30,length(q));

            % assign known states from actuated joints to the full state vectors trough their index
            q0(na_idx,:) = q;
            qd0(na_idx,:) = qdot;
            qdd0(na_idx,:) = qddot;
            
            % compute all unknown (unactuated) states, s.t. they satisfy the constraints
            [pos, vel, acc] = model.fillUnactJoints(model, q0, qd0, qdd0);

            % create zero torque matrix
            torq0 = zeros(30,length(q));

            % assign known torques to the full torque matrix through their index
            torq0(na_idx,:) = torq;

            % set torque matrix
            torq = torq0;

            % initialize observation Matrix
            W_ip = [];

            % compute observation matrix for the left/right leg, depending on the chosen indices
            for j = 1:i
                [Y_ip, ~, ~] = Yphi(model, pos(:,j),vel(:,j), acc(:,j));
                W_ip = [W_ip; Y_ip(sort([na_idx nu_idx]),p_idx)];
            end

            % make the observation matrix a double matrix
            W_ip = double(W_ip);

            % get the torque vector depending on the chosen constraint variant
            if constraintVariant == 1
                % torque vector for all joints of the left/right leg
                T = double(reshape(torq(sort([na_idx nu_idx]),:),[length(sort([na_idx nu_idx]))*length(q) 1]));
            elseif constraintVariant == 2
                % torque vector for actuated joints only of the left/right leg
                T = double(reshape(torq(na_idx,:),[length(na_idx)*length(q) 1]));
            end

            % data measurement vectors
            dataFull = double([pos vel acc]);
            data = dataFull(sort([na_idx nu_idx]),:);
    
        else

            % set position, velocity and acceleration variables
            pos = q;
            vel = qdot;
            acc = qddot;

            % compute observation matrix for the left/right arm
            for j = 1:i
                [Y_ip, ~, ~] = Yphi(model, pos(:,j),vel(:,j), acc(:,j));
                W_ip = [W_ip; Y_ip];
            end

            % make the observation matrix a double matrix
            W_ip = double(W_ip);

            % get the torque vector
            T = double(reshape(torq,[length(pos(:,1))*length(t) 1]));

            % data measurement vector
            data = double([pos vel acc]);
            dataFull = data;

        end
    
    end
    
end
%% Dataprocessing script

% This script processes the raw data. 

% Inputs:

% cfT = cutoff frqueqency for torque
% cfV = cutoff frqueqency for velocity
% cfA = cutoff frqueqency for acceleration
% orderT = order of Butterworth filter for torque
% orderV = order of Butterworth filter for velocity
% orderA = order of Butterworth filter for acceleration

% cutTimeBefore = remove the time before this value
% cutTimeAfter = remove the time after this value
% pointNumber = approximate number of wished data points

% time = raw time data
% pos = raw position data
% vel = raw velocity data
% torque = raw torque data

% Outputs:

% i = actual number of data points
% t = processed time data
% q = processed position data
% qdot = processed velocity data
% qddot = processed acceleration data
% torque = processed torque data

% Author: Daniel Haugk, 2024, University of Michigan

function [i,t,q,qdot,qddot,torq] = dataProcessing(time,pos,vel,torque,cfT,cfV,cfA,cutTimeBefore,cutTimeAfter,orderT,orderV,orderA,pointNumber)

    %empty data for usage of variables
    t = [];
    q = [];

    %clip time after/before cutTime seconds
    clipAfter = cutTimeAfter;
    clipBefore = cutTimeBefore;

    % number of joints
    nj = length(pos(:,1));    

    if clipBefore == 0
    % case if time starts at 0

        % remove time after "cutTimeAfter" seconds
        time(time>clipAfter) = [];

        % number of datapoints
        i = length(time)-1; 
    
        % downsample factor. Approximates the number of points someone wants to keep.
        ds = ceil(i/pointNumber);
    
        % time between 2 measurements. 1/dt is the frequency.
        dt = time(end)-time(end-1);
    
        % approximate acceleration and filter it.
        qddots = centralDifference(vel,dt,cfA,orderA);

        % set the initial acceleration to 0.
        qddots(:,1) = 0;
        
        % measurement frequency in Hz
        fs = 1/dt; 

        % low-pass filter coefficients for torque (T) and velocity (V)
        [b, a] = butter(orderT, cfT/(fs/2), 'low');
        [b2, a2] = butter(orderV, cfV/(fs/2), 'low');
    
        % apply the filter using filtfilt to achieve zero-phase filtering
        for k = 1:nj
    
            torque(k,:) = filtfilt(b, a, torque(k,:));
            vel(k,:) = filtfilt(b2, a2, vel(k,:));
    
        end

        % downsample the data by the factor ds.
        for k = 1:nj
            torq(k,:) = downsample(torque(k,1:i),ds);
            q(k,:) = downsample(pos(k,1:i),ds);
            qdot(k,:) = downsample(vel(k,1:i),ds);
            qddot(k,:) = downsample(qddots(k,1:i),ds);
            t = downsample(time(1:i),ds);
        end
        
    else
    % case if cutTimeBefore > 0

        % remove time after "cutTimeAfter" seconds
        time(time>clipAfter) = [];

        % number of points left after first removal
        i = length(time);
    
        % remove time before "cutTimeBefore" seconds
        time(time<clipBefore) = [];

        % number of points left after second removal
        n = length(time);
    
        % index of the first element after removing both time instances
        j = i - n;    
    
        % downsample factor. Approximates the number of points someone wants to keep.
        ds = ceil(n/pointNumber);
    
        % time between 2 measurements. 1/dt is the frequency.
        dt = time(end)-time(end-1);
    

        % approximate acceleration and filter it.
        qddots = centralDifference(vel,dt,cfA,orderA);

        % set the initial acceleration to 0.
        qddots(:,1) = 0;
        
        % measurement frequency in Hz
        fs = 1/dt; 

        % low-pass filter coefficients for torque (T) and velocity (V)
        [b, a] = butter(orderT, cfT/(fs/2), 'low');
        [b2, a2] = butter(orderV, cfV/(fs/2), 'low');
    
        % apply the filter using filtfilt to achieve zero-phase filtering
        for k = 1:nj
    
            torque(k,:) = filtfilt(b, a, torque(k,:));
            vel(k,:) = filtfilt(b2, a2, vel(k,:));
    
        end
        
        % downsample the data by the factor ds
        for k = 1:nj
            torq(k,:) = downsample(torque(k,j:i-1),ds);
            q(k,:) = downsample(pos(k,j:i-1),ds);
            qdot(k,:) = downsample(vel(k,j:i-1),ds);
            qddot(k,:) = downsample(qddots(k,j:i-1),ds);
            t = downsample(time(1:n),ds);
        end
    end

    % number of downsampled datapoints
    i = length(t);

end
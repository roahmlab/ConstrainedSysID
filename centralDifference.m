%% 4th Order Central difference
% Author: Daniel Haugk, 2024, University of Michigan

function dydt = centralDifference(y,dt,cfA,orderA)

% coefficients for 4th order central difference
CentCoeff = [1/12 -2/3 0 2/3 -1/12]';

% coefficients for 4th order forward difference
ForwCoeff = [-25/12 4 -3 4/3 -1/4]';

% coefficients for 4th order backward difference
BackwCoeff = flip(-ForwCoeff);

% number of datapoints
n = length(y(1,:));

% number of joints
nj = length(y(:,1));

% measurement frequency in Hz
fs = 1/dt;

% Design the filter
[b, a] = butter(orderA, cfA/(fs/2), 'low');

% approximation for each row from 1 to nq, of y seperately
for j = 1:nj
    
    % Apply the filter using filtfilt to achieve zero-phase filtering
    yfilt = filtfilt(b, a, y(j,:));

    for i=1:n
        switch i
            case 1
                % use FORWARD difference here for the first point
                y_ = yfilt(i:i+4);
                dydt(j,i) = y_*ForwCoeff/dt;
            case 2
                % use FORWARD difference here for the second point
                y_ = yfilt(i:i+4);
                dydt(j,i) = y_*ForwCoeff/dt;
            case n-1
                % use BACKWARD difference here for the second last point
                y_ = yfilt(i-4:i);
                dydt(j,i) = y_*BackwCoeff/dt;
            case n
                % use BACKWARD difference here for the last point
                y_ = yfilt(i-4:i);
                dydt(j,i) = y_*BackwCoeff/dt;
            otherwise
                % use CENTRAL difference for every other point
                y_ = yfilt(i-2:i+2);
                dydt(j,i) = y_*CentCoeff/dt;
        end
    end
end

end
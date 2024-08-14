%% SysID data script to gather all relevant data for the identification

% Author: Daniel Haugk, 2024, University of Michigan

%% get necessary data for data processing

% set data paths
datapaths = {'Digit_datap11.txt','Digit_datap21.txt'};

% set the start time for each data file that is used
cutTimeBefore = [5 5];

% set the end time for each data file that is used
cutTimeAfter = [10 14];

% set the wished number of data points
pointNumber = 2100;

% set the cutoff frequency for torque for each data file that is used
cfT = [33 33];

% set the cutoff frequency for velocity for each data file that is used
cfV = [18 18];

% set the cutoff frequency for acceleration for each data file that is used
cfA = [15 15];

% set the Butterworth filter order for torque
orderT = 4;

% set the Butterworth filter order for velocity
orderV = 4;

% set the Butterworth filter order for acceleration
orderA = 4;

%% data processing

% this part is robot specific and has to be changed accordingly! 
% specific parts which need to be changed depending on the robot and general parts which don't need to be changed are marked accordingly

% specify the body part of digit (digit specific)
bodypart = "LeftLeg";

% set variables for all data, if more than 1 data file is used (general part)
i = 0;
t = [];
torq = [];
qdot = [];
qddot = [];
q = [];

% get the data for each file seperately and then connect them (general part)
for k = 1:length(datapaths)

    % run step 1 of digits file (digit specific, this part has to be changed with a new file according to the robot that is used)
    [time,pos,vel,torque,J,na_idx,nu_idx,~,~,~,~] = SysIDDigit(datapaths{k},bodypart,1,[],[],[],[],[],constraintVariant);
    
    % process the obtained raw data (general part)
    [i_,t_,q_,qdot_,qddot_,torq_] = dataProcessing(time,pos,vel,torque,cfT(k),cfV(k),cfA(k),cutTimeBefore(k),cutTimeAfter(k),orderT,orderV,orderA,pointNumber);

    % connect the data from each data file (general part)
    i = i + i_;
    t = [t t_];
    torq = [torq torq_];
    qdot = [qdot qdot_];
    qddot = [qddot qddot_];
    q = [q q_];

end

% run step 2 of digits file (digit specific, this part has to be changed with a new file according to the robot that is used)
[~,~,~,~,~,~,~,W_ip,T,data,dataFull] = SysIDDigit(datapaths,bodypart,2,q,qdot,qddot,torq,i,constraintVariant);
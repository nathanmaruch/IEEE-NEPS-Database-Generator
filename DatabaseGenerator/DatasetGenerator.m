%% IEEE 39-bus New England power system fault simulation for ANN dataset structuring
%% Clearing workspace and global registers.
clear
clear global
clc

%% Simulation stage.

% Load PS data file into workspace.
datane_nathan 

% Identify number of bueses and lines in order to specify the number of
% simulations to be done.
n_bus = length(bus); 
n_line = length(line); 
fb = line(:,1); 
tb = line(:,2); 
num_samp = ceil(sw_con(5,1)/sw_con(1,7));
pre_samp = ceil(sw_con(2,1)/sw_con(1,7)); 
post_samp = num_samp - pre_samp;

% Sweep through all possible cases.
% No-fault simulation.
for l_fault = 1:n_line
    
    datane_nathan; %Load file into empty workspace
    
    f_type = 6;
    
    sw_con(2,2) = line(l_fault,1); %Specify fault bus
    sw_con(2,3) = line(l_fault,2); %Specify far bus
    sw_con(2,6) = f_type; %Specify fault type
    
    save('nedata.mat'); %Save workspace in order to load it into s_simu
    
    s_simu_nathan % Fault simulation script
    
    % Fault labeling
    f_options    = [0:3];
    f_label_pre  = repmat((f_options==3),pre_samp,1); 
    % Clustering the short-circuit faults.
    if f_type == 0 || f_type == 1 || f_type == 2 || f_type == 3
        f_label_post = repmat((f_options==0),post_samp,1);
    else
        f_label_post = repmat((f_options==f_type-3),post_samp,1);
    end
    f_label      = [f_label_pre; f_label_post];
    
    %PMUs are in buses 4, 8, 16, 28, 30 - 39
    v_bus = [bus_v(4,:);
             bus_v(8,:);
             bus_v(16,:);
             bus_v(28,:);
             bus_v(30:39,:)];
    v_bus = v_bus';
    % Noise insertion with SNR = 45 dB.
    v_bus = awgn(v_bus,45,'measured'); 
    v_bus = [v_bus f_label];
    % Downsample by a factor of 5.
    v_bus = downsample(v_bus,5);
    
    i_bus = [ilt(6,:)+ilf(8,:)+ilf(9,:);
             ilt(10,:)+ilt(14,:)+ilf(15,:);
             ilt(24,:)+ilf(25,:)+ilf(26,:)+ilf(27,:)+ilf(28,:);
             ilf(42,:)+ilf(44,:);
             ilt(5,:);
             ilf(46,:);
             ilt(19,:);
             ilt(31,:);
             ilt(33,:);
             ilt(36,:);
             ilt(38,:);
             ilt(40,:);
             ilt(45,:);
             ilt(2,:)+ilt(9,:)];
    i_bus = i_bus';
    i_bus = awgn(i_bus,45,'measured');
    i_bus = [i_bus f_label];
    i_bus = downsample(i_bus,5);
    
    
    filename = strcat('Caso',num2str(l_fault),num2str(f_type));
    save(filename,'v_bus','i_bus');
    
    clearvars -except l_fault f_type num_samp pre_samp post_samp;
end 

clear
clear global
clc

datane_nathan 

n_bus = length(bus); 
n_line = length(line); 
fb = line(:,1); 
tb = line(:,2); 
num_samp = sw_con(5,1)/sw_con(1,7);
pre_samp = ceil(sw_con(2,1)/sw_con(1,7)); 
post_samp = num_samp - pre_samp;

% Fault simulations. 
for b_fault = 1:n_bus
    
    datane_nathan;
    
    count_f = 0;
    count_t = 0;
    
    for numline = 1:n_line
        if ((fb(numline) == b_fault) & (count_f == 0))
            
            count_f = 1 + count_f;
            
            for f_type = 0:5
                
                datane_nathan;
                
                sw_con(2,2) = b_fault;
                sw_con(2,3) = tb(numline);
                sw_con(2,6) = f_type;
                
                save('nedata.mat');
                
                s_simu_nathan;
                
                f_options    = [0:3];
                f_label_pre  = repmat((f_options==3),pre_samp,1); 
                % Clustering the short-circuit faults.
                if f_type == 0 || f_type == 1 || f_type == 2 || f_type == 3
                    f_label_post = repmat((f_options==0),post_samp,1);
                else
                    f_label_post = repmat((f_options==f_type-3),post_samp,1);
                end
                f_label      = [f_label_pre; f_label_post];
                
                v_bus = [bus_v(4,:);
                         bus_v(8,:);
                         bus_v(16,:);
                         bus_v(28,:);
                         bus_v(30:39,:)];
                v_bus = v_bus';
                v_bus = awgn(v_bus,45,'measured');
                v_bus = [v_bus f_label]; 
                v_bus = downsample(v_bus,5);
                
                i_bus = [ilt(6,:)+ilf(8,:)+ilf(9,:);
                         ilt(10,:)+ilt(14,:)+ilf(15,:);
                         ilt(24,:)+ilf(25,:)+ilf(26,:)+ilf(27,:)+ilf(28,:);
                         ilf(42,:)+ilf(44,:);
                         ilt(5,:);
                         ilf(46,:);
                         ilt(19,:);
                         ilt(31,:);
                         ilt(33,:);
                         ilt(36,:);
                         ilt(38,:);
                         ilt(40,:);
                         ilt(45,:);
                         ilt(2,:)+ilt(9,:)];
                i_bus = i_bus';
                i_bus = awgn(i_bus,45,'measured');
                i_bus = [i_bus f_label];
                i_bus = downsample(i_bus,5);
                
                filename = strcat('Caso',num2str(b_fault),num2str(f_type));
                save(filename,'v_bus','i_bus');
                
                clearvars -except f_type numline b_fault n_bus n_line fb tb count_f count_t num_samp pre_samp post_samp;
            end
        elseif ((tb(numline) == b_fault) & (count_t == 0))
            
            count_t = 1 + count_t;
            
            for f_type = 0:5
                
                datane_nathan;
                
                sw_con(2,2) = b_fault;
                sw_con(2,3) = fb(numline);
                sw_con(2,6) = f_type;
                
                save('nedata.mat');
                
                s_simu_nathan;
                
                f_options    = [0:3];
                f_label_pre  = repmat((f_options==3),pre_samp,1); 
                % Clustering the short-circuit faults.
                if f_type == 0 || f_type == 1 || f_type == 2 || f_type == 3
                    f_label_post = repmat((f_options==0),post_samp,1);
                else
                    f_label_post = repmat((f_options==f_type-3),post_samp,1);
                end
                     
                % f_label_post = repmat((f_options==f_type),post_samp,1);
                f_label      = [f_label_pre; f_label_post];
                
                v_bus = [bus_v(4,:);
                         bus_v(8,:);
                         bus_v(16,:);
                         bus_v(28,:);
                         bus_v(30:39,:)];
                v_bus = v_bus';
                v_bus = awgn(v_bus,45,'measured');
                v_bus = [v_bus f_label];
                v_bus = downsample(v_bus,5);
                
                i_bus = [ilt(6,:)+ilf(8,:)+ilf(9,:);
                         ilt(10,:)+ilt(14,:)+ilf(15,:);
                         ilt(24,:)+ilf(25,:)+ilf(26,:)+ilf(27,:)+ilf(28,:);
                         ilf(42,:)+ilf(44,:);
                         ilt(5,:);
                         ilf(46,:);
                         ilt(19,:);
                         ilt(31,:);
                         ilt(33,:);
                         ilt(36,:);
                         ilt(38,:);
                         ilt(40,:);
                         ilt(45,:);
                         ilt(2,:)+ilt(9,:)];
                i_bus = i_bus';
                i_bus = awgn(i_bus,45,'measured');
                i_bus = [i_bus f_label]; 
                i_bus = downsample(i_bus,5);
                
                filename = strcat('Caso',num2str(b_fault),num2str(f_type));
                save(filename,'v_bus','i_bus');
                
                clearvars -except f_type numline b_fault n_bus n_line fb tb count_f count_t num_samp pre_samp post_samp;
            end
        end
    end
end

%% Data organization, Part I - Scrambling the dataset in blocks of W samples. 
clear
clear global
clc

delete('sim_fle.mat');
delete('nedata.mat');

mat = dir('Caso*.mat'); %Load all .mat files corresponding to all simulations

for q = 1:length(mat) %Structure the variable matrixes so that all simulations fit into a single (sample numbers)x15 matrix
    
    if q == 1
        
        Q = load(mat(q).name);
        v_bus = Q.v_bus; 
        i_bus = Q.i_bus;
        
    else
        K = load(mat(q).name); 
        v_bus = [v_bus; K.v_bus];
        i_bus = [i_bus; K.i_bus];
    end
    
end

clear K mat q Q; 

v_b = v_bus(:,(1:14)); %Values-only matrix
i_b = i_bus(:,(1:14));
v_idx = v_bus(:,(15:18)); %Fault type array
i_idx = i_bus(:,(15:18)); 

%Calculate the frequency of the PMU buses. Nominal frequency is 1 pu.  
mod_v = abs(v_b); 
ang_v = angle(v_b);
mod_i = abs(i_b); 
ang_i = angle(i_b);

freq_size = size(ang_v); 
freq = ones(freq_size); 
freq(2:end,:) = 1 + diff(ang_v); % The first value is the unchanged nominal frequency. 

%Convert rectangular to polar form so that values are more appropriate to a
%[0 1] interval upon normalization
mod_v = (mod_v - min(mod_v,[],'all'))/(max(mod_v,[],'all') - min(mod_v,[],'all'));
ang_v = (ang_v - min(ang_v,[],'all'))/(max(ang_v,[],'all') - min(ang_v,[],'all'));

freq =  (freq - min(freq,[],'all'))/(max(freq,[],'all') - min(freq,[],'all'));
 
mod_i = (mod_i - min(mod_v,[],'all'))/(max(mod_i,[],'all') - min(mod_i,[],'all'));
ang_i = (ang_i - min(ang_i,[],'all'))/(max(ang_i,[],'all') - min(ang_i,[],'all'));


clear i_b v_b v_bus i_bus freq_size

% PMU Selection
% Edit the pmu_idx array to select to PMUs based on their index in the
% voltage, frequency and current arrays.
% Comment to disable this feature. 

% pmu_idx = [1 2 3 4 6 8 14]; %Keeping 4, 8, 16, 28, 31, 33 and 39. 
% 
% mod_v = mod_v(:,pmu_idx);
% ang_v = ang_v(:,pmu_idx);
% freq = freq(:,pmu_idx);
% mod_i = mod_i(:,pmu_idx);
% ang_i = ang_i(:,pmu_idx);

% Dataset without previous instant in consideration.
dataset = [mod_v ang_v mod_i ang_i freq v_idx];

% Multiplying line and load contingency events
% datamod = dataset;

% H = 3; %Constant for dataset expansion
% for l = 1:H
%     for k = 1:length(dataset)
%         if(v_idx(k,2) == 1 || v_idx(k,3) == 1 || (v_idx(k,4) == 1 && l > 1))
%             datamod = [datamod; datamod(k,:)];
%         end
%     end
% end
% 
% dataset = datamod;

%Shuffle the matrix lines in windows of W samples.
%This is done so that even though the data order is randomized, the samples
%coincide both in timestamp and fault label.
W = 1; 
col_num = size(dataset,2);
dataset = randblock(dataset, [W col_num]);
filename = strcat(num2str(W),'SampleWindow_Freq_Downsampled_Results');
save(filename, 'dataset');

clear i_idx v_idx W col_num datamod
%% Data organization, Part II - Separating it as inputs and targets. 

inputs  = dataset(:,1:70); 
targets = dataset(:,71:74);

filename = 'InputLabel_PMU_IEEE39';
save(filename, 'inputs', 'targets');

clear ang_i ang_v mod_i mod_v dataset freq
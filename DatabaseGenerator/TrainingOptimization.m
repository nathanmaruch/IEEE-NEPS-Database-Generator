clear
clc

% Load database and set inputs and targets for the network. 
load('best_net.mat');

clearvars -except inputs targets outputs best_idx  

prec = zeros(50,4);
rec = zeros(50,4); 
f1 = zeros(50,4); 
f1g = zeros(50,1); 
acc = zeros(50,1);

perf = zeros(50,1); 
best_perf = 1;

% Create ANN for Pattern Recognition
hiddenLayerSize = best_idx; 
for k = 1:50
    
    net = patternnet(hiddenLayerSize); 

    % Set up division of data for Training, Validation and Testing
    net.divideParam.trainRatio = 70/100; 
    net.divideParam.valRatio = 15/100; 
    net.divideParam.testRatio = 15/100; 

    % Train the network 
    [net,tr] = train(net,inputs,targets); 
    
    outputs = net(inputs); 
    
    %Metrics
    [c, cm, ind, per] = confusion(targets,outputs); 
    
    confmat = cm';
    
    for m = 1:4
    
        prec(k,m) = confmat(m,m)/(sum(confmat(m,:)));
        rec(k,m) = confmat(m,m)/(sum(confmat(:,m)));
        f1(k,m) = 2*(prec(k,m)*rec(k,m))/(prec(k,m)+rec(k,m));
        
    end
    
    acc(k,1) = trace(confmat)/sum(confmat, 'all'); 
    f1g(k,1) = mean(f1(:,1));  
    
    % Get network performance metric
    perf(k,1) = perform(net,targets,outputs);
    
    if perf(k,1) < best_perf
        best_net = net;
        best_idx = k;
        best_perf = perf(k,1); 
    end
    
end

delete tgt

save perf; 
save best_net; 
save best_idx;  








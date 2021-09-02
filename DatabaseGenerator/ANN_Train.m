clear
clc

load('InputLabel_PMU_IEEE39_14PMUwn.mat'); 

inputs = inputs';
targets = targets';

best_perf = 1;

hiddenLayerSize = zeros(1,2);

for l = 1:20
    for k = 0:20
        
        if k == 0
            hiddenLayerSize = l; 
        else
            hiddenLayerSize = [l k];
        end
                
        net = patternnet(hiddenLayerSize);
        
        net.divideParam.trainRatio = 70/100; 
        net.divideParam.valRatio = 15/100; 
        net.divideParam.testRatio = 15/100; 
        
        [net, tr] = train(net,inputs,targets); 
        outputs = net(inputs); 
        
        perf = perform(net,targets,outputs); 
        
        if perf < best_perf
            best_net = net; 
            best_idx = [l k]; 
            best_perf = perf;
        end
    end
end

save perf; 
save best_net; 
save best_idx; 

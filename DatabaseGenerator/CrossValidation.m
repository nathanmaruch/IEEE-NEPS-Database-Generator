clear
clc

load('best_net.mat');

clearvars -except inputs targets outputs best_idx hiddenLayerSize

N = size(inputs,2);

n = size(outputs,1);

k = 10; 

cvFolds = 1 + mod(randperm(N),k);

prec = zeros(k,n);
rec = zeros(k,n); 
f1 = zeros(k,n); 
f1g = zeros(k,1); 
acc = zeros(k,1);

for i = 1:k
    
    testIdx = (cvFolds == i); 
    trainIdx = ~testIdx;  
    
    trInd = find(trainIdx); 
    tstInd = find(testIdx); 
    
    perf = zeros(10,1);
    best_perf = 1; 
    
    for q = 1:5
        
        net = patternnet(hiddenLayerSize); 
        
        net.divideParam.trainInd = trInd; 
        net.divideParam.testInd = tstInd; 
        
        [net, tr] = train(net,inputs,targets); 
        
        outputs = net(inputs);
        
        perf(q,1) = perform(net,targets,outputs);
        
        if perf(q,1) < best_perf
            
            best_net = net;
            best_idx = q;
            best_perf = perf(q,1); 
            
            [c, cm, ind, per] = confusion(targets,outputs); 
            
        end
    end
       
    confmat = cm';
    
    for m = 1:n
    
        prec(i,m) = confmat(m,m)/(sum(confmat(m,:)));
        rec(i,m) = confmat(m,m)/(sum(confmat(:,m)));
        f1(i,m) = 2*(prec(i,m)*rec(i,m))/(prec(i,m)+rec(i,m));
        
    end
    
    acc(i,1) = trace(confmat)/sum(confmat, 'all'); 
    f1g(i,1) = mean(f1(1,:));  
    
end

total_acc = mean(acc); 
total_prec = mean(prec); 
total_rec = mean(rec); 
total_f1 = mean(f1); 
total_f1g = mean(f1g); 

clearvars -except total_acc total_prec total_rec total_f1 prec rec f1 acc f1g total_f1g;
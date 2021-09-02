clear 
clc 

load('best_net.mat')

% Training Confusion Matrix 
outTrn = best_net(inputs(:,tr.trainInd)); 
tgtTrn = targets(:,tr.trainInd);

% Validation Confusion Matrix 
outVal = best_net(inputs(:,tr.valInd)); 
tgtVal = targets(:,tr.valInd);

% Testing Confusion Matrix 
outTst = best_net(inputs(:,tr.testInd)); 
tgtTst = targets(:,tr.testInd);

% Overall Confusion Matrix
out = best_net(inputs); 
tgt = targets; 
plotconfusion(tgtTrn, outTrn, 'Training', tgtVal, outVal, 'Validation', tgtTst, outTst, 'Testing', tgt,out,'Overall')


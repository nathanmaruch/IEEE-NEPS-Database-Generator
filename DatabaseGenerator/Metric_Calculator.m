clear 
clc

load('D:\Google Drive\Nathan Maruch - PIBIC 2020\5. IAA - Localização de Faltas\ann_results.mat');

tgt_val = results.target; 
output = results.output; 

k = length(output);
pred_val = zeros(2,k);

for q = 1:k
    if(output(1,q) > 0.5)
        pred_val(1,q) = 1;
    else
        if(output(1,q) < 0.5)
            pred_val(2,q) = 1;
        end
    end
end

accuracy = mean(double(pred_val == tgt_val));
acc_all = mean(double(0 == tgt_val));

actual_pos = sum(tgt_val == 1);
actual_neg = sum(tgt_val == 0);

true_pos = sum((pred_val == 1) & (tgt_val == 1));
false_pos = sum((pred_val == 1) & (tgt_val == 0));
true_neg = sum((pred_val == 0) & (tgt_val == 0));
false_neg = sum((pred_val == 0) & (tgt_val == 1));

precision = 0;
if ((true_pos + false_pos) > 0)
    precision = true_pos/(true_pos + false_pos);
end

recall = 0;
if((true_pos + false_neg) > 0)
    recall = true_pos/(true_pos + false_neg);
end

F1 = 0;
if((precision + recall) > 0)
    F1 = 2*((precision*recall)/(precision + recall));
end











function [WW_bvp,AA_acc]=windetection(BVP,Acc,W_bvp,overlap)
%% BVP: bvp data trials read directly from the input csv file 
%% Acc: accel data trials read directly from the input csv file (both should have the same fs and they should be synchronized)
%% W_bvp: number of windows that you want to divide the dataset with a 50% of overlap per each window.
%% overlap: percentage of overlap (if overlap==1 is consired it won't be applied any overlap)
%% WW_bvp and AA_acc are the cell arrays with the resulting window segments
num_samples=floor(length(BVP)/W_bvp);
if (overlap~=1)
    for(i=1:W_bvp)
      if (round(num_samples*(i-1)-(num_samples*(i-1)/(1/overlap)))<=length(BVP) && round(num_samples*(i)+num_samples*(i)/(1/overlap))<=length(Acc))
        WW_bvp{i}=BVP(round(num_samples*(i-1)-(num_samples*(i-1)/(1/overlap)))+1:round(num_samples*(i)+(num_samples*(i)/(1/overlap))));
        AA_acc{i}=Acc(round(num_samples*(i-1)-(num_samples*(i-1)/(1/overlap)))+1:round(num_samples*(i)+(num_samples*(i)/(1/overlap))));
      else
        WW_bvp{i}=BVP(round(num_samples*(i-1)-(num_samples*(i-1)/(1/overlap)))+1:length(BVP));
        AA_acc{i}=Acc(round(num_samples*(i-1)-(num_samples*(i-1)/(1/overlap)))+1:length(Acc));  
      end;
    end;
else
    for(i=1:W_bvp)
      if (round(num_samples*(i-1))<=length(BVP) && round(num_samples*(i))<=length(Acc))
        WW_bvp{i}=BVP(round(num_samples*(i-1))+1:round(num_samples*(i)));
        AA_acc{i}=Acc(round(num_samples*(i-1))+1:round(num_samples*(i)));
      else
        WW_bvp{i}=BVP(round(num_samples*(i-1))+1:length(BVP));
        AA_acc{i}=Acc(round(num_samples*(i-1))+1:length(Acc));  
      end;
    end;
end;

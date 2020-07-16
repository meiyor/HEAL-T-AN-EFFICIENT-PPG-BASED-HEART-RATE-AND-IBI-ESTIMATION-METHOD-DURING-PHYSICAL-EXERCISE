function [WW,AA]=windetection(BVP,Acc,W)
%% data windowing for HEAL-T second method
num_samples=round(length(BVP)/W);
for(i=1:W)
    WW{i}=BVP(num_samples*(i-1)+1:num_samples*(i));
    AA{i}=Acc(round(num_samples*(i-1)/2)+1:round(num_samples*(i)/2));
end;

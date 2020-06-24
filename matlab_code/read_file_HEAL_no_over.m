function read_file_HEAL_no_over(inpath,outpath,subt,tses,type,tamt,W,seld,selk,fss,ini_filt,inc,overlap,s_sel)
%% subt:subject index that you want to explore inside /Data folder
%% tses: session that you desire to characterize per subject.
%% type: is a cell that represent the prefix of each included per session.
%% fs: is a vector of the same size of type cell, that represent the sampling frequency defined by the Empatica handheld
%% in: is the factor of downsampling that should be taking into account to calculate the interpolated IBI.
%% p_inter: this is the smoothing factor to interpolate the estimated IBI (between 0 and 1).
%% tam: User defined size for the current simulation (source BVP signal, PPG).
%% W: the number of windows to select the less variable signal segments from PPGsignal.
%% seld: selector for the input of the complete size of your dataset,
%% selk: selector (0) add the thirdparty path for (eeglab especially) (1) remove thirdparty path 
%% fss: sampling frequency for the bvp input data (Hz)
%% ini_filt: initial filter bandwidth defined by the user
%% inc: constant increment for BHW bandwidth defined by the user
%% overlap: percentage of overlap defined by the user.
%% s_sel: selector for type of spline 0 -> cubic spline , 1 -> moving average filter.
close all;
addpath(genpath([pwd '/functions_HR/' ]));
P=0;
if (seld==0)
 tam=tamt*60*60*fss; %% select hours time delimitator;
else
 tam=tamt;
end;
ttp(1)=0;
h1=figure(1);
for (j=1:length(type))
  temp=load([inpath '/' subt '/' tses '/' type{j} '.txt']);
  time{j}=temp(:,1)-temp(1,1);
if (strcmp(type{j},'ibi') || strcmp(type{j},'gsr') || strcmp(type{j},'bvp') || strcmp(type{j},'temp'))
  D{j}=temp(:,2);
 else
   D{j}=temp(:,2:4);  
end;
subplot(length(type),1,j);
plot(time{j}(1:length(time{j})-1),D{j}(1:length(time{j})-1,:));
grid on;
end;
time{1}=downsample(time{1},2);
D{2}=sqrt(D{2}(:,1).^2+D{2}(:,2).^2+D{2}(:,3).^2)';
if (length(D{1})>=1000)%% assure that the data is big enough to process
[WW,AA]=windetection(downsample(D{1},2),D{2},W,overlap);
[IBI1,IBI2,T1,T2,BVPn]=eval_winHR_no_over(WW,AA,W,time,ini_filt,inc,0.6054,selk,fss,0,s_sel);
if (exist([ outpath '/'])==0)
    mkdir([ outpath '/']);
end;
pid=feature('getpid');
if (selk==0)
   f1=fopen([ outpath '/' num2str(pid)  '_ibi.txt' ],'w');
   f2=fopen([ outpath '/' num2str(pid)  '_iibi.txt'],'w');
   f3=fopen([ outpath '/' num2str(pid) '_bvpn.txt'],'w'); 
else
 f1=fopen([ outpath '/' num2str(pid)  '_ibibase.txt'],'w');
 f2=fopen([ outpath '/' num2str(pid)  '_iibibase.txt'],'w');                                                                  
 f3=fopen([ outpath '/' num2str(pid) '_bvpnbase.txt'],'w');  
end;
if (length(IBI1)<=length(T1))  
  for (k=1:length(IBI1))
      fprintf(f1,'%f,%f\n',[T1(k) IBI1(k)]);
   end;
  fclose(f1);
  for (k=1:length(IBI2))
     fprintf(f2,'%f\n',IBI2(k));
  end;
  fclose(f2);
 for (k=1:length(BVPn))
    fprintf(f3,'%f,%f\n',[T1(k) BVPn(k)]);
 end;
 fclose(f3);
else
   for (k=1:length(T1))                                                                                                                                                         
      fprintf(f1,'%f,%f\n',[T1(k) IBI1(k)]);                                                                                                                                      
   end;                                                                                                                                                                           
  fclose(f1);                                                                                                                                                                     
  for (k=1:length(T2))                                                                                                                                                          
     fprintf(f2,'%f\n',IBI2(k));                                                                                                                                      
  end;                                                                                                                                                                            
  fclose(f2);
  for (k=1:length(T1))
    fprintf(f3,'%f,%f\n',[T1(k) BVPn(k)]);
  end;
  fclose(f3);   
end;
end;
disp(['Process ' num2str(pid) ' has ended succesfully']);

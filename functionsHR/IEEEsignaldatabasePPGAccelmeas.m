function [tt]=IEEEsignaldatabasePPGAccelmeas(S,type,fss,seld,baseo,oversel,sel_tr_t,Wk,selp)
% S: Subject index
%type: exercise type routine
% fss: is the sampling frequency for this TROIKA data fssHz
% seld: this selector is giving to eval certain path nomenclature
% baseo: this is to evaluate the baseline and HEAL-T method 0-> HEAL-T
% method, 1 -> baseline
% oversel: this parameter controls when the approach is giving by averaging
% or by continuous sampling using ECG. 0 -> averaging, 1-> continous
% sel_tr_t: chooses the training 0, or the test set 1 from the 
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.6178;
opts.Normalize = 'on';
close all
addpath(genpath([pwd '/Training_data']));
addpath(genpath([pwd '/Extra_TrainingData']));
addpath(genpath([pwd '/TestData']));
addpath(genpath([pwd '/TrueBPM']));
if(baseo==0)
  addpath(genpath([pwd '/thirdparty/eeglab12_0_2_5b']));
else
  rmpath(genpath([pwd '/thirdparty/eeglab12_0_2_5b']));  
end;
if (seld==0)
if (S>=0 && S<=9)
 if(sel_tr_t==0)   
    F=load(['DATA_0' num2str(S) '_TYPE0' num2str(type) '.mat']);
    Q=load(['DATA_0' num2str(S) '_TYPE0' num2str(type) '_BPMtrace.mat']);
 else
    F=load(['TEST_S0' num2str(S) '_T0' num2str(type) '.mat']);
    Q=load(['True_S0' num2str(S) '_T0' num2str(type) '.mat']);
 end;
else
  F=load(['DATA_' num2str(S) '_TYPE0' num2str(type) '.mat']);
  Q=load(['DATA_' num2str(S) '_TYPE0' num2str(type) '_BPMtrace.mat']);
end;
else
  F=load(['DATA_S0' num2str(S) '_T0' num2str(type) '.mat']);
  Q=load(['BPM_S0' num2str(S) '_T0' num2str(type) '.mat']);
end;
if (sel_tr_t==0)
    [HRVcont,Tdata]=calcECGdistdef(F.sig(1,:),fss);
end;
if (sel_tr_t==1)
    F.sig(2:6,:)=F.sig(1:5,:);
end;
PPG1=F.sig(2,:);
PPG2=F.sig(3,:);
Accel=F.sig(4:6,:);
Data=[PPG1 ; PPG2 ; F.sig(4:6,:)];
tic;
if (baseo==0)
   % [As,Ws]=runica(Data,'sphering','off');
   % Datanew=icaproj(Data,As*Ws,[1 2]);
   % PPG1=Datanew(1,:);
   % PPG2=Datanew(2,:);
   % Accel=sqrt(Accel(1,:).^2+Accel(2,:).^2+Accel(3,:).^2);
else
    PPG1=Data(1,:);
    PPG2=Data(2,:);
    Accel=F.sig(4,:); %% Kim States that following a cylindrical coordinates only x-axis are related really with motion
end;
%[PPG1,r,vr]=ssa(PPG1,10);
%[PPG2,r,vr]=ssa(PPG2,10);
figure(1);
plot(linspace(0,(1/fss)*length(F.sig),length(F.sig)),F.sig(2:3,:));
grid on
figure(2);
plot(linspace(0,(1/fss)*length(F.sig),length(F.sig)),F.sig(4:6,:));
grid on
figure(3)
plot(Q.BPM0);
grid on;
%PPG1=downsample(PPG1,10);
%PPG2=downsample(PPG2,10);
%Accel=downsample(PPG2,10);
if (oversel==0)
         [W,Acnew,timeW,sW,P]=detectcepstwin(PPG1,Accel,fss,Wk,floor(Wk/3));
elseif (oversel==2)
         [W,Acnew,timeW,sW,P]=detectcepstwin(PPG1,Accel,fss,Wk,floor(Wk/3));
          Wn=detectcepstwin(F.sig(1,:),Accel,fss,Wk,floor(Wk/3));
else
       [W,Acnew,timeW,sW,P]=detectcepstwin(PPG1,Accel,fss,Wk,Wk/3);
end;
profile on;
if (baseo==0 && selp==0)
 [As,Ws]=runica([W{1};Acnew{1}],'sphering','off');
 Datanew=icaproj([W{1};Acnew{1}],As*Ws,[1]);
 W{1}=Datanew(1,:);
 Acnew{1}=Datanew(2,:);
end;
%profile on;
%[Htest,llg]=evalwinHR(W,Acnew,P-1,Q,timeW,sW,0,0.5978,baseo,fss,oversel,selp);
tt=evalwinmeas(W,Acnew,P-1,Q,timeW,sW,0,0.5978,baseo,fss,oversel,selp);
%profile viewer
%p = profile('info');
%profsave(p,'profile_results')
% temp=[];
% if (oversel==2)
% for (k=1:length(W))
%     [SW{k},T]=calcECGdistdef(Wn{k},fss);
%     if (k>=2)
%          temp=[temp SW{k-1}(1:249) (SW{k-1}(249:249+749)+SW{k}(1:750))/2]; %% HR is already defined inside this function
%     end;
% end; 
%     HRVcont=conv(temp,(1/900).*ones([1 900])); %% this is assignment is only used when you want to use an ovelap evaluation
%     Tdata=linspace(0,(1/fss)*length(PPG1),length(HRVcont));
% end;
% if (oversel==1 || oversel==2)
%    [fitresult, gof]=fit([1:1:llg]',Htest(700/2:llg+700/2-1)',ft,opts);
%     smoothpp=fitresult(1:llg);
%     figure(5)
%     plot(Tdata(1:llg),Htest(700/2:llg+700/2-1));
%     hold on;
%     plot(Tdata(1:llg),HRVcont(1:llg),'r');
%     plot(Tdata(1:llg),smoothpp,'g');
%     grid on;
%     xlabel('Time [s]');
%     ylabel('HRV [BPM]');
%     Error1=mean(abs(Htest(700/2:llg+700/2-1)-HRVcont(900/2:llg+900/2-1)))
%     Error2=mean(abs(Htest(700/2:llg+700/2-1)-HRVcont(900/2:llg+900/2-1))./HRVcont(900/2:llg+900/2-1))
%     Error3=mean(abs(smoothpp'-HRVcont(900/2:llg+900/2-1)))
%     Error4=mean(abs(smoothpp'-HRVcont(900/2:llg+900/2-1))./HRVcont(900/2:llg+900/2-1))
% end;

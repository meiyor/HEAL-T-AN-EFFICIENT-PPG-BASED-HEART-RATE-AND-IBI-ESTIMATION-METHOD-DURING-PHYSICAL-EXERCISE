function [Htest,llg,sBVP,TTBVP,smpp]=IEEEsignaldatabasePPGAccel(S,type,fss,seld,baseo,oversel,sel_tr_t,selp,selq,ini_filt,inc_filt,W_num,overlap)
% S: Subject index
%type: exercise type routine
% fss: is the sampling frequency for this TROIKA data fssHz
% seld: this selector is giving to eval certain path nomenclature
% baseo: this is to evaluate the baseline and HEAL-T method 0-> HEAL-T
% method, 1 -> baseline
% oversel: this parameter controls when the approach is giving by averaging
% or by continuous sampling using ECG. 0 -> averaging, 1-> continous
% sel_tr_t: chooses the training 0, or the test set 1 from the 
% selq : this is to select the evaluation of TROIKA dataset or the evaluation of HR based on Arindam's data
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
if (selq==0)
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
else
 if (length(S)<=2)   
    PPG1=S{1}./max(S{1});
    PPG2=S{1}./max(S{1});
    S{2}(:,1)=S{2}(:,1)./max(S{2}(:,1));
    S{2}(:,1)=S{2}(:,1)-mean(S{2}(:,1));
    S{2}(:,2)=S{2}(:,2)./max(S{2}(:,2));                                                         
    S{2}(:,2)=S{2}(:,2)-mean(S{2}(:,2));
    S{2}(:,3)=S{2}(:,3)./max(S{2}(:,3));                                                         
    S{2}(:,3)=S{2}(:,3)-mean(S{2}(:,3));
    Accel=S{2};
 else %% only if your experiment has two PPG channels
    PPG1=S{1}./max(S{1});
    PPG2=S{2}./max(S{2}); 
    S{3}(:,1)=S{3}(:,1)./max(S{3}(:,1));
    S{3}(:,1)=S{3}(:,1)-mean(S{3}(:,1));
    S{3}(:,2)=S{3}(:,2)./max(S{3}(:,2));                                                         
    S{3}(:,2)=S{3}(:,2)-mean(S{3}(:,2));
    S{3}(:,3)=S{3}(:,3)./max(S{3}(:,3));                                                         
    S{3}(:,3)=S{3}(:,3)-mean(S{3}(:,3));
    Accel=S{3};
 end;
 if (length(PPG1)>=length(Accel))
   PPG1=PPG1(1:length(PPG1)-(length(PPG1)-length(Accel))); %./max(PPG1(1:length(PPG1)-(length(PPG1)-length(Accel))));
   PPG2=[diff(PPG2(1:length(PPG2)-(length(PPG2)-length(Accel)))) ; 0]; %./max(PPG1(1:length(PPG1)-(length(PPG1)-length(Accel))))) ; 100]; 
else
  Accel=Accel(1:length(Accel)-(length(Accel)-length(PPG1)),:);
 end;
 Data=[PPG1' ; PPG2' ; Accel'];  
 Accel=Accel';
end;
tic;
if (baseo==0)
   if (selq==0)
    [As,Ws]=runica(Data,'sphering','off');
    val_s=std(Data');
    h_s=fliplr(sort(val_s));
    Datanew=icaproj(Data,As*Ws,[find(val_s==h_s(1)) find(val_s==h_s(2))]);
    PPG1=Datanew(1,:);
    PPG2=Datanew(2,:);
   else
    if (length(S)<=2)   
        [As,Ws]=runica(Data(2:end,:),'sphering','off');                                                                                                                                    
        Datanew=icaproj(Data(2:end,:),As*Ws,[1]);
    else
        [As,Ws]=runica(Data(:,:),'sphering','off');                                                                                                                                    
        Datanew=icaproj(Data(:,:),As*Ws,[1 2]);
    end;
    PPG1=Datanew(1,:);                 
   end;
    Accel=sqrt(Accel(1,:).^2+Accel(2,:).^2+Accel(3,:).^2);
else
    PPG1=Data(1,:);
    PPG2=Data(2,:);
    Accel=F.sig(4,:); %% Kim States that following a cylindrical coordinates only x-axis are related really with motion
end;

figure(1);
if (selq==0)
   plot(linspace(0,(1/fss)*length(F.sig),length(F.sig)),F.sig(2:3,:));
   grid on
   figure(2);
   plot(linspace(0,(1/fss)*length(F.sig),length(F.sig)),F.sig(4:6,:));
   grid on
   figure(3)
   plot(Q.BPM0); 
   grid on;
else
  plot(linspace(0,(1/fss)*length(PPG1),length(PPG1)),PPG1);
   grid on
   figure(2);
   plot(linspace(0,(1/fss)*length(Accel),length(Accel)),Accel);
   grid on
end;

if (oversel==0)
        % [W,Acnew,timeW,sW,P]=detectcepstwin2(Data,fss,100,50);
      if (selq==0)
         [W,Acnew,timeW,sW,P]=detectcepstwin(PPG1,Accel,fss,W_num,overlap); %% simple detect based on the prelabeling and where's the current window located
      else  
	[W,Acnew,timeW,sW,P,sal]=detectcepstwin3(real(PPG1),Accel,fss,W_num,overlap);
        if (sal==1)
           error('Error: Window size is greater than the entire trial size!!');
           return;
        end;
      end;
     elseif (oversel==2)
      if (selq==0)
	    [W,Acnew,timeW,sW,P]=detectcepstwin(real(PPG1),Accel,fss,W_num,overlap);
         Wn=detectcepstwin(F.sig(1,:),Accel,fss,8,6);
      else
	[W,Acnew,timeW,sW,P,sal]=detectcepstwin3(real(PPG1),Accel,fss,W_num,overlap);
        if (sal==1)
	  error('Error: Window size is greater than the entire trial size!!');
          return;
        end;
       end;
else
  [W,Acnew,timeW,sW,P]=detectcepstwin(real(PPG1),Accel,fss,W_num,overlap);
end;
if (selq==1)
    Q.BPM0=ones([length(W) 1]);
end;
%profile on;
if (selq==0)
   [Htest,llg]=evalwinHR(W,Acnew,P-1,Q,timeW,sW,0,0.5978,baseo,fss,oversel,selp);
  else
   [Htest,llg,sBVP,TTBVP,smpp]=evalwinHR3(W,Acnew,P-1,Q,ini_filt,inc_filt,sW,0,0.5978,baseo,fss,oversel,selp,W_num,overlap);
 end;
%tt=evalwinmeas(W,Acnew,P-1,Q,timeW,sW,0,0.5978,baseo,fss,oversel,selp);
%profile viewer
%p = profile('info');
%profsave(p,'profile_results')
temp=[];
if (oversel==2 && selq==0)
for (k=1:length(W))
    [SW{k},T]=calcECGdistdef(Wn{k},fss);
     if (k>=2)
       temp=[temp  (SW{k-1}(249:249+749)+SW{k}(1:750))/2]; %% HR is already defined inside this function
      else
       temp=SW{1}(1:249);
      end;
end; 
    HRVcont=conv(temp,(1/900).*ones([1 900])); %% this is assignment is only used when you want to use an ovelap evaluation
    Tdata=linspace(0,(1/fss)*length(PPG1),length(HRVcont));
end;
if (oversel==1 || oversel==2 && selq==0)
   [fitresult, gof]=fit([1:1:llg]',Htest(700/2:llg+700/2-1)',ft,opts);
    smoothpp=fitresult(1:llg);
    figure(5)
    plot(Tdata(1:llg),Htest(700/2:llg+700/2-1));
    hold on;
    plot(Tdata(1:llg),HRVcont(1:llg),'r');
    plot(Tdata(1:llg),smoothpp,'g');
    grid on;
    xlabel('Time [s]');
    ylabel('HRV [BPM]');
    Error1=mean(abs(Htest(700/2:llg+700/2-1)-HRVcont(900/2:llg+900/2-1)))
    Error2=mean(abs(Htest(700/2:llg+700/2-1)-HRVcont(900/2:llg+900/2-1))./HRVcont(900/2:llg+900/2-1))
    Error3=mean(abs(smoothpp'-HRVcont(900/2:llg+900/2-1)))
    Error4=mean(abs(smoothpp'-HRVcont(900/2:llg+900/2-1))./HRVcont(900/2:llg+900/2-1))
end;

function readfileHEAL(subt,session,type,fs,in,p_inter,tam,W,namesb,Mc,varc)
%% subt:subject index that you want to explore inside /Data folder
%% session: session that you desire to characterize per subject.
%% type: is a cell that represent the prefix of each included per session.
%% fs: is a vector of the same size of type cell, that represent the sampling frequency defined by the Emphatica handheld
%% in: is the factor of downsampling that should be taking into account to calculate the interpolated IBI.
%% p_inter: this is the smoothing factor to interpolate the estimated IBI (between 0 and 1).
%% tam: User defined size for the current simulation (source BVP signal, PPG).
%% W: the number of windows to select the less variable signal segments from PPGsignal.
%% Mc: order for moving average filter stage
close all;
addpath(genpath([pwd '/eeglab12_0_2_5b']));
P=0;
ttp(1)=0; 
%% Read the session File to Synchronize the start of the data session 
%% per each stream
f=fopen([pwd '/Data/' subt '/' session '/session.txt'],'r');
temp=fscanf(f,'%s',1);
start=str2num(temp(1:17));
temp=fscanf(f,'%s',1);
kk=temp(1:17);
ennd=str2num(temp(1:17));
fclose(f);
h1=figure(1);
for (j=1:length(type))
temp=0;
if(j>=2)
    fclose(afile);
end;
i=1;
S=kk;
afile=fopen([pwd '/Data/' subt '/' session '/' type{j} '.txt'],'r');
%% Only 1 channel data on session 
if (strcmp(type{j},'ibi') || strcmp(type{j},'gsr') || strcmp(type{j},'bvp') || strcmp(type{j},'temp'))
   while(~feof(afile) && i<=tam && str2num(S(1:17))-ennd<=0)
       S=fscanf(afile,'%s',1);
       time{j}(i)=str2num(S(1:17))-start;
       D{j}(i)=str2num(S(19:length(S)));
       temp=time{j}(i);
       i=i+1;
   end;
elseif (strcmp(type{j},'acc')) %% if it is accel divided it in 3 file positions
    while(~feof(afile) && i<=tam/2 && str2num(S(1:17))-ennd<=0)
        S=fscanf(afile,'%s',1);
        time{j}(i)=str2num(S(1:17))-start; %% all times relative with session starting time
        temk=S(18:length(S));
        cmm=find(temk==',');
        D{j}(i,1)=str2num(temk(cmm(1)+1:cmm(2)-1)); % X axis
        D{j}(i,2)=str2num(temk(cmm(2)+1:cmm(3)-1)); % Y axis
        D{j}(i,3)=str2num(temk(cmm(3)+1:length(temk))); % Z axis
        i=i+1;
    end;
end;
subplot(length(type),1,j);
if (type{j}=='ibi')
    yy=interp1(time{j},D{j},linspace(min(time{j}),max(time{j}),1000));
     plot(linspace(min(time{j}),max(time{j}),1000),yy,'r'); %% plot each data accordingly with data asked on input parameters
     hold on 
     plot(time{j}(1:length(time{j})-1),D{j}(1:length(time{j})-1));
     if (j==1)
        P=[min(time{j}),max(time{j})]; %% when you want to synch with empatica ibi select ibi as first ibi
     end;
 else
  if (length(P)~=2) 
      if (strcmp(type{j},'acc')) 
        plot(time{j}(1:length(time{j})-1),D{j}(1:length(time{j})-1,:));  
      else
        plot(time{j}(1:length(time{j})-1),D{j}(1:length(time{j})-1));
      end;
   else
        pos1=max(find(time{j}<=P(1))); %% plot from synchronized time with empatica IBI for any signal
        pos2=min(find(time{j}>=P(2)));
        plot(time{j}(pos1:pos2),D{j}(pos1:pos2));
   end;
end;
grid on;
xlabel('Time [s]');
saveas(h1,[pwd '/folder_images/timesignals.fig']);
end;
j=1; %% fix the signal index to capture BVP only (always called BVP signal as first index)
h2=figure(2);
F{j}=abs(fft(D{j})); %% Spectrum of th BVP signal 
plot(linspace(0,fs(j)/2,length(F{j})/2),F{j}(1:length(F{j})/2));
grid on;
saveas(h2,[pwd '/folder_images/' namesb{1} '.fig']);
h3=figure(3);
% Filters and Movement removal
T=conv(D{j},(1/Mc)*ones([1 Mc])); %% mean average filter
Tpp=T;
T=compacc(T(round(Mc/2):tam+round(Mc/2)-1),D{j+1},W,varc); %% selection Filter based on accelerometer magnitude
%T=compacc(D{j},D{j+1},W); %% selection Filter based on accelerometer magnitude
%T=conv(T,(1/80)*ones([1 80])); %% mean average filter
%T=T(40:tam+39);
Tk=T;
plot(time{j}(1:length(time{j})),T);
grid on;
saveas(h3,[pwd '/folder_images/' namesb{2} '.fig']);
SNRratio=abs((mean(T)/std(T))/(mean(D{j})/std(D{j})))
Dp=T;
p=1;
temp1=0; %% reset signal temps
temp2=0;
tl=0;
n=1;
%% RR- Distance detection based on calculated time
for(k=2:length(Dp)-1)
    if(Dp(k)>=Dp(k-1) && Dp(k)>=Dp(k+1) && Dp(k)>0) %%IBI estimated
        G(p)=Dp(k); %% Signal R-peak detected
        tt(p)=time{j}(k);
        temp1=tt(p);
        tl=temp1-temp2; %% time distance calculation with temp1 and temp2
        if(p>=2)
            temp2=tt(p);
        end;
        p=p+1;
    end;
    if(tl>=0.1 && tl<=1.4) %% time amplitude filter to avoud another remanent noise
     ttp(n)=tl;%% assign the amplitude filter on a new variable ttp
     n=n+1;
    end;
    if (n==1)
          ttd(k)=ttp(n);
    else
          ttd(k)=ttp(n-1);
    end;
    timem(k)=time{j}(k);
end;
hold on
plot(tt,G,'r*');
xlabel('Time [s]');
ylabel('Signal Level');
h4=figure(4);
plot(timem,ttd,'g-');
hold on
%% Smoothing spline calculation
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = p_inter;
[fitresult, gof] = fit(downsample(timem,in)', downsample(ttd,in)', ft, opts ); %% downsample the D-IBI to calculate better the smooth and not have discontinuities 
Dinter=fitresult(1:length(downsample(ttd,in)));
plot(downsample(timem,in),Dinter','r-');
xlabel('Time [s]');
ylabel('RR interval [s]');
grid on
%% Writing CSV files
D1{1}=time{j};
D1{2}=D{j};
D1{3}=T;
D2{1}=downsample(timem,in);
D2{2}=Dinter';
writecvsfile('BVP_data.csv','IBI_data.csv',D1,D2);
saveas(h4,[pwd '/folder_images/' namesb{3} '.fig']);
h5=figure(5);
Fn{j}= (abs(fft((Dinter-mean(Dinter)),128*length(Dinter)))); %% Spectrum of th IIBI signal 
plot(linspace(0,fs(j)/(2*in),length(Fn{j})/2),Fn{j}(1:length(Fn{j})/2));
grid on;
saveas(h5,[pwd '/folder_images/' namesb{4} '.fig']);
h6=figure(6);
Qn{j}=cceps(xcorr(Dinter),128*length(Dinter)); %% Spectrum of th IIBI signal 
ylabel('Signal Level');
xlabel('Cepstral Sample');
plot(Qn{j});
grid on;
saveas(h6,[pwd '/folder_images/' namesb{5} '.fig']);
compareIBImnumeric(time{3},timem,D{3},ttd,Dinter,in,D{4},time{4});
save(['datasignals_' subt '_' session '.mat'],'D','Dinter','time','ttd','timem','T','G','tt');
A=1;

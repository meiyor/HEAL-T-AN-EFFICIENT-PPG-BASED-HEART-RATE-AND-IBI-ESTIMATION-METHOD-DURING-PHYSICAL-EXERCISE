function readfileHEALver2(subt,session,type,tamt,W,seld)
%% subt:subject index that you want to explore inside /Data folder
%% session: session that you desire to characterize per subject.
%% type: is a cell that represent the prefix of each included per session.
%% fs: is a vector of the same size of type cell, that represent the sampling frequency defined by the Emphatica handheld
%% in: is the factor of downsampling that should be taking into account to calculate the interpolated IBI.
%% p_inter: this is the smoothing factor to interpolate the estimated IBI (between 0 and 1).
%% tam: User defined size for the current simulation (source BVP signal, PPG).
%% W: the number of windows to select the less variable signal segments from PPGsignal.
close all;
addpath(genpath([pwd '/eeglab12_0_2_5b']));
P=0;
if (seld==0)
 tam=tamt*60*60*64; %% select hours time delimitatorl;
else
 tam=tamt;
end;
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
    temp=load([pwd '/Data/' subt '/' session '/' type{j} '.txt']);
    time{j}=temp(:,1)-start;
    if (strcmp(type{j},'ibi') || strcmp(type{j},'gsr') || strcmp(type{j},'bvp') || strcmp(type{j},'temp'))
        D{j}=temp(:,2);
    else
        D{j}=temp(:,2:4);  
    end;
    subplot(length(type),1,j);
    plot(time{j}(1:length(time{j})-1),D{j}(1:length(time{j})-1,:));
    grid on;
% temp=0;
% if(j>=2)
%     fclose(afile);
% end;
% i=1;
% S=kk;
% afile=fopen([pwd '/Data/' subt '/' session '/' type{j} '.txt'],'r');
% %% Only 1 channel data on session 
% if (strcmp(type{j},'ibi') || strcmp(type{j},'gsr') || strcmp(type{j},'bvp') || strcmp(type{j},'temp'))
%    while(~feof(afile) && i<=tam && str2num(S(1:17))-ennd<=0)
%        S=fscanf(afile,'%s',1);
%        time{j}(i)=str2num(S(1:17))-start;
%        D{j}(i)=str2num(S(19:length(S)));
%        temp=time{j}(i);
%        i=i+1;
%    end;
% elseif (strcmp(type{j},'acc')) %% if it is accel divided it in 3 file positions
%     while(~feof(afile) && i<=tam/2 && str2num(S(1:17))-ennd<=0)
%         S=fscanf(afile,'%s',1);
%         time{j}(i)=str2num(S(1:17))-start; %% all times relative with session starting time
%         temk=S(18:length(S));
%         cmm=find(temk==',');
%         D{j}(i,1)=str2num(temk(cmm(1)+1:cmm(2)-1)); % X axis
%         D{j}(i,2)=str2num(temk(cmm(2)+1:cmm(3)-1)); % Y axis
%         D{j}(i,3)=str2num(temk(cmm(3)+1:length(temk))); % Z axis
%         i=i+1;
%     end;
% end;
% subplot(length(type),1,j);
% if (type{j}=='ibi')
%     yy=interp1(time{j},D{j},linspace(min(time{j}),max(time{j}),1000));
%      plot(linspace(min(time{j}),max(time{j}),1000),yy,'r'); %% plot each data accordingly with data asked on input parameters
%      hold on 
%      plot(time{j}(1:length(time{j})-1),D{j}(1:length(time{j})-1));
%      if (j==1)
%         P=[min(time{j}),max(time{j})]; %% when you want to synch with empatica ibi select ibi as first ibi
%      end;
%  else
%   if (length(P)~=2) 
%       if (strcmp(type{j},'acc')) 
%         plot(time{j}(1:length(time{j})-1),D{j}(1:length(time{j})-1,:));  
%       else
%         plot(time{j}(1:length(time{j})-1),D{j}(1:length(time{j})-1));
%       end;
%    else
%         pos1=max(find(time{j}<=P(1))); %% plot from synchronized time with empatica IBI for any signal
%         pos2=min(find(time{j}>=P(2)));
%         plot(time{j}(pos1:pos2),D{j}(pos1:pos2));
%    end;
% end;
% grid on;
% xlabel('Time [s]');temp=0;
% if(j>=2)
%     fclose(afile);
% end;
% i=1;
% S=kk;
% afile=fopen([pwd '/Data/' subt '/' session '/' type{j} '.txt'],'r');
% %% Only 1 channel data on session 
% if (strcmp(type{j},'ibi') || strcmp(type{j},'gsr') || strcmp(type{j},'bvp') || strcmp(type{j},'temp'))
%    while(~feof(afile) && i<=tam && str2num(S(1:17))-ennd<=0)
%        S=fscanf(afile,'%s',1);
%        time{j}(i)=str2num(S(1:17))-start;
%        D{j}(i)=str2num(S(19:length(S)));
%        temp=time{j}(i);
%        i=i+1;
%    end;
% elseif (strcmp(type{j},'acc')) %% if it is accel divided it in 3 file positions
%     while(~feof(afile) && i<=tam/2 && str2num(S(1:17))-ennd<=0)
%         S=fscanf(afile,'%s',1);
%         time{j}(i)=str2num(S(1:17))-start; %% all times relative with session starting time
%         temk=S(18:length(S));
%         cmm=find(temk==',');
%         D{j}(i,1)=str2num(temk(cmm(1)+1:cmm(2)-1)); % X axis
%         D{j}(i,2)=str2num(temk(cmm(2)+1:cmm(3)-1)); % Y axis
%         D{j}(i,3)=str2num(temk(cmm(3)+1:length(temk))); % Z axis
%         i=i+1;
%     end;
% end;
% subplot(length(type),1,j);
% if (type{j}=='ibi')
%     yy=interp1(time{j},D{j},linspace(min(time{j}),max(time{j}),1000));
%      plot(linspace(min(time{j}),max(time{j}),1000),yy,'r'); %% plot each data accordingly with data asked on input parameters
%      hold on 
%      plot(time{j}(1:length(time{j})-1),D{j}(1:length(time{j})-1));
%      if (j==1)
%         P=[min(time{j}),max(time{j})]; %% when you want to synch with empatica ibi select ibi as first ibi
%      end;
%  else
%   if (length(P)~=2) 
%       if (strcmp(type{j},'acc')) 
%         plot(time{j}(1:length(time{j})-1),D{j}(1:length(time{j})-1,:));  
%       else
%         plot(time{j}(1:length(time{j})-1),D{j}(1:length(time{j})-1));
%       end;
%    else
%         pos1=max(find(time{j}<=P(1))); %% plot from synchronized time with empatica IBI for any signal
%         pos2=min(find(time{j}>=P(2)));
%         plot(time{j}(pos1:pos2),D{j}(pos1:pos2));
%    end;
% end;
% grid on;
% xlabel('Time [s]');
% saveas(h1,[pwd '/folder_images/timesignals.fig']);
% 
% saveas(h1,[pwd '/folder_images/timesignals.fig']);
end;
D{2}=sqrt(D{2}(:,1).^2+D{2}(:,2).^2+D{2}(:,3).^2)';
%[WW,AA]=windetection(D{1},D{2},W);
[WW,AA]=windetection(D{1},D{2},W);
[IBI1,IBI2,T1,T2,BVPn]=evalwinHRver2(WW,AA,D{3},time,W,0.6014,0,64,0,128);
f1=fopen([pwd '/outdir/' [subt '_' session] 'ibi.txt'],'w');
f2=fopen([pwd '/outdir/' [subt '_' session] 'iibi.txt'],'w');
f3=fopen([pwd '/outdir/' [subt '_' session] 'bvpn.txt'],'w');
for (k=1:length(T1))
    fprintf(f1,'%f,%f\n',[T1(k) IBI1(k)]);
end;
fclose(f1);
for (k=1:length(T2))
    fprintf(f2,'%f,%f\n',[T2(k) IBI2(k)]);
end;
fclose(f2);
for (k=1:length(T1))
    fprintf(f3,'%f,%f\n',[T1(k) BVPn(k)]);
end;
fclose(f3);

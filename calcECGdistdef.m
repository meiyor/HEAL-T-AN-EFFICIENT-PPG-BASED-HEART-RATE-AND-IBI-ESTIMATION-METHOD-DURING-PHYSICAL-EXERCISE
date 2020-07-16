function [HRVcont,Tdata]=calcECGdistdef(ECGin,fss)
%% to calculate this separately
ECG=cwt(ECGin,3,'db2'); %% 2 levels and db3 for subject #12
[ECGval,ECGpos]=findpeaks(ECG,'MINPEAKHEIGHT',50,'MINPEAKDISTANCE',30); 
Tdata=linspace(0,(1/fss)*length(ECGin),length(ECGin));
%%HRV ECG
ttd(1)=0;
tempn=1;
n=1;
for(p=1:length(ECGpos))
  ttl(p)=Tdata(ECGpos(p));
end;
for(l=1:length(ECG)-1)
   if (l==ECGpos(n) && n<length(ttl))
       ttd(n)=ttl(n+1)-ttl(n);
       tempn=ttd(n);
       n=n+1;
   end;
   ttr(l)=tempn;
 %ttr=1;
end;
HRVcont=60./ttr;
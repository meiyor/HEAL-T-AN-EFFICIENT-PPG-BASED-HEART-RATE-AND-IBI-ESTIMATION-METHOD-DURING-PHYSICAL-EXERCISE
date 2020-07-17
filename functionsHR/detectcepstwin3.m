function [W,Acnew,timeW,sW,p,sal]=detectcepstwin3(BVP,Acc,fss,sec_base,overlap_sec)
timebase_sec=sec_base;
time_overlap_sec=overlap_sec;
time_max=length(BVP)*(1/fss);
samplesbase=round((length(BVP)/time_max).*timebase_sec);
samplesover=round((length(BVP)/time_max).*time_overlap_sec);
n=1;
p=1;
c=0;
av=1;
if (samplesbase<=length(BVP) && all(isnan(BVP)~=1))
%% run the ICA here before do the segments to seperate the entire trials ICs. IF you want to segment it plese run ICA for each window as your suitability
while (n<=length(BVP))
   if (p==1) 
     W{p}=BVP(1:samplesbase);   
     Acnew{p}=Acc(1:samplesbase);
   else
     W{p}=BVP(samplesbase*(p-1)-samplesover*(p-1)+1:samplesbase*(p)-samplesover*(p-1));
     Acnew{p}=Acc(samplesbase*(p-1)-samplesover*(p-1)+1:samplesbase*(p)-samplesover*(p-1));
   end
    %% Group each window on different exercise methodology %% define here the whatever ind  you will like to use. For our baseline approch use ind=1 %% if you want to
    %% use an external accelerometer classifier you can create your ind based on its decision and create your own ranges. 
    ind=1;
    if (av==ind)
        c=c+1;
    else
        sW(av)=c;
        c=1;
    end;
    timeW{ind,c}=p;
    av=ind;
    n=samplesbase*(p+1)-samplesover*(p);
    sW(p)=ind;
    p=p+1;
    sal=0;
end;
 else
   W=BVP;
   Acnew=Acc;
   timeW=0;
   sW=0;
   p=0;
   sal=1;
end;
A=1;

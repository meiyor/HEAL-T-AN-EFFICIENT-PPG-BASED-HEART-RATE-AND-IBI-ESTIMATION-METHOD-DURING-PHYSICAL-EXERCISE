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
while (n<=length(BVP))
   if (p==1) 
     W{p}=BVP(1:samplesbase);   
     Acnew{p}=Acc(1:samplesbase);
   else
     W{p}=BVP(samplesbase*(p-1)-samplesover*(p-1)+1:samplesbase*(p)-samplesover*(p-1));
     Acnew{p}=Acc(samplesbase*(p-1)-samplesover*(p-1)+1:samplesbase*(p)-samplesover*(p-1));
   end
    %% Group each window on different exercise methodology
    timecheck=(1/fss)*(samplesbase*(p)-samplesover*(p-1));
    Fm{p}=abs(fft(W{p},100000));
   fl{p}=linspace(0,16,length(Fm{p})/2);
   Fm{p}=Fm{p}(1:length(Fm{p})/2);
   Ym{p}=sort(Fm{p});
   posm(p)=max(find(Fm{p}==Ym{p}(length(Ym{p}))));
   refp(p)=fl{p}(posm(p));
   if (refp(p)>=0 && refp(p)<=0.75) %% rest
       ind=1;
   elseif (refp(p)>0.75 && refp(p)<=1.5) %% 6km first phase
      ind=2;
   elseif (refp(p)>1.5 && refp(p)<=3.5) %% 12km first phase
      ind=3;
   elseif (refp(p)>3.5 && refp(p)<=5.5) %% 6km second phase
      ind=4;
   elseif (refp(p)>5.5 && refp(p)<=7.5) %% 12km second phase
      ind=5;
    else % rest second stage
         ind=6;
    end;    
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

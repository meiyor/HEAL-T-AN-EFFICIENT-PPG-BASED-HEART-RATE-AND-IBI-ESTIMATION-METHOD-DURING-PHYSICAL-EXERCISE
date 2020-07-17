function [W,Acnew,timeW,sW,p]=detectcepstwin2(BVP,fss,sec_base,overlap_sec)
timebase_sec=sec_base;
time_overlap_sec=overlap_sec;
time_max=length(BVP(1,:))*(1/fss);
samplesbase=round((length(BVP(1,:))/time_max).*timebase_sec);
samplesover=round((length(BVP(1,:))/time_max).*time_overlap_sec);
n=1;
p=1;
c=0;
av=1;
while (n<=length(BVP))
    if (p==1)
      [As,Ws]=runica(BVP(:,1:samplesbase),'sphering','off');
      BVPn=icaproj(BVP(:,1:samplesbase),As*Ws,[1 2]);
    else
      [As,Ws]=runica(BVP(:,samplesbase*(p-1)-samplesover*(p-1)+1:samplesbase*(p)-samplesover*(p-1)),'sphering','off');
      BVPn=icaproj(BVP(:,samplesbase*(p-1)-samplesover*(p-1)+1:samplesbase*(p)-samplesover*(p-1)),As*Ws,[1 2]);
    end;
    W{p}=BVPn(:,1);
Acnew{p}=sqrt(BVPn(3,:).^2+BVPn(4,:).^2+BVPn(5,:).^2);
    %% Group each window on different exercise methodology
    timecheck=(1/fss)*(samplesbase*(p)-samplesover*(p-1));    
     if (timecheck>0 && timecheck<=30) %% rest
		       ind=1;
     elseif (timecheck>30 && timecheck<=90) %% 6km first phase
                       ind=2;
     elseif (timecheck>90 && timecheck<=150) %% 12km first phase
                       ind=3;
     elseif (timecheck>150 && timecheck<=210) %% 6km second phase
                       ind=4;
     elseif (timecheck>210 && timecheck<=270) %% 12km second phase
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
   p=p+1;
end;
sW(av)=c;
A=1;

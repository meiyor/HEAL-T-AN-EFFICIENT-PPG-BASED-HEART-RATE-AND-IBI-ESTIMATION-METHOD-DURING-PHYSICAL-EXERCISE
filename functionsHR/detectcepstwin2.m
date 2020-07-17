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
    %% run the ICA here before do the segments to seperate the entire trials ICs. IF you want to segment it plese run ICA for each window as your suitability
    if (p==1)
      [As,Ws]=runica(BVP(:,1:samplesbase),'sphering','off');
      BVPn=icaproj(BVP(:,1:samplesbase),As*Ws,[1 2]);
    else
      [As,Ws]=runica(BVP(:,samplesbase*(p-1)-samplesover*(p-1)+1:samplesbase*(p)-samplesover*(p-1)),'sphering','off');
      BVPn=icaproj(BVP(:,samplesbase*(p-1)-samplesover*(p-1)+1:samplesbase*(p)-samplesover*(p-1)),As*Ws,[1 2]);
    end;
    W{p}=BVPn(:,1);
    Acnew{p}=sqrt(BVPn(3,:).^2+BVPn(4,:).^2+BVPn(5,:).^2);
    %% Group each window on different exercise methodology %% define here the ind will you like to use for our baseline approch use ind=1 %% if you want to
    %% use an external accelerometer classifier you can create you ranges here and divide them 
         ind=1;
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

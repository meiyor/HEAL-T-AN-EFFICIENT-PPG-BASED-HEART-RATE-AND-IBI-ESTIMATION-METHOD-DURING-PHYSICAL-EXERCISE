function [MR,pmaxr,MRc,pi,index1,Dp,G,tt]=probereconst(X,Acc,fss,M,mm,pr,ss,ind,kin,sel,selk,W_n)
%% X: BVP signal current window index
%% Acc: Accelerometer current window  index
%% kk: optional - cepstral indexes to be extracted from the orignal signal
%% fss: BVP signal original sampling frequency
%% M: Moving average filter order
%% mm: number of samples to calculate the maximum BVP points on the final BVP filtered signal (R-peaks)
%% pr: updated FIR filter frequency range for dynamic variation
%% ss: plotting selector (e.g. BVP and R-peak detection)
%% ind: HR-peak detection process activition-indicator
%% kin: spectrum peaks selection index for HR-peak detection process.
%% sel,selk: method selectors
%% W_n: length of the window
close all;
[Xn,N]=cceps(X);
Y=icceps(Xn,N);
if (ss==1)
figure(1);
plot(linspace(0,8,length(X)),X);
grid on;
figure(2);
plot(linspace(0,8,length(Y)),Y,'r');
grid on;
end;
p=1;
pi=0;
%% Moving average filter
Yf=conv(Y,(1/M).*ones([1 M]));
Yf=Yf(round(M/2):length(Y)+round(M/2)-1);
Yf=filtsig(Yf,0,fss,700,pr(1),pr(2)); 
Yf=Yf(length(Yf)/2-length(Y)/2+1:length(Yf)/2+length(Y)/2);
Af=filtsig(Acc,0,fss,700,pr(1),pr(2)); 
Af=Af(length(Yf)/2-length(Acc)/2+1:length(Yf)/2+length(Acc)/2);
An=(abs(fft((xcorr(Yf)),500000)));
Ap=(abs(fft((xcorr(Af)),500000)));
pmaxr=1;
if (sel==0 || sel==1) 
    if (selk==1)
        Qex=createxpmat(length(Yf),3*length(Yf));
        Qexa=createxpmat(length(Af),3*length(Yf));
        An=MFOCUSS(Qex,Yf'./max(Yf),0.1,'p',0.1,'prune_gamma',1e-4,'max_iters',1000,'epsilon',1e-8,'print',0);
        Ap=MFOCUSS(Qexa,Af'./max(Af),0.1,'p',0.1,'prune_gamma',1e-4,'max_iters',1000,'epsilon',1e-8,'print',0);
        An=abs(An);
        Ap=abs(Ap);
    end;
    [pq,pmaxr,pi,index1]=indexlookspect(An./max(An),Ap./max(Ap),linspace(0,fss/2,length(An)/2),[pr(1) pr(2)],ind,kin,0);
end;
Yf=filtsig(Yf,0,fss,700,pq(1),pq(2)); 
Yf=Yf(length(Yf)/2-length(Y)/2+1:length(Yf)/2+length(Y)/2);
Dp=Yf./max(Yf);
temp2=0;
tl=0;
ttp(1)=0;
ttd(1)=0;
n=1;
[par1,pcomplete]=findpeaks(Yf,'MINPEAKHEIGHT',0.2*max(Yf),'MINPEAKDISTANCE',30);
tcomplete=pcomplete.*(1/fss);
for(k=2:length(Dp)-1)
    if(Dp(k)>=Dp(k-1) && Dp(k)>=Dp(k+1) && Dp(k)>0 && std([Dp(k-1) Dp(k) Dp(k+1)])>=0.00001) %%IBI estimated
      f1=0; %ack for incremental slope
      f2=0; %ack for decremental slope
       if(k>=mm && k<=length(Dp)-mm)
          %%slope variation
          b=0; %% base incremental slope
          for(l=k-(mm-1):k)
             if (Dp(l)<b)
                 f1=1;
                 break;
             end;
             b=Dp(l);
          end;
          c=Dp(k); %% base decremental slope
          for(l=k:k+(mm-1))
             if (Dp(l)>c)
                 f2=1;
                 break;
             end;
             c=Dp(l);
          end;  
      end;
      if (f1==0 && f2==0);
        G(p)=Dp(k); %% Signal R-peak detected
        tt(p)=(k-1)/fss;
        %% time corrected
        pos1=min(find(tcomplete>=tt(p)));
        pos2=max(find(tcomplete<=tt(p)));
        if(abs(tcomplete(pos1)-tt(p))>=abs(tcomplete(pos2)-tt(p)))
            tt(p)=tt(p)-abs(tcomplete(pos2)-tt(p));
        else
           if (~isempty(pos1))
            tt(p)=tt(p)+abs(tcomplete(pos1)-tt(p));
           else
             tt(p)=tt(p)-abs(tcomplete(pos2)-tt(p));
           end;
        end;
        temp1=tt(p);
        tl=temp1-temp2; %% time distance calculation with temp1 and temp2
        if(p>=2)
            temp2=tt(p);
        end;
        p=p+1;
      end;
    end;
    if(tl>=0.22 && tl<=1.7) %% time amplitude filter to avoud another remanent noise
     ttp(n)=tl;%% assign the amplitude filter on a new variable ttp
     n=n+1;
    end;
    if (n==1)
          ttd(k)=ttp(n);
    else
          ttd(k)=ttp(n-1);
    end;
    timem(k)=(k-1)./fss;
end;
if (ss==1)
    figure(3)
    plot(linspace(0,W_n-1,length(Dp)),Dp);
    hold on
    t_po=linspace(0,W_n-1,length(Dp));
    plot(t_po(tt*fss),G,'r*');
    grid on;
    figure(4)
    posm=max(find(ttd==0));
    if(posm==999)
       ttd(posm+1)=0.5; 
    end;
    plot(timem(posm+1:length(ttd)),60./ttd(posm+1:length(ttd)));
    grid on;
end;
posm=max(find(ttd==0));
if(posm==999)
   ttd(posm+1)=0.5; 
end;
MR=pmaxr*60; %% calculated by frequency
MRc=[ttd(1) ttd(2:length(ttd))]; %% calculated by time

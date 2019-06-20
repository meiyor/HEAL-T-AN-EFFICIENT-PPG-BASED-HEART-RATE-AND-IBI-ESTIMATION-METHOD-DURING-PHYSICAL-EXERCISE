function [MR,pmaxr,Dp,ppn]=recon_filt(BVP_values,Acc,fss,M,mm,pr,ss,ind,kin,sel)
close all;
%% cepstral reconstruction first (it doesn't change the signal representation)
[BVP_valuesn,N]=cceps(BVP_values);
v=ones([1 round(length(BVP_valuesn)/2)]);
Y=icceps(BVP_valuesn,N);
if (ss==1) %% only if you want to plot PPG signal in the middle.
    figure(1);
    plot(linspace(0,8,length(BVP_values)),BVP_values);
    grid on;
    figure(2);
    plot(linspace(0,8,length(Y)),Y,'r');
    grid on;
end;
p=1;
%% Moving average filter
Yf=conv(Y,(1/M).*ones([1 M]));
Yf=Yf(round(M/2):length(Y)+round(M/2)-1);
pmaxr=1;
if (sel==0) %% application of the moving average filter
    Yf=filtsig(Yf,0,fss,500,pr(1),pr(2)); 
    Yf=Yf(length(Yf)/2-length(Y)/2+1:length(Yf)/2+length(Y)/2);
    Af=filtsig(Acc,0,fss,500,pr(1),pr(2)); 
    Af=Af(length(Af)/2-length(Acc)/2+1:length(Af)/2+length(Acc)/2); 
    if ((iscolumn(Yf) && iscolumn(Af)) || (isrow(Yf) && isrow(Af)))
        if (length(Yf)<=length(Af))
             An=(abs(fft((xcorr(Yf)-xcorr(Af(1:length(Yf)))),10000)));                                               
        else                                                                                                      
             An=(abs(fft((xcorr(Yf(1:length(Af)))-xcorr(Af)),10000)));
        end;
    else
        if (length(Yf)<=length(Af))
            An=(abs(fft((xcorr(Yf')-xcorr(Af(1:length(Yf)))),10000)));
        else 
             An=(abs(fft((xcorr(Yf(1:length(Af))')-xcorr(Af)),10000)));  
        end;
    end;
    [pq,pmaxr,ppn]=indexlookspect(hamming(length(An))'.*(abs(hilbert(An./max(An)))),linspace(0,fss/2,length(An)/2),[pr(1) pr(2)],ind,kin,ss);
    Yf=filtsig(Yf,0,fss,500,pq(1),pq(2)); 
    Yf=Yf(length(Yf)/2-length(Y)/2+1:length(Yf)/2+length(Y)/2);
end;
Dp=Yf./max(Yf);
temp1=0; %% reset signal temps
temp2=0;
tl=0;
ttp(1)=0;
ttd(1)=0;
n=1;
[par1,pcomplete]=findpeaks(Yf,'MINPEAKHEIGHT',20,'MINPEAKDISTANCE',10); %% 40,30 for 10 window
tcomplete=pcomplete.*(1/fss);
for(k=2:length(Dp)-1)
    if(Dp(k)>=Dp(k-1) && Dp(k)>=Dp(k+1) && Dp(k)>0 && std([Dp(k-1) Dp(k) Dp(k+1)])>=0.0001) %%IBI estimated
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
%         %% time corrected
%         pos1=min(find(tcomplete>=tt(p)));
%         pos2=max(find(tcomplete<=tt(p)));
%         if(abs(tcomplete(pos1)-tt(p))>=abs(tcomplete(pos2)-tt(p)))
%             tt(p)=tt(p)-abs(tcomplete(pos2)-tt(p));
%         else
%            if (~isempty(pos1))
%             tt(p)=tt(p)+abs(tcomplete(pos1)-tt(p));
%            else
%              tt(p)=tt(p)-abs(tcomplete(pos2)-tt(p));
%            end;
%         end;
        temp1=tt(p);
        tl=temp1-temp2; %% time distance calculation with temp1 and temp2
        if(p>=1)
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
plot(linspace(0,8,length(Dp)),Dp);
hold on
plot(tt,G,'r*');
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
MR=[ttd(1) ttd(1:length(ttd))];

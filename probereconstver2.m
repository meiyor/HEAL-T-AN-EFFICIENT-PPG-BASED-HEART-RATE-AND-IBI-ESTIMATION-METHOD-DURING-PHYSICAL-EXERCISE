function [HR_val,pmaxr,Dp]=probereconstver2(BVP,Acc,kk,fss,pp,M,mm,pr,plot_ind,ind,kin,sel)
close all;
[BVPn,N]=cceps(BVP);
v=ones([1 round(length(BVPn)/2)]);
[Qn Nn]=cceps(Acc);
ACCQ=sort(Qn(kk:round(length(Qn)/2)));
d=[kk(1:pp)];   
BVPn(d)=zeros([1 length(d)]);
Y=icceps(BVPn,N);
if (plot_ind==1)
figure(1);
plot(linspace(0,8,length(BVP)),BVP);
grid on;
figure(2);
plot(linspace(0,8,length(Y)),Y,'r');
grid on;
end;
p=1;
%% Moving average filter
Sal_Filter=conv(Y,(1/M).*ones([1 M]));
Sal_Filter=Sal_Filter(round(M/2):length(Y)+round(M/2)-1);
pmaxr=1;
if (sel==0) 
    Sal_Filter=filtsig(Sal_Filter,0,fss,500,pr(1),pr(2)); 
    Sal_Filter=Sal_Filter(length(Sal_Filter)/2-length(Y)/2+1:length(Sal_Filter)/2+length(Y)/2)';
    Af=filtsig(Acc,0,fss,500,pr(1),pr(2)); 
    Af=Af(length(Af)/2-length(Acc)/2+1:length(Af)/2+length(Acc)/2);
    An=(abs(fft((xcorr(Sal_Filter(1:length(Af)))-xcorr(Af)),10000)));
    [pq,pmaxr]=indexlookspect(hamming(length(An))'.*(abs(hilbert(An./max(An)))),linspace(0,fss/2,length(An)/2),[pr(1) pr(2)],ind,kin,plot_ind);
    Sal_Filter=filtsig(Sal_Filter,0,fss,500,pq(1),pq(2)); 
    Sal_Filter=Sal_Filter(length(Sal_Filter)/2-length(Y)/2+1:length(Sal_Filter)/2+length(Y)/2);
end;
Dp=Sal_Filter./max(Sal_Filter);
temp1=0; %% reset signal temps
temp2=0;
tl=0;
ttp(1)=0;
ttd(1)=0;
n=1;
[par1,pcomplete]=findpeaks(Sal_Filter,'MINPEAKHEIGHT',40,'MINPEAKDISTANCE',20);
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
    if(tl>=0.12 && tl<=2.1) %% time amplitude filter to avoud another remanent noise
     ttp(n)=tl;%% aplot_indign the amplitude filter on a new variable ttp
     n=n+1;
    end;
    if (n==1)
          ttd(k)=ttp(n);
    else
          ttd(k)=ttp(n-1);
    end;
    timem(k)=(k-1)./fss;
end;
if (plot_ind==1)
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
HR_val=[ttd(1) ttd(1:length(ttd))];
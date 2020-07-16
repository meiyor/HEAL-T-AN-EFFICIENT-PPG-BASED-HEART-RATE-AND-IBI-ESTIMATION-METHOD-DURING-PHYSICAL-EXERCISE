function [MR,llg,sBVP,TT,smoothpp]=evalwinHR3(W,AA,P,Q,ini_filt,inc_filt,ssW,selg,pinter,ssel,fss,selff,selp,W_n,overlap)
addpath([pwd '/BlandAltman']);
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = pinter;
temp=[];
for(i=1:P)
  %% finding zeros 
  pos=find(W{i}==0);
  if(length(pos)~=0)
      for(k=1:length(pos))
          W{i}(pos(k))=0.000001;
      end;    
  end;
  profile on;
  if (ssel==1) %% baseline implementation
     PPGf=filtsig(W{i},0,fss,500,0.7,4.5); % initial filter
     %PPGf=PPGf(length(PPGf)/2-length(W{i})/2+1:length(PPGf)/2+length(W{i})/2);
     Hs=dsp.SignalSource(PPGf,'SamplesPerFrame',1,'SignalEndAction','Cyclic repetition');
     Hn=dsp.SignalSource(AA{i},'SamplesPerFrame',1,'SignalEndAction','Cyclic repetition');
     Hd = dsp.FIRFilter('Numerator',fir1(31,[0.0412 0.321]));
     Hadapt = dsp.LMSFilter('Length',32,'StepSize', 0.0001);
    for k = 1:150
     s = step(Hs);
     d = step(Hd,s);   
     n = step(Hn); % Noise Accel
     [y,e]  = step(Hadapt,n',d');
     d=d-y';
    end   
    W{i}=d;
    posmat=1;
    dacc=1;
  else
  n=0;
  PQ1=ini_filt{1}; 
  PQ2=ini_filt{2}; 
  if (ssW(i)>=1 && ssW(i)<=2)
     [MR(i),pqmax(i),MRc{i},index1,sBVP{i},TT{i}]=probereconst(W{i},AA{i},fss,20,10,PQ1,selg,0,0,ssel,selp,W_n)
     PQQ=PQ1;
     per=0.10;
  elseif (ssW(i)>=3 && ssW(i)<=4)
     [MR(i),pqmax(i),MRc{i},index1,sBVP{i},TT{i}]=probereconst(W{i},AA{i},fss,20,10,PQ1,selg,0,0,ssel,selp,W_n)
     PQQ=PQ1;   
     per=0.10;
  else
     [MR(i),pqmax(i),MRc{i},index1,sBVP{i},TT{i}]=probereconst(W{i},AA{i},fss,20,10,PQ2,selg,0,0,ssel,selp,W_n)
     PQQ=PQ2;
     per=0.10;
  end;
  inv_p=0;
  %% We have to use 2 previous indexes for a continuous HR regression (based on ECG peaks) and 4 previous peak per each 8s windows approach  
  if (i>=5)
      pi=0;
      while((pqmax(i)>=mean(pqmax(i-4:i-1))+per*mean(pqmax(i-4:i-1))))
          if (n==1)      
               PQQ(2)=pqmax(i)-0.001;
               PQQ(1)=PQQ(1)+0.3;
           if (length(index1)<=10)     
                      per=0.10;
              else
                    per=0.05;
              end;
          end;
          [MR(i),pqmax(i),MRc{i},pi,index1]=probereconst(W{i},AA{i},fss,20,10,PQQ,selg,1,n,ssel,selp,W_n)
          if (pi==1)
               break;
           end
          n=n+1;
          inv_p=1;
      end;
      %% std(pqmax(i-4:i-1))
       while((pqmax(i)<=mean(pqmax(i-4:i-1))-per*mean(pqmax(i-4:i-1))))
           if (n==1)   
              PQQ(1)=pqmax(i)+0.001;
              PQQ(2)=PQQ(2)-0.3;
           if(length(index1)<=10)     
                      per=0.10;
              else
                    per=0.05;
              end;
           end;
           [MR(i),pqmax(i),MRc{i},pi,index1]=probereconst(W{i},AA{i},fss,20,10,PQQ,selg,2,n,ssel,selp,W_n)
           if (pi==1)
               break;
           end
           n=n+1;
          inv_p=2;
       end;
  end;
  if(i>=2)
      posz=find(MRc{i}==0);
      MRc{i}(posz)=MRc{i-1}(length(MRc{i-1})).*ones([1 length(posz)]);
  end;
  if (selff==0) %% no window segmentatio continuous signal reconstruction
      temp=[temp MRc{i}];
  else
      if (i>=2)
        temp=[temp  (MRc{i-1}((W_n-overlap)*fss:W_n*fss-1)+MRc{i}(1:overlap*fss))/2]; %% for IBI overlap
        temp=[temp MRc{i}(overlap*fss:(W_n-overlap)*fss-1)]; 
      else
        temp=MRc{i}(1:(W_n-overlap)*fss-1);
      end;
      A=1;
  end;
  Q.BPM0(i);
  A=1;
 end;
end;
poszz=max(find(temp==0));
llg=length(temp(poszz+1:length(temp)));
HRVcontesti=conv(60./temp(poszz+1:length(temp)),(1/700).*ones([1 700]));
[fitresult, gof]=fit([1:1:P]',MR',ft,opts);
smoothpp=fitresult(1:P);
MRn=conv(MR,(1/5)*ones([1 5]));
MRq=MRn(5:length(MR)+4);
timee=toc
figure(1)
plot(MR);
hold on
plot(MRq,'k');
plot(Q.BPM0,'r');
plot(fitresult(1:P),'g');
grid on;
xlabel('Windows Index [W]');
ylabel('Heart Rate [BPM]');
if (selff==0)
    figure(2)
    plot(abs(MR-Q.BPM0'));
    hold on;
    plot(abs(MR-Q.BPM0')./Q.BPM0','r');
    plot(abs(smoothpp-Q.BPM0),'g');
    plot(abs(smoothpp-Q.BPM0)./Q.BPM0,'k');
    grid on;
    Error1=(1/P).*sum(abs(MR-Q.BPM0'))
    Error2=(1/P).*sum(abs((MR-Q.BPM0'))./Q.BPM0')
    Error3=(1/P).*sum(abs(smoothpp-Q.BPM0))
    Error4=(1/P).*sum(abs((smoothpp-Q.BPM0))./Q.BPM0)
end;
A=1;
%figure(6) %% to plot and calculate Blandman plot per each subject
%BlandAltman(MR',Q.BPM0,{'Predicted HR','Ground-Truth (ECG) HR'});
%grid on;

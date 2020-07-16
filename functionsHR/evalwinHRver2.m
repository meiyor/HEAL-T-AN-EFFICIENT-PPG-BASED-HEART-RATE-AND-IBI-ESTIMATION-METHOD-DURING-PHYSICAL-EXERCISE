function [temp,smoothpp,T1,T2,tempbvp]=evalwinHRver2(W,AA,IBI,timem,P,pinter,ssel,fss,selg,F)
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
%opts.Normalize = 'on';
opts.SmoothingParam = pinter;
temp=[];
tempbvp=[];
%smoothpp=[];
for(i=1:P)
  %% finding zeros 
  pos=find(W{i}==0);
  if(length(pos)~=0)
      for(k=1:length(pos))
          W{i}(pos(k))=0.000001;
      end;    
  end;
  if (ssel==0)
  [Xn{i},N]=cceps(W{i});
  [Qn{i} Nn]=cceps(AA{i});
  Y{i}=abs(fft(xcorr(W{i})));
  G{i}=abs(fft(xcorr(AA{i})));
  figure(3);
  plot(linspace(0,fss/2,length(Y{i})/2),Y{i}(1:length(Y{i})/2)./max(Y{i}(1:length(Y{i})/2)));
  hold on
  plot(linspace(0,fss/2,length(G{i})/2),G{i}(1:length(G{i})/2)./max(G{i}(1:length(G{i})/2)),'r');
  grid on;
  %%Cepstrum Accel sort
  ACCs{i}=sort(Qn{i}(10:length(Qn{i})/2));
  for (k=1:100)
   dacc(k)=max(find(Qn{i}(10:length(Qn{i})/2)==ACCs{i}(k)))+9; 
  end;
  %%Cepstrum BVP sort
  BVPs{i}=sort(Xn{i}(10:length(Xn{i})/2));
  for (k=1:100)
   dbvp(k)=max(find(Xn{i}(10:length(Xn{i})/2)==BVPs{i}(k)))+9; 
  end;
  posmat=find(ismember(dacc,dbvp)==1);
    %%Cepstrum Accel sort
  ACCs1{i}=sort(Qn{i}(length(Qn{i})/2+100:length(Qn{i})-10));
  for (k=1:100)
   dacc1(k)=max(find(Qn{i}(length(Qn{i})/2+100:length(Qn{i})-10)==ACCs1{i}(k)))+length(Qn{i})/2+100; 
  end;
  %%Cepstrum BVP sort
  BVPs1{i}=sort(Xn{i}(length(Xn{i})/2+100:length(Xn{i})-10));
  for (k=1:100)
   dbvp1(k)=max(find(Xn{i}(length(Xn{i})/2+100:length(Xn{i})-10)==BVPs1{i}(k)))+length(Xn{i})/2+100; 
  end;
  posmat1=find(ismember(dacc1,dbvp1)==1);
  else
     PPGf=filtsig(W{i},0,fss,500,0.7,4.5); 
     PPGf=PPGf(length(PPGf)/2-length(W{i})/2+1:length(PPGf)/2+length(W{i})/2);
     Hs=dsp.SignalSource(PPGf,'SamplesPerFrame',1,'SignalEndAction','Cyclic repetition');
     Hn=dsp.SignalSource(AA{i},'SamplesPerFrame',1,'SignalEndAction','Cyclic repetition');
     Hd = dsp.FIRFilter('Numerator',fir1(31,[0.0412 0.321]));
     Hadapt = dsp.LMSFilter('Length',32,'StepSize', 0.0001);
    for k = 1:50
     s = step(Hs);
     d = step(Hd,s);   
     n = step(Hn); % Noise Accel
     [y,e]  = step(Hadapt,n',d');
     d=d-y';
    end   
    W{i}=d;
    posmat=1;
    posmat1=1;
    dacc=1;
    dacc1=1;
 end;
  n=1;
  PQ1=[0.7 3.5]; 
  PQ2=[0.7 4.5];
  if (var(AA{i})<=256)
    [MR{i},pqmax(i),BVPP{i}]=probereconstver2(downsample(W{i},2),AA{i},dacc(posmat),fss/2,0,20,5,PQ1,selg,0,0,ssel)
     PQQ=PQ1;
  else
    [MR{i},pqmax(i),BVPP{i}]=probereconstver2(downsample(W{i},2),AA{i},dacc(posmat),fss/2,0,20,5,PQ2,selg,0,0,ssel)      
    PQQ=PQ2;
  end;
  if (i>=2)
      while(pqmax(i)>=pqmax(i-1)+0.22*pqmax(i-1))
         [MR{i},pqmax(i),BVPP{i}]=probereconstver2(downsample(W{i},2),AA{i},dacc(posmat),fss/2,0,20,5,PQQ,selg,1,n,ssel)
         n=n+1;
      end;
      while(pqmax(i)<=pqmax(i-1)-0.22*pqmax(i-1))
         [MR{i},pqmax(i),BVPP{i}]=probereconstver2(downsample(W{i},2),AA{i},dacc(posmat),fss/2,0,20,5,PQQ,selg,2,n,ssel)
         n=n+1;
      end;
 end;
 if(i>=2)
      posz=find(MR{i}==0);
      MR{i}(posz)=MR{i-1}(length(MR{i-1})).*ones([1 length(posz)]);
  end;
  temp=[temp MR{i}];
  tempbvp=[tempbvp BVPP{i}];
  %[fitresult, gof]=fit(downsample(timem{1}(1:length(MR{i})),F/2)',downsample(MR{i}',F/2),ft,opts);
  %smoothk=fitresult(1:length(downsample(MR{i}',F/2)));
  %smoothpp=[smoothpp smoothk'];
end;
%[fitresult, gof]=fit(downsample(timem{1},F)',downsample(temp',F/2),ft,opts);
%smoothpp=fitresult(1:length(downsample(timem{1},F)));
ssq=length(temp);
smoothpp=conv(temp,(1/100).*ones([1 100]));
smoothpp=smoothpp(50:ssq+49);
plot(downsample(timem{1},2)',temp(1:length(downsample(timem{1},2))));
hold on;
post=max(find(timem{3}<=max(downsample(timem{1},2))));
plot(timem{3}(1:post),IBI(1:post),'r');
plot(downsample(timem{1},2)',smoothpp(1:length(downsample(timem{1},2))),'g');
%plot(downsample(timem{1}(1:length(timem{1})),F),smoothpp(1:length(smoothpp)),'g');
grid on;
for(k=1:post)
    TT(k)=max(find(downsample(timem{1},2)<=timem{3}(k)));
    TTn(k)=max(find(downsample(timem{1},F)<=timem{3}(k)))
    Error(k)=abs(IBI(k)-temp(TT(k)));
    Errorn(k)=abs(IBI(k)-smoothpp(TTn(k)));
end;
figure(2)
plot(timem{3}(1:post),Error);
hold on;
plot(timem{3}(1:post),Errorn,'g');
grid on;
T1=downsample(timem{1},2);
T2=downsample(timem{1},F);
A=1;

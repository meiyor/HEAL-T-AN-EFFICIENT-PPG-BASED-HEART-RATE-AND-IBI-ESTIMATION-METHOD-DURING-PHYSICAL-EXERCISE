function [ssp]=evalwinmeas(W,AA,P,Q,Wtime,ssW,selg,pinter,ssel,fss,selff,selp)
addpath([pwd '/BlandAltman']);
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
%opts.Normalize = 'on';
opts.SmoothingParam = pinter;
temp=[];
for(i=1:1)
  %% finding zeros 
  pos=find(W{i}==0);
  if(length(pos)~=0)
      for(k=1:length(pos))
          W{i}(pos(k))=0.000001;
      end;    
  end;
%  profile on;
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
   dacc1(k)=1; %max(find(Qn{i}(length(Qn{i})/2+100:length(Qn{i})-10)==ACCs1{i}(k)))+length(Qn{i})/2+100; 
  end;
  %%Cepstrum BVP sort
  BVPs1{i}=sort(Xn{i}(length(Xn{i})/2+100:length(Xn{i})-10));
  for (k=1:100)
   dbvp1(k)=1; %max(find(Xn{i}(length(Xn{i})/2+100:length(Xn{i})-10)==BVPs1{i}(k)))+length(Xn{i})/2+100; 
  end;
  posmat1=find(ismember(dacc1,dbvp1)==1);
  else
     PPGf=filtsig(W{i},0,fss,500,0.7,4.5); 
     PPGf=PPGf(length(PPGf)/2-length(W{i})/2+1:length(PPGf)/2+length(W{i})/2);
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
    posmat1=1;
    dacc=1;
    dacc1=1;
 end;
  n=0;
  PQ1=[0.9 2.5]; %% 1.5 for subject #12
  PQ2=[1.7 3.5]; %% 1.2 for subject #2, 1.8 for subject #12
  %if (i>2 && MR(i-1)>=160)
   %  PQ1(1)=PQ1(1)+1;
   %  PQ2(1)=PQ2(1)+1;
 % end;
  if (i>=Wtime{1,1} && i<=Wtime{1,ssW(1)})
     [MR(i),pqmax(i),MRc{i},index1]=probereconst(W{i},AA{i},dacc(posmat),fss,0,20,10,PQ1,selg,0,0,ssel,selp)
     PQQ=PQ1;
     per=0.05;
  elseif (i>=Wtime{2,1} && i<=Wtime{2,ssW(2)})
     [MR(i),pqmax(i),MRc{i},index1]=probereconst(W{i},AA{i},dacc(posmat),fss,0,20,10,PQ2,selg,0,0,ssel,selp)
     PQQ=PQ2;   
     per=0.05;
  else
     [MR(i),pqmax(i),MRc{i},index1]=probereconst(W{i},AA{i},dacc(posmat),fss,0,20,10,PQ2,selg,0,0,ssel,selp)
     PQQ=PQ2;
     per=0.05;
  end;
  inv_p=0;
  %% We have to use 2 previous indexes for a continuous HR regression (based on ECG peaks) and 4 previous peak per each 8s windows approach  
  if (i>=5)
      pi=0;
      %std(pqmax(i-4:i))
      while((pqmax(i)>=mean(pqmax(i-4:i-1))+per*mean(pqmax(i-4:i-1))))
          if (n==1)      
               PQQ(2)=pqmax(i)-0.001;
               PQQ(1)=PQQ(1)+0.001;
           if (length(index1)<=10)     
                      per=0.10;
              else
                    per=0.05;
              end;
          end;
          [MR(i),pqmax(i),MRc{i},pi,index1]=probereconst(W{i},AA{i},dacc(posmat),fss,0,20,10,PQQ,selg,1,n,ssel,selp)
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
              PQQ(2)=PQQ(2)-0.001;
           if(length(index1)<=10)     
                      per=0.10;
              else
                    per=0.05;
              end;
           end;
           [MR(i),pqmax(i),MRc{i},pi,index1]=probereconst(W{i},AA{i},dacc(posmat),fss,0,20,10,PQQ,selg,2,n,ssel,selp)
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
  if (selff==0)
      temp=[temp MRc{i}];
  else
      if (i>=2)
        temp=[temp MRc{i-1}(1:249) (MRc{i-1}(249:249+749)+MRc{i}(1:750))/2]; %% for IBI overlap
      end;
      A=1;
  end;
  Q.BPM0(i)
  profile viewer
  p = profile('info');
  %profsave(p,'profile_results');
  profile clear;
  profile off;
  Aqq=p.FunctionTable;
  ssp=0;
  for (l=1:length(Aqq))
      B(l)=Aqq(l).TotalTime;
  end;
  B=sort(B);
  ssp=B(length(Aqq))+B(length(Aqq)-1)+B(length(Aqq)-2);
  A=1;
end;

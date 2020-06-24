function [temp,smoothpp,T1,T2,tempbvp]=eval_winHR_no_over(W_bvp,AA_acc,P,timem,ini_filt,inc_filt,pinter,ssel,fss,selg,s_spline)
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = pinter;
temp=[];
tempbvp=[];

for(i=1:P) %% iterating along all the segments.
  %% finding zeros 
  pos=find(W_bvp{i}==0);
  if(length(pos)~=0)
      for(k=1:length(pos))
          W_bvp{i}(pos(k))=0.000001; %% To avoid having blank spots in the segments and discontinuities 
      end;    
  end;
  
  PPGf=filtsig(W_bvp{i},0,fss,500,0.7,4.5); %% initial filter 
  %PPGf=PPGf(length(PPGf)/2-length(W_bvp{i})/2+1:length(PPGf)/2+length(W_bvp{i})/2); %% only filter PPG
  if (length(PPGf)>=length(AA_acc{i}))
    PPGf=PPGf(1:length(AA_acc{i}));
  else
    AA_acc{i}=AA_acc{i}(1:length(PPGf));
  end;
  
  if (ssel==1)
  %% Baseline Implementation
  if (isrow(PPGf))
     Hs=dsp.SignalSource(PPGf,'SamplesPerFrame',1,'SignalEndAction','Cyclic repetition');
  else
     Hs=dsp.SignalSource(PPGf','SamplesPerFrame',1,'SignalEndAction','Cyclic repetition'); 
  end;
     Hn=dsp.SignalSource(AA_acc{i},'SamplesPerFrame',1,'SignalEndAction','Cyclic repetition');
     Hd = dsp.FIRFilter('Numerator',fir1(31,[0.0412 0.321])); %% Settings from Han et. al
     Hadapt = dsp.LMSFilter('Length',32,'StepSize', 0.000000021);
    for k = 1:50
     s = step(Hs);
     d = step(Hd,s);   
     n = step(Hn); % Noise Accel model based on Han et. al
     [y,e]  = step(Hadapt,n',d');
     d=d-y';
    end   
    W_bvp{i}=d;
    posmat=1;
    posmat1=1;
    dacc=1;
    dacc1=1;
 
  else
  %% implementation of HEAL-T core
  
  n=1;
  PQ1=ini_filt{1}; %% initial BHW bandwidths 
  if (length(ini_filt)>=2)
    PQ2=ini_filt{2};
  else
    PQ2=ini_filt{1};  
  end;
  if (ssel==0)
  if (var(AA_acc{i})<=256)
    [HR_values{i},pqmax(i),BVPn{i}]=recon_filt(W_bvp{i},AA_acc{i},fss,20,5,PQ2,selg,0,0,ssel)
     PQQ=PQ2;
  else
    [HR_values{i},pqmax(i),BVPn{i}]=recon_filt(W_bvp{i},AA_acc{i},fss,20,5,PQ1,selg,0,0,ssel)      
    PQQ=PQ1;
  end;
  if (i>=2)
      while(pqmax(i)>=pqmax(i-1)+0.15*pqmax(i-1))
	  [HR_values{i},pqmax(i),BVPn{i},pi]=recon_filt(W_bvp{i},AA_acc{i},fss,20,5,PQQ,selg,1,n,ssel)
	  PQQ(2)=pqmax(i)-inc_filt(2);
	  PQQ(1)=PQQ(1)+inc_filt(1);
	  if (pi==1)
		break;
           end
          n=n+1;
      end;
      while(pqmax(i)<=pqmax(i-1)-0.15*pqmax(i-1))
	  [HR_values{i},pqmax(i),BVPn{i},pi]=recon_filt(W_bvp{i},AA_acc{i},fss,20,5,PQQ,selg,2,n,ssel)
	  PQQ(1)=pqmax(i)+0.001;
	  if (pi==1)
	        break;
           end
          PQQ(2)=PQQ(2)-0.3;
         n=n+1;
      end;
   end;
   else
     if (var(AA_acc{i})<=256)
       [HR_values{i},pqmax(i),BVPn{i},pi]=recon_filt(W_bvp{i},AA_acc{i},fss,20,5,PQ2,selg,0,0,ssel)
	    PQQ=PQ2;
     else
       [HR_values{i},pqmax(i),BVPn{i},pi]=recon_filt(W_bvp{i},AA_acc{i},fss,20,5,PQ1,selg,0,0,ssel)
	   PQQ=PQ1;
	end;
		if (i>=2)
		    while(pqmax(i)>=pqmax(i-1)+0.15*pqmax(i-1))
		      [HR_values{i},pqmax(i),BVPn{i},pi]=recon_filt(W_bvp{i},AA_acc{i},fss,20,5,PQQ,selg,1,n,ssel)
		      PQQ(2)=pqmax(i)-0.001;
		      if (pi==1)
                break;	   
              end;
             PQQ(1)=PQQ(1)+0.3;               
             n=n+1;
		 end;
		 while(pqmax(i)<=pqmax(i-1)-0.15*pqmax(i-1))
		   [HR_values{i},pqmax(i),BVPn{i},pi]=recon_filt(W_bvp{i},AA_acc{i},fss,20,5,PQQ,selg,2,n,ssel)
		   PQQ(1)=pqmax(i)+0.001;
		    if (pi==1)
			   break;
                      end
                   PQQ(2)=PQQ(2)-0.3;
                   n=n+1;
		end;
	end;
 end;
 if(i>=2)
      posz=find(HR_values{i}==0);
      HR_values{i}(posz)=HR_values{i-1}(length(HR_values{i-1})).*ones([1 length(posz)]);
  end;
  temp=[temp HR_values{i}];
  if (iscolumn(BVPn{i}))
     tempbvp=[tempbvp BVPn{i}'];
  else
    tempbvp=[tempbvp BVPn{i}];
  end;
end;
end;
ssq=length(temp);
if (s_spline==0) %% cubic spline 
    [fitresult, gof]=fit(timem{1},temp,ft,opts);
    smoothpp=fitresult(1:length(downsample(timem{1},F)));
else
    smoothpp=conv(temp,(1/100).*ones([1 100])); %% using a moving average filter for smoothing
    smoothpp=smoothpp(50:ssq+49);
end;
Tk=timem{1};
if (length(temp)<=length(Tk))
     plot(Tk(1:length(temp)),temp);
else
    plot(Tk,temp(1:length(Tk))); 
end;
hold on;
Tq=timem{1}; 
if (length(smoothpp)<=length(Tq))
  plot(Tq(1:length(smoothpp)),smoothpp,'g');
else
  plot(Tq,smoothpp(1:length(Tq)),'g'); 
end;
grid on;
T1=Tk;
T2=Tq;

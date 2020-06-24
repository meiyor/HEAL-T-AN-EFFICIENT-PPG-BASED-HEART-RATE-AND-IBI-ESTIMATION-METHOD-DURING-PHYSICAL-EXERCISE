function [pp,pprt,pqq]=indexlookspect(BVP,FBVP,prange,ind,k,sel)
%% BVP: 
pos1=max(find(FBVP<=prange(1)));
pos2=min(find(FBVP>=prange(2)));
temp=sort(BVP(pos1:pos2));
pqq=0;
[bb index1]=findpeaks(BVP(pos1:pos2),'MINPEAKHEIGHT',0.00001,'MINPEAKDISTANCE',1);
NN=sort(BVP(pos1+index1));
bb=find(BVP(pos1+index1)==NN(length(NN)));
if (sel==1)
    figure(5);
    plot(FBVP(pos1:pos2),(BVP(pos1:pos2)));
    hold on
    plot(FBVP(pos1+index1),BVP(pos1+index1),'*r');
    %plot(FBVP(pos1+index2),BVP(pos1+index2),'*g');
    grid on;
end;
if (temp(length(temp))~=NN(length(NN)))
    tkind=NN(length(NN));
else
    tkind=temp(length(temp));
end;
if (ind==0)
    ppref=find(BVP(pos1:pos2)==tkind);
elseif(ind==1)
    if (bb-k>=1)
        ppref=find(BVP(pos1:pos2)==BVP(pos1+index1(bb-k)));
    else
      if (bb-k+1>=1)
        ppref=find(BVP(pos1:pos2)==BVP(pos1+index1(bb-k+1)));
      else
	 ppref=find(BVP(pos1:pos2)==tkind);
         pqq=1;   
      end; 
     end;
   else
    if (bb+k<=length(index1))
        ppref=find(BVP(pos1:pos2)==BVP(pos1+index1(bb+k)));
    else
      if (bb+k-1<=length(index1))
        ppref=find(BVP(pos1:pos2)==BVP(pos1+index1(bb+k-1)));
      else
	ppref=find(BVP(pos1:pos2)==tkind);
        pqq=1; 
      end;
      end;
end;
Fnew=FBVP(pos1:pos2);
pprt=Fnew(ppref);
pp(1)=Fnew(ppref)-0.2;
pp(2)=Fnew(ppref)+0.2;

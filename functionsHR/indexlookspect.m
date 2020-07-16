function [pp,pprt,pav,index1]=indexlookspect(BVP,QBVP,FBVP,prange,ind,k,sel)
%QBVP=abs(hilbert(BVP));
pos1=max(find(FBVP<=prange(1)));
pos2=min(find(FBVP>=prange(2)));
%pav=0;
temp=sort(BVP(pos1:pos2));
tempn=sort(QBVP(pos1:pos2));
RBVP=BVP-QBVP;
[bb index1]=findpeaks(BVP(pos1:pos2),'MINPEAKHEIGHT',0.00001,'MINPEAKDISTANCE',1);
[bq index2]=findpeaks(QBVP(pos1:pos2),'MINPEAKHEIGHT',0.00001,'MINPEAKDISTANCE',1);
[bt index3]=findpeaks(RBVP(pos1:pos2),'MINPEAKHEIGHT',0.00001,'MINPEAKDISTANCE',1);
NN=sort(BVP(pos1+index1-1));
bb=find(BVP(pos1+index1-1)==NN(length(NN)));
QQ=sort(QBVP(pos1+index2-1));
if (length(QQ)~=0)   
   qq=find(QBVP(pos1+index2-1)==QQ(length(QQ)));
 else
   qq=find(QBVP(pos1+index2-1)==1);    
end;
if (sel==1)
    figure(5);
    plot(FBVP(pos1:pos2),(BVP(pos1:pos2)));
    hold on
    plot(FBVP(pos1+index1-1),BVP(pos1+index1-1),'*r');
    %plot(FBVP(pos1+index2),BVP(pos1+index2),'*g');
    grid on;
    figure(6)
    plot(FBVP(pos1:pos2),(QBVP(pos1:pos2)));
    hold on
    plot(FBVP(pos1+index2-1),QBVP(pos1+index2-1),'*r');
    %plot(FBVP(pos1+index2),BVP(pos1+index2),'*g');
    grid on;
    figure(7)
    plot(FBVP(pos1:pos2),(RBVP(pos1:pos2)));
    hold on
    plot(FBVP(pos1+index3-1),RBVP(pos1+index3-1),'*r');
    %plot(FBVP(pos1+index2),BVP(pos1+index2),'*g');
    grid on;
end;
if (temp(length(temp))~=NN(length(NN))) %% BVP spectrum
    tkind=NN(length(NN));
else
    tkind=temp(length(temp));
end;
if (length(QQ)~=0 && tempn(length(tempn))~=QQ(length(QQ))) %% Accel spectrum
    tkinda=QQ(length(QQ));
else
    tkinda=tempn(length(tempn));
end;
pqq1=find(BVP(pos1:pos2)==tkind);
pqq2=find(QBVP(pos1:pos2)==tkinda);
Fnew=FBVP(pos1:pos2);
if (length(index1)-(k-1)>=1 && k<=length(index1))
if (ind==0)
    ppref=find(BVP(pos1:pos2)==tkind);
    pav=0;
elseif(ind==1)
    if (all(abs(Fnew(pqq1)-FBVP(pos1+index2))>=0.1))
         ppref=find(BVP(pos1:pos2)==tkind);
         pav=1;
    else    
        if (length(index1)-k>=1)
            ppref=find(BVP(pos1:pos2)==BVP(pos1-1+index1(length(index1)-k)));
        else
            ppref=find(BVP(pos1:pos2)==BVP(pos1-1+index1(length(index1)-(k-1))));
        end;
        pav=0;
    end;
else
   if (all(abs(Fnew(pqq1)-FBVP(pos1-1+index2))>=0.1))
         ppref=find(BVP(pos1-1:pos2)==tkind);
         pav=1;
    else    
         if (k+1<=length(index1))
             ppref=find(BVP(pos1:pos2)==BVP(pos1-1+index1(k+1)));
         else
             ppref=find(BVP(pos1:pos2)==BVP(pos1-1+index1(k-1)));
        end;
        pav=0;
   end;
 end;
else
    if (k>length(index1))
       ppref=find(BVP(pos1:pos2)==BVP(pos1-1+index1(k-1)));
    end;
    if (length(index1)-k<=1)
       ppref=find(BVP(pos1:pos2)==BVP(pos1-1+index1(length(index1)-(k-1)+1)));
    end;
    pav=1;
end;
pprt=Fnew(ppref);
pp(1)=Fnew(ppref)-0.05;
pp(2)=Fnew(ppref)+0.05;

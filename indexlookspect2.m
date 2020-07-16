function [pp,pprt,pav]=indexlookspect2(BVP,QBVP,FBVP,prange,ind,k,sel)
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
NN=sort(BVP(pos1+index1));
bb=find(BVP(pos1+index1)==NN(length(NN)));
QQ=sort(QBVP(pos1+index2));
qq=find(QBVP(pos1+index2)==QQ(length(QQ)));
if (sel==1)
    figure(5);
    plot(FBVP(pos1:pos2),(BVP(pos1:pos2)));
    hold on
    plot(FBVP(pos1+index1),BVP(pos1+index1),'*r');
    %plot(FBVP(pos1+index2),BVP(pos1+index2),'*g');
    grid on;
    figure(6)
    plot(FBVP(pos1:pos2),(QBVP(pos1:pos2)));
    hold on
    plot(FBVP(pos1+index2),QBVP(pos1+index2),'*r');
    %plot(FBVP(pos1+index2),BVP(pos1+index2),'*g');
    grid on;
    figure(7)
    plot(FBVP(pos1:pos2),(RBVP(pos1:pos2)));
    hold on
    plot(FBVP(pos1+index3),RBVP(pos1+index3),'*r');
    %plot(FBVP(pos1+index2),BVP(pos1+index2),'*g');
    grid on;
end;
if (temp(length(temp))~=NN(length(NN))) %% BVP spectrum
    tkind=NN(length(NN));
else
    tkind=temp(length(temp));
end;
if (tempn(length(tempn))~=QQ(length(QQ))) %% Accel spectrum
    tkinda=QQ(length(QQ));
else
    tkinda=tempn(length(tempn));
end;
pqq1=find(BVP(pos1:pos2)==tkind);
pqq2=find(QBVP(pos1:pos2)==tkinda);
Fnew=FBVP(pos1:pos2);
pprt=Fnew(ppref);
pp(1)=Fnew(ppref)-0.05;
pp(2)=Fnew(ppref)+0.05;

function  runsubjectnew_heal_t(inpath,outpath,subt,tses,type,fss,ini_filt,inc_filt,W_num,overlap)
cd([pwd '/functions_HR']);
 
%% with this cycles we can read the current folder file that doesn't have the desired file prefix

for (j=1:length(type))
   temp=load([inpath '/' subt '/' tses '/' type{j} '.txt']);
   time{j}=temp(:,1)-temp(1,1);
  if (strcmp(type{j},'ibi') || strcmp(type{j},'gsr') || strcmp(type{j},'bvp') || strcmp(type{j},'temp'))
      D{j}=temp(:,2);
     else
      D{j}=temp(:,2:4);
   end;
%% given bvp and accel D{} would have two positions BVP and accel
end;
D{1}=downsample(D{1},2);
[HR,llg,sBVP,TT,smpp]=IEEEsignaldatabasePPGAccel(D,2,fss,0,0,0,0,0,1,ini_filt,inc_filt,W_num,overlap);
[HRn,llgn,sBVPn,TTn,smppn]=IEEEsignaldatabasePPGAccel(D,2,fss,0,0,2,0,0,1,ini_filt,inc_filt,W_num,overlap); %% calculated from point by point approach
Gsvp=cell2mat(TT');
templ=Gsvp(1,1:overlap*fss-1);
for(k=2:length(sBVP))
    temp=(Gsvp(k-1,(W_num-overlap)*fss:W_num*fss-1)+Gsvp(k,1:overlap*fss))/2;
    templ=[templ temp];
    templ=[templ Gsvp(k,overlap*fss:(W_num-overlap)*fss-1)]; 
end;
templ=[templ Gsvp(length(sBVP),(W_num-overlap)*fss+1:W_num*fss-1)];
close all;
Tdata=linspace(0,(1/fss)*length(D{1}),length(D{1}));
Tdatap=linspace(0,(1/fss)*length(D{2}),length(D{2}));
Tdatan=linspace(0,(1/fss)*length(D{1}),length(HR));
T_dataind=[1:1:length(HR)];
Tdatak=linspace(0,(1/fss)*length(D{1}),length(HRn));
Tdataq=linspace(0,(1/fss)*length(templ),length(templ));
h1=figure(1);
subplot(311)
plot(Tdatap,D{2},'LineWidth',2);
grid on;
xlabel('Time [s]');
ylabel('Quantized Signal Level');
title([subt '/' tses '/']);
legend('X axis','Y axis','Z axis');
subplot(312);
plot(Tdata,D{1}./max(D{1}));
hold on;
plot(Tdataq,templ,'r','LineWidth',2);
legend('raw BVP','cleansed BVP');
grid on;
ylabel('Normalized Signal Level')
xlabel('Time [s]');
subplot(313);
plot(Tdatan,HR,'b','LineWidth',2);
hold on
plot(Tdatan,smpp,'g','LineWidth',2);
plot(Tdatak,smppn,'k','LineWidth',2);
xlabel('Time [s]');
ylabel('Heart-rate [BPM]');
legend('HR','Smoothed HR','HR IBI smoothed');
grid on
%saveas(h1,'graphdata1.fig')
if (exist([ outpath '/'])==0)            
  mkdir([ outpath '/']);                                                                
end;                                                                                                                                                                          
pid=feature('getpid'); 
saveas(h1,[ outpath '/' num2str(pid) '_salfig.fig'])
saveas(gcf,[ outpath '/' num2str(pid) '_salfig.eps'],'eps2c');
saveas(gcf,[ outpath '/' num2str(pid) '_salfig.jpg'],'jpg');
f1=fopen([ outpath '/' num2str(pid) '_ibipeakselHEAL.txt'],'w');                                                     
  for (k=1:length(HR))                                                                                                                 
     fprintf(f1,'%i,%f,%f\n',[T_dataind(k) HR(k) smpp(k)]);                                                                                
   end;                                                                                                                                   
  fclose(f1); 
f2=fopen([ outpath '/' num2str(pid) '_ibiIBIselHEAL.txt'],'w');                                                  
 for (k=1:length(HRn))
    fprintf(f2,'%f,%f,%f\n',[Tdatak(k) HRn(k) smppn(k)]);
 end;
fclose(f2);num2str(pid);
f1=fopen([ outpath '/' 'ibipeakselHEAL.txt'],'w');
for (k=1:length(HR))
  fprintf(f1,'%i,%f,%f\n',[T_dataind(k) HR(k) smpp(k)]);
end;
fclose(f1);
f2=fopen([ outpath '/' 'ibiIBIselHEAL.txt'],'w');
for (k=1:length(HRn))
  fprintf(f2,'%f,%f,%f\n',[Tdatak(k) HRn(k) smppn(k)]);
end;
fclose(f2);num2str(pid)
disp(['Process ' num2str(pid) ' has ended succesfully']);
cd('..');

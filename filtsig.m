function Y=filtsig(T,sel,fs,M,fp1,fp2)
F=log10(abs(fft(T)));
if (sel==1)
    figure(1);
    title('Raw Signal');
    subplot(211);
    semilogy(linspace(0,fs/2,length(F)/2),F(1:length(F)/2));
    xlabel('Frequency [Hz]');
    grid on;
    subplot(212)
    plot(linspace(0,(1/fs)*length(T),length(T)),T,'r');
    grid on;
    xlabel('Time [S]');
end;
%% Filter To let pass and reduce noise.
%% FIR Filter between 2-40Hz to cover the most critical EEG signals
f=linspace(0,pi,1000);
Ap=0;
As=-60;
xp1=(fp1*pi/(fs));
xp2=(fp2*pi/(fs));
Hd=rectpuls(f,2*abs(xp1-xp2));
pos=max(find(f<=xp1));
Hd=circshift(Hd',[pos 1]);
Hd=Hd';
pos1=min(find(Hd==1));
pos2=max(find(Hd==1));
Hd=cat(2,Hd(1:pos1)+10^(As/10),cat(2,Hd(pos1+1:pos2)*10^(Ap/10),Hd(pos2+1:length(Hd))+10^(As/10)));
h=real(ifft(Hd));
h=cat(2,fliplr(h(1:length(h)/2)),fliplr(h(length(h)/2:length(h))))./max(h).*blackman(length(h)+1)';
if (sel==1)
 figure(2);
 freqz(h,[1])
end;
if (sel==1 || sel==0)
Y=conv(T,h);
else
 Y=h;   
end;
%% Averaging Filter M=input
if (sel==1)
    figure(3)
    freqz((1/M).*ones([1 M]),1);
end;
%%Y=conv(Y,(1/M).*ones([1 M]));
if (sel==1)
    F=(abs(fft(Y)));
    figure(4);
    title('Filtered Signal');
    subplot(211);
    plot(linspace(0,fs/2,length(F)/2),F(1:length(F)/2));
    xlabel('Frequency [Hz]');
    grid on;
    subplot(212)
    plot(linspace(0,(1/fs)*length(Y),length(Y)),Y,'r');
    grid on;
    xlabel('Time [S]');
end;
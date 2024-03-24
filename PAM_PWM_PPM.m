%%%%%%%%%%%%%%%%%%%%%%%%
%MATLAB CODE FOR PULSE MODULATIONS
%Author hassan saker
%%
%signal definitions
A0=1;      %%Amplitude
fm=800;    %%frequency of signal
fs=10*fm;   %%clock frequency
fd=100*fs;    %%sampling frequency
t=0:1/fd:5/fm-1/fd;
signal=A0*cos(2*pi*t*fm);
figure('Name','Signal');
plot(t,signal);
title('information signal');
xlabel('Time');
ylabel('Amplitude');
clock=0.5*square(2*pi*fs*t,10)+0.5;
figure('Name','Clock');
plot(t,clock);
title('Clock');
xlabel('Time');
ylabel('Amplitude');
axis([0 5/fs 0  1.2]);
difc=[max(diff(clock),0) 0];
difc(1)=1;
figure ('Name','Dirac');
plot(difc);
title('Dirac ');
xlabel('time');
ylabel('difc');
%%
%PAM
sample=0;
sampl_hold=zeros(1,length(signal));
for i=1:length(signal)
    if difc(i)>0
        sample=signal(i);
    end
    sampl_hold(i)=sample;
end
figure ('Name','sample');
plot(t,signal,'r');
hold on;
plot(t,  sampl_hold);
title('sample and hold');
xlabel('time');
ylabel('Amplitude');
legend('signal','sample and hold signal');
PAM=sampl_hold.*clock;
figure ('Name','PAM');
plot(PAM);
title('PAM ');
xlabel('time');
ylabel('Amplitude');
%%
%fourier transform for PAM 
N=length(PAM);
f=[-fd/2:fd/N:(fd/2)-(fd/N)];
PAM_f=fftshift(fft(PAM));
figure ('Name','fourier PAM');
plot(f,abs(PAM_f./max(PAM_f))); %% normalized to 1 by dividing to maximum value
title('F PAM ');
xlabel('frequency');
ylabel('Amplitude');
%%
%demodulation of PAM and that done by using LPF
N=6;
fc=fs/(2*(fd/2));
[b,a]=butter(N,Wc);  
freqz(b,a);
%restoring the signal
y=filter(b,a,PAM);
figure ('Name','PAM');
plot(y);
title('PAM demodulation ');
xlabel('time');
ylabel('Amplitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%PWM MODULATION
A=1.2;
sw=A*sawtooth(2*pi*fs*t);
figure ('Name','sawtooth signal with inforrmation signal');
subplot(2,1,1)
plot(t,sw,t,signal);
title('sawtooth signal with inforrmation signal ');
xlabel('time');
ylabel('Amplitude');
axis([0 10/fs -A-0.5  A+0.5]);
% the PWM modulated signal;
PWM=zeros(1,length(signal));
for i=1:length(signal)
    if sw(i)<= signal(i)
        PWM(i)=1;
    end
    if sw(i)>signal(i)
        PWM(i)=0;
    end
end
subplot(2,1,2);
plot(t,PWM);
title('PWM');
xlabel('time');
ylabel('Amplitude');
axis([0 10/fs 0 1.2]);
%%
%fourier transform for PWM 
N=length(PWM);
f=[-fd/2:fd/N:(fd/2)-(fd/N)];
PWM_f=fftshift(fft(PWM));
figure ('Name','fourier PWM');
plot(f,abs(PWM_f./max(PWM_f))); %% normalized to 1 by dividing to maximum value
title('F PWM ');
xlabel('frequency');
ylabel('Amplitude');
%%
%Demodulation of PWM
Ts=1/fs;
Td=1/fd;
Tm=1/fm;
soep=Ts/Td; %% number of samples over each period
for k=1:floor((5*Tm-Td)/(Ts))
    sum=0;
    for j=(k-1)*(Ts/Td):k*(Ts/Td)
    sum=sum+PWM(j+1);
    end
    pd(k)=(sum-(soep/2))/(soep);
end
figure ('Name','detect PWM');
plot(pd);
title('Demod PWM ');
xlabel('time');
ylabel('Amplitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%PPM MODULATION
W=0.1*soep;  %% width of pulse 
N=length(PWM);
PPM=zeros(1,N);
i=1;
while i<=N-1
    if (PWM(i)==1 && PWM(i+1)==0)
    for j=i:W+i
        if(j<=N)
        PPM(j)=1;
        else
            break;
        end
    end
    i=W+i-1;
    end
    i=i+1;
end
plot(PWM,'LineWidth',1.5);
hold on;
plot(PPM,'r','LineWidth',1.5);
title('Modulation PPM ');
xlabel('time');
ylabel('Amplitude');
legend('PWM','PPM');
axis([0 1000 0 1.2]);
%%
%fourier transform for PPM 
N=length(PPM);
f=[-fd/2:fd/N:(fd/2)-(fd/N)];
PMM_f=fftshift(fft(PPM));
figure ('Name','fourier PPM');
plot(f,abs(PMM_f./max(PMM_f))); %% normalized to 1 by dividing to maximum value
title('F PPM ');
xlabel('frequency');
ylabel('Amplitude');
%%
%PWM MODULATION
W=0.1*soep;  %% width of pulse 
N=length(PPM);
det_PWM=zeros(1,N);
%detect the fallen edge
i=1;
k=1;
while i<=N-1
    if (PPM(i)==0 && PPM(i+1)==1)
    fallen_edge(k)=i+1;
    k=k+1;
    end
    i=i+1;
end
%forming PWM by assign 1 till I met the fallen edge
k=1;
v=1;
while k<=N
    if(k==fallen_edge(v))
             for j=k:soep*v
         det_PWM(j)=0;
            end
         k=v*soep;
         v=v+1;
         
    end
    det_PWM(k)=1;
    k=k+1;
end

plot(det_PWM,'LineWidth',1.5);
hold on;
plot(PPM,'r','LineWidth',1.5);
title('demod PPM ');
xlabel('time');
ylabel('Amplitude');
legend('det PWM','PPM');
axis([0 1000 0 1.2]);
%Demodulation of PWM
Ts=1/fs;
Td=1/fd;
Tm=1/fm;
soep=Ts/Td; %% number of samples over each period
for k=1:floor((5*Tm-Td)/(Ts))
    sum=0;
    for j=(k-1)*(Ts/Td):k*(Ts/Td)
    sum=sum+det_PWM(j+1);
    end
    pd(k)=(sum-(soep/2))/(soep);
end
figure ('Name','detect PWM');
plot(pd);
title('Demod PWM ');
xlabel('time');
ylabel('Amplitude');

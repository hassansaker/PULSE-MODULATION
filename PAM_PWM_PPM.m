%%%%%%%%%%%%%%%%%%%%%%%%
%MATLAB CODE FOR PULSE MODULATIONS
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

%% 
%PWM MODULATION
A=1.2;
sw=A*sawtooth(2*pi*fs*t+pi);
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
figure ('Name','fourier PAM');
plot(f,abs(PWM_f./max(PWM_f))); %% normalized to 1 by dividing to maximum value
title('F PAM ');
xlabel('frequency');
ylabel('Amplitude');
%%
%Demodulation of PWM


clear all

Fs=8000;
A=1;
F0=440;
phi=0;
len=Fs*0.008;

t=(0:len-1)/Fs;
f=(0:len-1)*(Fs/len);

x=A*sin(2*pi*F0*t+phi);
y=fft(x);

y

figure(1)
clf
plot(f, abs(y))
grid on
title('Amplitude spectrum')
xlabel('Frequency[Hz]')
ylabel('Amplitude')
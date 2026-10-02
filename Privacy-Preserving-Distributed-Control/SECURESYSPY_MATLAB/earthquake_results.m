clc;clear;
close all


%% elcentro
ax1=figure();
ax1.Color='white';
load ('earthquake_results.mat');
subplot(3,1,1);
plot(t_el,el_acc,'black','Linewidth',1.5);
xlim([0 35]);
ylabel('Base Acceleration (g)','Fontsize',14,'interpreter', 'Latex')
xlabel('Time (sec)','Fontsize',14,'interpreter', 'Latex')
grid on
subplot(3,1,2);
plot(el_freq,abs(el_four),'black','Linewidth',1.5);
ylabel('Fourier Amplitude','Fontsize',14,'interpreter', 'Latex')
xlabel('Frequency (Hz)','Fontsize',14,'interpreter', 'Latex')
xlim([0,max(abs(el_freq))])
grid on
%el_four_max=max(abs(el_four));
yf_el=20*log10(abs(el_four));
subplot(3,1,3);
plot(el_freq,yf_el,'black','Linewidth',1.5);
ylabel('Fourier Amplitude (dB)','Fontsize',14,'interpreter', 'Latex')
xlabel('Frequency (Hz)','Fontsize',14,'interpreter', 'Latex')
xlim([0,50])
grid on
%savefig('EL CENTRO.fig')
%% TABAS
ax2=figure;
ax2.Color='white';
subplot(3,1,1);
plot(t_TAB,TAB_acc,'black','Linewidth',1.5);
ylabel('Base Acceleration (g)','Fontsize',14,'interpreter', 'Latex')
xlabel('Time (sec)','Fontsize',14,'interpreter', 'Latex')
xlim([0 35]);
grid on
subplot(3,1,2);
plot(TAB_freq,abs(TAB_four),'black','Linewidth',1.5);
ylabel('Fourier Amplitude','Fontsize',14,'interpreter', 'Latex')
xlabel('Frequency (Hz)','Fontsize',14,'interpreter', 'Latex')
xlim([0,max(abs(TAB_freq))])
grid on
%TAB_four_max=max(abs(TAB_four));
yf_TAB=20*log10(abs(TAB_four));
subplot(3,1,3);
plot(TAB_freq,yf_TAB,'black','Linewidth',1.5);
ylabel('Fourier Amplitude (dB)','Fontsize',14,'interpreter', 'Latex')
xlabel('Frequency (Hz)','Fontsize',14,'interpreter', 'Latex')
xlim([0,25])
grid on
%savefig('TABAS.fig')
%% Loma Prieta
ax3=figure;
ax3.Color='white';
subplot(3,1,1);
plot(t_Loma,Loma_acc,'black','Linewidth',1.5);
ylabel('Base Acceleration (g)','Fontsize',14,'interpreter', 'Latex')
xlabel('Time (sec)','Fontsize',14,'interpreter', 'Latex')
xlim([0 35]);
grid on
subplot(3,1,2);
plot(Loma_freq,abs(Loma_four),'black','Linewidth',1.5);
ylabel('Fourier Amplitude','Fontsize',14,'interpreter', 'Latex')
xlabel('Frequency (Hz)','Fontsize',14,'interpreter', 'Latex')
xlim([0,max(abs(Loma_freq))])
grid on
%Loma_four_max=max(abs(Loma_four));
yf_Loma=20*log10(abs(Loma_four));
subplot(3,1,3);
plot(Loma_freq,yf_Loma,'black','Linewidth',1.5);
ylabel('Fourier Amplitude (dB)','Fontsize',14,'interpreter', 'Latex')
xlabel('Frequency (Hz)','Fontsize',14,'interpreter', 'Latex')
xlim([0 100])
grid on
%savefig('LOMA PRIETA.fig')
%% Frequency response/sampling rate or sampling frequency
load ('results.mat');
fs=1/0.008; %Hz %1/Ts
tlength=length(taim);
df=fs/tlength; %hertz per sample
f = -fs/2:df:fs/2-df + (df/2)*mod(tlength,2);
disp1_fft=fft(disp1,4375);
disp1c_fft=fft(disp1c,4375);
disp1cs_fft=fft(disp1cs,4375);
ax4=figure();
subplot(3,1,1);
plot(f,abs(disp1_fft),'black','Linewidth',1.5);
xlim([0 inf]);
grid on
subplot(3,1,2);
plot(f,abs(disp1c_fft),'black','Linewidth',1.5);
xlim([0 inf]);
grid on
subplot(3,1,3);
plot(f,abs(disp1cs_fft),'black','Linewidth',1.5);
xlim([0 inf]);
grid on
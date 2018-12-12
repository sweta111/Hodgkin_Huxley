%THIS PROGRAM DEMONSTRATES HODGKIN HUXLEY MODEL IN CURRENT CLAMP EXPERIMENTS AND SHOWS ACTION POTENTIAL PROPAGATION
%Time is in secs, voltage in mvs, conductances in m mho/mm^2, capacitance in uF/mm^2

% threshold value of current is 0.0223
clear all;
close all;
clc;
k=1;
istep=0.0001;
for ImpCur=0:istep:0.8 
%TimeTot=input('enter the time for which stimulus is applied in milliseconds');
gkmax=.36;
vk=-77;
gnamax=1.20;
vna=50;
gl=0.003;
vl=-54.387;
cm=.01; 

dt=0.01;
niter=100000;%change from 10,000 to 100,000
t=(1:niter)*dt;
iapp=ImpCur*ones(1,niter);
%for i=1:100
% iapp(1,i)=ImpCur;
%end;
v=-64.9964;
m=0.0530;
h=0.5960;
n=0.3177;

gnahist=zeros(1,niter);
gkhist=zeros(1,niter);
vhist=zeros(1,niter);
mhist=zeros(1,niter);
hhist=zeros(1,niter);
nhist=zeros(1,niter);


for iter=1:niter
gna=gnamax*m^3*h;
gk=gkmax*n^4;
gtot=gna+gk+gl;
vinf = ((gna*vna+gk*vk+gl*vl)+ iapp(iter))/gtot;
tauv = cm/gtot;
v=vinf+(v-vinf)*exp(-dt/tauv);
alpham = 0.1*(v+40)/(1-exp(-(v+40)/10));
betam = 4*exp(-0.0556*(v+65));
alphan = 0.01*(v+55)/(1-exp(-(v+55)/10));
betan = 0.125*exp(-(v+65)/80);
alphah = 0.07*exp(-0.05*(v+65));
betah = 1/(1+exp(-0.1*(v+35)));
taum = 1/(alpham+betam);
tauh = 1/(alphah+betah);
taun = 1/(alphan+betan);
minf = alpham*taum;
hinf = alphah*tauh;
ninf = alphan*taun;
 m=minf+(m-minf)*exp(-dt/taum);
 h=hinf+(h-hinf)*exp(-dt/tauh);
 n=ninf+(n-ninf)*exp(-dt/taun);
vhist(iter)=v; mhist(iter)=m; hhist(iter)=h; nhist(iter)=n;
gnahist(iter) = gna;
gkhist(iter) = gk;
end;

j=1;
real_peaks=zeros;
[peaks, locs]=findpeaks(vhist);%Storing the porential value for each peak for each value of I_ext

for temp=1:length(peaks)
if peaks(temp) >=12 % (Thresholding it by 12) minimum potential value at which a waveform is considered AP.
real_peaks(j)=peaks(temp);
j=j+1;
end;
end;

if real_peaks ~= 0
no_of_peaks(k)=length(real_peaks);
else
no_of_peaks(k)=0;
end;

k=k+1
end;

figure(1)
%subplot(2,1,1)
plot(t,vhist)
title('Voltage vs Time')

figure(2)
%subplot(2,1,2)
plot(t,mhist,'y-', t,hhist,'g.',t,nhist,'b-')
legend('m','h','n')

figure(3)
gna=gnamax*(mhist.^3).*hhist;
gk=gkmax*nhist.^4;
clf
plot(t,gna,'r');
hold on
plot(t,gk,'b');
legend('gna','gk')
hold off

%ploting the behaviour of no. of action potential for the given external current
figure(4);
X=0:istep:0.8;
plot(X,no_of_peaks*1000/(niter/100));
xlabel('External Current (I_{ext})');
xlim([0 1.2]);
ylabel('Firing (Discharge) Rate')
hold on;

for l=2:length(no_of_peaks) %to define I1, I2, I3.
if no_of_peaks(l)>0 && no_of_peaks(l-1)==0%getting the I_ext I1 till where no action potential is observed at starting
 I1=(l-1)*istep
end;
if no_of_peaks(l)>no_of_peaks(l-1)+5%getting the I_ext I2 till where the finite no. of action potential is observed
 I2=(l-1)*istep
end;
if no_of_peaks(l)<no_of_peaks(l-1)-5%getting the I_ext I3 till where the infinite no. of action potential is observed and after that no AP
 I3=(l-1)*istep
end;
end;

I1=0*no_of_peaks*1000/(niter/100)+I1;
plot(I1,no_of_peaks*1000/(niter/100),'c');
text(I1(1),-3,'I1');

I2=0*no_of_peaks*1000/(niter/100)+I2;
plot(I2,no_of_peaks*1000/(niter/100),'y');
text(I2(1),-3,'I2');

I3=0*no_of_peaks*1000/(niter/100)+I3;
plot(I3,no_of_peaks*1000/(niter/100),'m');
text(I3(1),-3,'I3');

text(0.6,90,'Before I1 => No A.P');
text(0.6,85,'From I1 to I2 => Finite no. of A.P.s');
text(0.6,80,'From I2 to I3 => Infinite no.of A.P.s');
text(0.6,75,'After I3 => No A.P.s');
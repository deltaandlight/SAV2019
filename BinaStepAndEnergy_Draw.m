filename='BinaStepAndEnergy';
dat1=importdata(filename,' ',4);
N=32;
A=dat1;
dtime=A(:,1);
energy=A(:,2);
time=dtime+0;
plot(time,energy);
saveas(gcf,['StepAndEnergy','.jpg']);
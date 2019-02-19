filename='BinaStepAndEnergy';
dat1=importdata(filename,' ',1);
N=3200;
A=dat1.data;
dtime=A(:,1);
energy=A(:,2);
time=dtime+0;
for i=1:N
    time(i)=i*dtime(1);
end
plot(time,energy);
saveas(gcf,['StepAndEnergy','.jpg']);
L=128;
%L=256;
x=0:(1/(L-1)):1;
y=0:(1/(L-1)):1;
hx=2*pi/L;
hy=2*pi/L;

filename='ASCIsol_0016';
dat1=importdata(filename,' ',4);
A=reshape(dat1.data,L,L);
image(A,'CDataMapping','scaled');
saveas(gcf,['myfig_image',16,'.jpg']);

filename='ASCIsol_0032';
dat1=importdata(filename,' ',4);
A=reshape(dat1.data,L,L);
image(A,'CDataMapping','scaled');
saveas(gcf,['myfig_image',32,'.jpg']);
    


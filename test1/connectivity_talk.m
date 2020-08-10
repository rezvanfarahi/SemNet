bpFilt = designfilt('bandstopfir','FilterOrder',1000, ...
         'CutoffFrequency1',13,'CutoffFrequency2',31, ...
         'SampleRate',1000);
fvtool(bpFilt)
dataIn = x12;
x12nb = filter(bpFilt,dataIn);

bpFilt = designfilt('bandpassfir','FilterOrder',1000, ...
         'CutoffFrequency1',13,'CutoffFrequency2',31, ...
         'SampleRate',1000);
fvtool(bpFilt)
dataIn = x1;
x1b = filter(bpFilt,dataIn);

x12wb=x12nb+x1b;%circshift(x1b,100);

bpFilt = designfilt('bandpassfir','FilterOrder',1000, ...
         'CutoffFrequency1',13,'CutoffFrequency2',31, ...
         'SampleRate',1000);
fvtool(bpFilt)
dataIn = x12wb;
x12b = filter(bpFilt,dataIn);

x1rb=x1b(randperm(length(x1b)));
mscohere(x1b,x12b,hanning(1024),512,1024,1000)
corrcoef(x1b,x12b)
plot(x1b(5000:10000))
hold on, plot(x12b(5000:10000),'r')

x1=rds(5000:10000,1);
plot(x1(2:end),x1(1:end-1),'.')
[p,S,mu]=polyfit(x1(1:end-4),x1(5:end),1);
plot(x1(1:end-4),x1(5:end),'.')
[y1,delta]=polyval(p,x1(1:end-4),S,mu);
hold on
plot(x1(1:end-4),y1,'r')
hold on
plot(x1(1:end-4),y1,'r')
ylabel('x(t)')
xlabel('x(t-\tau)')
x2=rds(5000:10000,4);
figure, plot3(x1(1:end-4),x2(1:end-4),x1(5:end),'.')
xlabel('x(t-\tau)')
ylabel('y(t-\tau)')
zlabel('x(t)')
grid on
[b,bint,r]=regress(x1(5:end),[x1(1:end-4),x2(1:end-4)]);
[X,Y] = meshgrid(-600:1:600);
z=b(1)*X+b(2)*Y;
C=zeros(1601,1601,3);
C(:,:,1)=1;
hold on
mesh(X,Y,z,C)


subplot(5,1,1), plot(100:0.01:105,x1(5000:10:10000))
ylabel('x1')
subplot(5,1,2), plot(100:0.01:105,x4(5000:10:10000))
ylabel('x2')
subplot(5,1,3), plot(100:0.01:105,x8(5000:10:10000))
ylabel('x3')
subplot(5,1,4), plot(100:0.01:105,x12(5000:10:10000))
ylabel('x4')
x1f=50-x1.^2/1000+x12;
subplot(5,1,5), plot(100:0.01:105,x1f(5000:10:10000))
ylabel('x5')
xlabel('time (s)')

figure,
subplot(4,1,1), plot(x1(5000:10:10000),x4(5000:10:10000),'.')
[p,S,mu]=polyfit(x1,x4,1);
[y1,delta]=polyval(p,x1,S,mu);
hold on
subplot(4,1,1),plot(x1(5000:10:10000),y1(5000:10:10000),'r')
ylabel('x2')

subplot(4,1,2), plot(x1(5000:10:10000),x8(5000:10:10000),'.')
[p,S,mu]=polyfit(x1,x8,1);
[y1,delta]=polyval(p,x1,S,mu);
hold on
subplot(4,1,2),plot(x1(5000:10:10000),y1(5000:10:10000),'r')
ylabel('x3')

subplot(4,1,3), plot(x1(5000:10:10000),x12(5000:10:10000),'.')
[p,S,mu]=polyfit(x1,x12,1);
[y1,delta]=polyval(p,x1,S,mu);
hold on
subplot(4,1,3),plot(x1(5000:10:10000),y1(5000:10:10000),'r')
ylabel('x4')


subplot(4,1,4), plot(x1(5000:10:10000),x1f(5000:10:10000),'.')
[p,S,mu]=polyfit(x1,x1f,1);
[y1,delta]=polyval(p,x1,S,mu);
hold on
subplot(4,1,4),plot(x1(5000:10:10000),y1(5000:10:10000),'r')
ylabel('x5')
xlabel('x1')


t=0.01:0.01:12*pi;
x11=0.2*sin(t);
x21=0.2*sin(t+pi/6);
figure, subplot(2,1,1), plot(linspace(0,2,length(t)),x11)
ylabel('ROI1')
subplot(2,1,2), plot(linspace(0,2,length(t)),x21)
ylabel('ROI2')
xlabel('time (s)')

t=0.01:0.01:12*pi;
x11=0.2*sin(t+pi/4);
x21=0.2*sin(t+5*pi/12);
figure, subplot(2,1,1), plot(linspace(0,2,length(t)),x11)
ylabel('ROI1')
subplot(2,1,2), plot(linspace(0,2,length(t)),x21)
ylabel('ROI2')
xlabel('time (s)')

x11=0.2*sin(t-pi/7);
x21=0.2*sin(t);
figure, subplot(2,1,1), plot(linspace(0,2,length(t)),x11)
ylabel('ROI1')
subplot(2,1,2), plot(linspace(0,2,length(t)),x21)
ylabel('ROI2')
xlabel('time (s)')


t=0.01:0.01:12*pi;
x1=0.5*sin(t);
x2=0.78*sin(2*t+pi/6);
x3=0.9*sin(7*t+pi/2);
x4=1.07*sin(3*t-pi/3);
x5=0.33*sin(0.5*t-pi/4);
subplot(5,1,1), plot(linspace(0,2,length(t)),x5)
subplot(5,1,2), plot(linspace(0,1.5,length(t)),x1)
subplot(5,1,3), plot(linspace(0,1,length(t)),x2)
subplot(5,1,4), plot(linspace(0,1,length(t)),x4)
subplot(5,1,5), plot(linspace(0,1,length(t)),x3)
xlabel('time (s)')
figure, plot(linspace(0,2,length(t)),x1+x2+x3+x4+x5)
xlabel('time (s)')


subplot(3,1,1), plot(t*10,x21)
ylabel('ROI1')
subplot(3,1,2), plot(t*10,x11)
ylabel('ROI2')
subplot(3,1,3), plot(t*10,x31)
ylabel('ROI3')
xlabel('time (ms)')
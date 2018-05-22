close all;

sampleRate = 12; %average
start = 42;
length = 100;
x = (start:1/sampleRate:(start+length))';
x = x + (1/(2*sampleRate+1))*rand(size(x));%randomization of sample rate with monotonically increasing garuntee
y = .1*sin(2*pi*x/100) + .2*sin(2*pi*x/10) + sin(2*pi*x) + .4*randn(size(x));
ts = [x y];

am = spAgeModelFitter(ts);
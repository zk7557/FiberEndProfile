% units are in um 
c = 2.9979*10^14;
a = 4.5;
lambda = 3;
k0 = 2 * pi / lambda;
n = 1.5078;
n2 = 1.4945;
kappa = k0 * sqrt(n^2-n2^2);

dr = 0.01;
r1 = 0:dr:a-dr;
r2 = a:dr:2*a-dr;
r = [r1 r2];
rN = 2*a/dr;

pFiber= [besselj(0, kappa * r1) (sqrt(a./r2) .* ... 
       besselj(0, kappa * a) .* exp(-(r2-a)/sqrt(2)))];

% figure(1); plot(r, pFiber);

dTheta = pi/1000;
theta = -pi/4:dTheta:pi/4;

dx = 0.001;
x = -2*a:dx:2*a;
x2 = x.^2;
N = 4*a/dx +1;
R = x2 .* ones(N, 1) + x2' .* ones(1, N);
R = sqrt(R);

% figure(2); imagesc(R);
boundary = besselj(0, kappa * a);
intensity = zeros(N);
for i = 1:N
    for j = 1:N
%         intensity(i,j) = piecewise(R(i,j)<a, besselj(0, kappa * r1), ...
%             R(i,j)>=a, sqrt(a/R(i,j)) * boundary * exp(-(r2-a)/0.7));
        if R(i,j) < a
            intensity(i,j) = besselj(0, kappa * R(i,j));
        else
            intensity(i,j) = sqrt(a/R(i,j))*boundary*exp(-(R(i,j)-a)/0.7);
        end
    end
end
intensityX = sum(intensity);
sinTheta = sin(theta);
phase = 2i * pi / lambda * (x .* sinTheta');
field = intensityX .* exp(phase);
profile = sum(field, 2);
profile = abs(profile').^2;

figure(3); plot(theta, profile);

% compare with gaussian beam 
w = lambda/pi/a;
figure(4); plot(theta, exp(-(theta.^2)/2/(w^2)));

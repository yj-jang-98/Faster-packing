clear all;
clc;
close all;

% four tank system
A1 = 28; A2 = 32; A3 = 28; A4 = 32;
a1 = 0.071; a2 = 0.057; a3 = 0.071; a4 = 0.057;
kc = 0.5;
g = 981;
h10 = 12.4;
h20 = 12.7;
h30 = 1.8;
h40 = 1.4;
v10 = 3;
v20 = 3;
k1 = 3.33;
k2 = 3.35;
g1 = 0.7;
g2 = 0.6;

A0 = [-a1/A1 * 2*g / (2*sqrt(2*g*h10)), 0 , a3/A1 * 2*g / (2*sqrt(2*g*h30)), 0;...
      0, -a2/A2 * 2*g / (2*sqrt(2*g*h20)), 0 , -a4/A2 * 2*g / (2*sqrt(2*g*h40));...
      0, 0, -a3/A3 * 2*g / (2*sqrt(2*g*h30)), 0 ;...
      0, 0, 0, -a4/A4 * 2*g / (2*sqrt(2*g*h40))];
B0 = [g1*k1/A1, 0;...
      0, g2*k2/A2;...
      0, (1-g2)*k2/A3;...
      (1-g1)*k1/A4,0];

C = [kc, 0, 0, 0;...
     0, kc, 0, 0];

Ts = 0.1; % Sampling time

sysC = ss(A0,B0,C,[]);
sysD = c2d(sysC, Ts);

A = sysD.A
B = sysD.B

% find K and L
Q = eye(4);
R1 = eye(2);
R2 = eye(2);
[~, K, ~] = idare(A,B,Q,R1,[],[]);
% K = place(A,B,[0.25,0.35+0.05i, 0.35-0.05i, 0.45, 0.49])
K = -K;
abs(eig(A+B*K))
eig(A+B*K)
[~, L, ~] = idare(A.', C.', Q, R2, [], []);
L = L.';
abs(eig(A-L*C))

F = A+B*K-L*C;
G = L;
H = K;

R = place(F.',H.',[0,1,2,-1]);
R = R.';
eig(F-R*H)

sys = ss(F-R*H, G, H, []);
[csys,T] = canon(sys, 'modal');
F_ = T*(F-R*H)/T
R_ = T*R
G_ = T*G
H_ = H/T

%% For simulation check
Fnorm = norm(F_,inf)
Gnorm = norm([G_,R_], inf)
Hnorm = norm(H_, inf)
Bnorm = norm(B, inf)
Acl = [A, B*H_ ; G_*C, F_+R_*H_];
abs(eig(Acl))

r = 1e-2; s = 1e-4; L = 1e-2;
M = 1;
lambda = 0.98;
sigma = 19.2;
N = 2^12;
nu = 2^7;
d = 8;
sigma_mult = 2*d*N*sigma*nu;
alpha = r*s*L*12*sigma_mult + r * Gnorm/2 + r*L*Gnorm*sigma
beta = r*s^2 * L * 5 * sigma_mult
gamma = r*s*L*sigma
eta = 1 + gamma + (Bnorm*beta + alpha)/(1-lambda)
qbound = 2*max(eta/(r*s*L), (Hnorm * eta + beta)/r*s^2*L)
epsilonBound = max(2*beta, (eta-1)*2*Hnorm)

alpha_ = alpha + r*s*L*5*(Fnorm + Gnorm/s)*3*sigma_mult
beta_ = beta + r*s*L*2*Hnorm*3*sigma_mult
eta_ = 1 + gamma + (Bnorm*beta_ + alpha_)/(1-lambda)
qbound_ = 2*max(eta_/(r*s*L), (Hnorm * eta_ + beta_)/r*s^2*L)
epsilonBound_ = max(2*beta_, (eta_-1)*2*Hnorm)


%%
scale = 10000;
G_ = scale * G_;
R_ = scale * R_;
H_ = scale * H_;

iter = 500;
% xplant = rand(5,1);
% xcont = rand(5,1);
xplant = [1;1;1;1];
xcont = [1.5;-1;2.5;0]*scale;

Y = [];
Xplant = [];
Xcont = [];
U = [];

for i = 1:iter
    y = C*xplant;
    u = H_*xcont/(scale^2);
    xplant = A*xplant+B*u;
    xcont = F_*xcont + G_*y + R_*u;
    % xcont =F_*xcont;
    % xcont = G_*y;
    % xcont = R_*u;
    Y = [Y, y];
    U = [U, u];
    Xplant = [Xplant,xplant];
    Xcont = [Xcont, xcont/scale];
end

figure(1)
subplot(221);
plot(Xplant)
subplot(222);
plot(Xcont)
subplot(223);
plot(Y)
subplot(224);
plot(U)












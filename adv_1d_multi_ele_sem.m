clear;
clc;

close all

N = 5;
E = 10;
Lx = 2;
Le = Lx/E;
J = 2/Le;

[z, w, p_gll, Dh] = gen_zwd(N);

c = 1;

Bh = diag(w);

% Initial conditions
xinit = std_element_mapping(E, N, Lx, z);
xinit = unique(xinit);

nu = 0.01;
init_coeff = 1/(4*nu);
uL = exp(-init_coeff*(xinit-0.5).^2);
uL = uL(1:end-1);

u_analytic_o = uL;
u_analytic = exp(-init_coeff*(xinit-1.5).^2);
u_analytic = u_analytic(1:end-1);

t = 0;
tE = 1;
dt = 1e-5;

Ie = eye(E);
[R, Rp] = gen_r(E, N);
Q = Rp;
A = kron(Ie, Bh*Dh);
A = Q'*A*Q*J;

M = kron(Ie, Bh);
M = Q'*M*Q;

A = M\A;
% A = A/2;

while t < tE + dt
    t = t + dt;
    
    dudt = -c*A*uL;

    uL = uL + dudt*dt;


end
plot(xinit(1:end-1),uL)

norm1 = norm(u_analytic - uL,1)/(E*N)

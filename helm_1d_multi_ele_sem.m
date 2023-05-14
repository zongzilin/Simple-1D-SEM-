clear;
clc;

close all

N =9;
E = 50;
Lx = 2;
Le = Lx/E;
J = 2/Le;

[z, w, p_gll, D] = gen_zwd(N);
%     [Ahh,Bhh,Chh,Dhh,zz,ww] =  semhat(N);
%     D = Dh;
Q = semq(E, N, 99);

x = std_element_mapping(E, N, Lx, z);
x = x - 1;

lam = 1;
f = (pi^2 + lam)*cos(pi.*x);
% f = ((1-x).^2).*((1+x).^2).*exp(2*x);

A = D'*diag(w)*D.*J^2 + lam*diag(w);

Ie = eye(E);
Ag = kron(Ie, A);
Bg = kron(Ie, diag(w));

Fg = Bg*f;
Fg = diag(Fg);

Ag = Q'*Ag*Q;
Fg = Q'*Fg*Q;
Fg = diag(Fg);

Fg(1) = cos(-pi);
Fg(end) = cos(pi);
Ag(1,:) = 0;
Ag(1,1) = 1;
Ag(end,:) = 0;
Ag(end,end) = 1;

u = Ag\Fg;
u = Q*u;
realU = cos(pi.*x);
sda = norm(u - realU)/(E*(N+1))

plot(x,realU,'-ko','MarkerIndices',1:25:length(realU),'Markersize',10,'linewidth',1.5)
hold on
plot(x,u,'-kx','MarkerIndices',10:20:length(u),'Markersize',10,'linewidth',1.5)

grid on
grid minor

xlabel('$x$','interpreter','latex')
ylabel('$u$','interpreter','latex')

% legend('Analytical','Spectral-elemet','interpreter','latex')

set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',25)
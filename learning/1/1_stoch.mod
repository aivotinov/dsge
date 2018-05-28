var c k n y w z tau i;
varexo e_tau e_a;

parameters alpha beta theta delta tauHat rho r phiSS omegaSS muSS;

model;
1/c = beta * (1 + r)/c(+1) * (alpha * k^(alpha - 1) * exp((1-alpha)*z(+1)) * n(+1)^(1-alpha) + 1 - delta);
c * theta / ((1 - n)*(1- tau)) = (1-alpha) * k(-1)^(alpha) * exp((1-alpha)*z) * n^(-alpha);
c + i = y;
y = k(-1)^(alpha)*(exp(z)* n)^(1-alpha);
i = k-(1-delta)*k(-1);
z = rho * z(-1) + e_a;
tau = tauHat + e_tau;
w = (1-alpha) * k(-1)^(alpha) * (exp(z)*n)^(1-alpha);
end;

alpha   = 0.4;
beta    = 0.99;
delta   = 0.02388;
theta     = 1.75;
tauHat = 0.13;
rho     = 0.95;
sigma   = (0.007/(1-alpha));
r = 1/0.99-1;
phiSS = ((1/beta-1+delta)/alpha)^(1/(1-alpha));
omegaSS = (phiSS)^(1-alpha)-delta;
muSS = (1-alpha)/(theta*(1-tauHat))*phiSS^(-alpha);

initval;
  k = muSS/(omegaSS + muSS * phiSS);
  c = omegaSS * k;
  n = phiSS * k;
  z = 0;
  e_tau = 0;
  e_a = 0;
end;
steady;


shocks;
var e_tau = sigma^2;
end;


check;

stoch_simul(order = 1);
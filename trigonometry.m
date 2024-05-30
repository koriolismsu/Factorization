clear all

theta_in = rand*2*pi;
phi = rand*2*pi;
psi = rand*2*pi;


val1 = sin(psi)*cos(phi+theta_in)-cos(theta_in)*sin(psi+phi);

val2 = sin(theta_in)*cos(psi+phi)-sin(theta_in+phi)*cos(psi);

delta = val2-val1
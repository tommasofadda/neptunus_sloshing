
function D = drag(U)

global Ab Cd rho_atm N

D = -1/2*rho_atm*N*Ab*Cd.*U.^2;

end
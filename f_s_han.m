
function f = f_s_han(Rn,Rb,Ln)

k = 110;
f = @(x) (Rn+(Rb-Rn)./(1+exp(-k*(x-Ln/2)))-(Rb-Rn)/(1+exp(k*Ln/2))-Rn*(1-(1/(1+exp(-k*Ln/2))-1/(1+exp(k*Ln/2)))))/(1./(1+exp(-k*Ln/2))-1./(1+exp(k*Ln/2)));
end

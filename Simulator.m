% Simulator 2.0 - First complete version available

global Ab Cd rho_atm An rho_w glob_A N

%% Fluids characteristics (water and atmoshphere)

Cv = 718;
Cp = 1005;
k = Cp/Cv;
R = Cp - Cv;

rho_atm = 1.225;          % air density [kg/m^3]
rho_w   = 1000;           % water density [kg/m^3]
g       = 9.80665;        % gravity Acc. [m/s^2]
Patm    = 1 * 101325;     % atmospheric pressure [kg/(m*s^2)]

%% Propulsive unit configuration --> Choose number of boosters and pressure

% Number of boosters
N  = 2;

% Starting pressure
P0a     = 7 * 101325;        % air pressure inside the bottle before launch
% Starting temperature
T0a     = 20 + 273.15;       % air temperature before launch = ambient temp.
rho0a   = P0a/(R*T0a);       % density of the compressed air (considering perfect gas law)

%% Propellant bottle characteristics --> Choose dimension and filling

% Bottle capacity     [1000ml 1500ml 1750ml 2000ml 3000ml]
%                     [  1      2      3      4      5   ]
bottle_capacity =     4;

V_bottle_span   =     [1.000, 1.500, 1.750, 2.000, 3.000];
M_bottle_span   =     [0.030, 0.040, 0.050, 0.054, 0.080];

Vb              = V_bottle_span(bottle_capacity) * 10^(-3); % bottle volume [m^3]
Mb              = M_bottle_span(bottle_capacity);           % bottle mass   [kg]

% Filling factor (from 0 to 1)
f_fill  = 0.16;             % filling factor

V0w     = Vb*f_fill;         % inital volume occupied by water
V0a     = Vb-V0w;            % inital volume occupied by air

% Bottle dimensions calculated w.r.t. to 1.5L bottle
Rb      = 0.0428*(Vb/(1.5*10^(-3)))^(1/3);   % bottle radius
Rn      = 0.011;                             % nozzle exit section, viz. nozzle radius
Ln      = 0.10*(Vb/(1.5*10^(-3)))^(1/3);     % nozzle length
Ab      = (pi*Rb^2)*N;                       % bottles area seen by the air

% Bottle length calculation from chosen volume, considering the nozzle
% sigmoid shape
f_nozzle_aux = f_s_han(Rn,Rb,Ln);       %Sigmoid shape for the nozzle
r_auxil      = @(x) f_nozzle_aux(x);
A_auxil      = @(x) pi*r_auxil(x)*r_auxil(x);
Vnozzle      = integral(@(x) A_auxil(x), 0, Ln,'ArrayValued',true);

%Bottle length
Lb           = (Vb-Vnozzle)/Ab+Ln;

%% Rocket masses --> Choose structure mass (dry weight without propellant bottles empty weight)

% Structure mass
Ms  = 1;                  % mass overall structure

M0w = rho_w*V0w;         % initial water mass inside the single bottle
M0a = rho0a*V0a;         % initial air mass inside the single bottle

SloshingMass = N*M0w;    % sloshing mass equal to the total propellant mass
if SloshingMass < 0.5
    SloshingMass = 0.5;
end

% Total Rocket mass
M0  = Ms + (Mb+M0w+M0a)*N + SloshingMass; % initial mass of the s/c

%% Launch tube dimensions --> choose launch tube lenght (lower than bottle length)

TL  = 0.2;                % tube length
TOD = 0.95*(2*Rn);       % tube outer diameter
TID = TOD - 2*0.002;     % tube inner diameter
TOA = pi*(TOD/2)^2;      % tube outer area
TIA = pi*(TID/2)^2;      % tube inner area

LauncherGasV = 0.001;    % additional gas inside launch tube but outside the bottle

%% Drag coefficient - to be calculated through CFD or statistics
Cd = 0.5;                % drag coefficient

%% Some formulas for bottle dimensions to use during thrust calculation
%Bottle nozzle shape
f_nozzle = f_s_han(Rn,Rb,Ln);
%Bottle radius at a given axial position
r = @(x) (x(1)>=0 & x(1)<=Ln)*f_nozzle(x) + (x(1)>Ln & x(1)<=Lb)*Rb;

% Plot bottle shape
% figure('Name','Bottle shape')
% fplot(y,'b')
% hold on
% fplot(-y,'b')
% axis('equal')
% set(gcf,'position',[489 343 860 320])
% xlim([-0.015 Lb*1.05]); ylim([-Rb*1.3 Rb*1.3])
% xhead = [0 0];
% yhead = [-Rn Rn];
% line(xhead,yhead,'Color','blue','LineStyle','-')
% xbottom = [Lb Lb];
% ybottom = [-Rb Rb];
% line(xbottom,ybottom,'Color','blue','LineStyle','-')
% hold off
% clear xhead xbottom yhead ybottom

A        = @(x) pi*r(x(1))*r(x(1));          % overall area function
glob_A   = A;
% Ainv   = @(x) ( pi.*r(x(1)).^2 ).^(-1);    % inverse of area function
Ainv     = @(x) A(x).^-1;
An       = feval(A,0);
Vw       = @(x) integral(@(x) A(x), 0, x(1),'ArrayValued',true);
%findH0   = @(x) integral(@(x) A(x), 0, x(1),'ArrayValued',true)- V0w;

%% Phase 0: Launch - Tube Phase

t0P0   = 0;
tMaxP0 = 2;
dtP0   = 1e-5;
tP0    = t0P0:dtP0:tMaxP0;
tOutP0 = 0;

ViP0   = LauncherGasV + Vb - V0w - TL*(TOA - TIA);
nP0    = length(tP0);

ThrustP0       = zeros(nP0,1);
PaP0           = zeros(nP0,1);
RHOaP0         = zeros(nP0,1);
TaP0           = zeros(nP0,1);
RocketMTotalP0 = zeros(nP0,1);
RocketAccP0    = zeros(nP0,1);
RocketVP0      = zeros(nP0,1);
RocketYP0      = zeros(nP0,1);

for i = 1:nP0
    if i ~= 1
        RHOaP0(i) = rho0a * ViP0 / (ViP0 + RocketYP0(i-1) * (TOA-TIA));    % TOA instead of TIA
    else
        RHOaP0(i) = rho0a;
    end
    PaP0(i)           = P0a * ( RHOaP0(i) / rho0a ) .^ k;
    TaP0(i)           = T0a * ( RHOaP0(i) / rho0a ) .^ (k-1);
    
    ThrustP0(i)          = N * ( PaP0(i) - Patm ) * TOA;   % TOA instead of TIA
    RocketMTotalP0(i)    = M0;
    
    if i ~= 1
        RocketAccP0(i)  = ( ThrustP0(i) + drag(RocketVP0(i-1)) ) /  RocketMTotalP0(i) - g;
        RocketVP0(i)    = RocketVP0(i-1) + RocketAccP0(i) * dtP0;
        RocketYP0(i)    = RocketYP0(i-1) + RocketVP0(i) * dtP0;
    else
        RocketAccP0(i)  = ThrustP0(i) / RocketMTotalP0(i) - g;
    end
    
    if RocketYP0(i) >= TL
        break
    end
    
    tOutP0 = tOutP0 + 1;
end

ThrustP0       = ThrustP0(1:tOutP0);
PaP0           = PaP0(1:tOutP0);
RHOaP0         = RHOaP0(1:tOutP0);
TaP0           = TaP0(1:tOutP0);
RocketMTotalP0 = RocketMTotalP0(1:tOutP0);
RocketAccP0    = RocketAccP0(1:tOutP0);
RocketVP0      = RocketVP0(1:tOutP0);
RocketYP0      = RocketYP0(1:tOutP0);

Pa01       = PaP0(end);
Rhoa01     = RHOaP0(end);
RocketV01  = RocketVP0(end);
RocketY01  = RocketYP0(end);
RocketAcc01 = RocketAccP0(end);

%% Phase 1: Water - Impulse Phase

% options = optimset;            % Optimization Toolbox needed
% options.Display = 'off';
% options.TolFun = 1e-12;
% options.MaxFunEvals = 10000;

a_bis=0;
b_bis=Lb;
fa_bis = integral(@(x) A(x), 0, a_bis,'ArrayValued',true)- V0w;
fb_bis = integral(@(x) A(x), 0, b_bis,'ArrayValued',true)- V0w;

tol_bis=1e-12;
maxit=ceil(log((b_bis-a_bis)/tol_bis)/log(2));
it=0;

% Bisection method to find H0 - Initial water height inside bottle
while it<=maxit && abs(b_bis-a_bis)>=tol_bis+eps*max([abs(a_bis) abs(b_bis)])
    it=it+1;
    xk(it)=a_bis+(b_bis-a_bis)*0.5;
    fxk= integral(@(x) A(x(1)), 0, xk(it),'ArrayValued',true)- V0w;
    if fxk==0
        break;
    elseif sign(fxk)*sign(fa_bis)>0
        a_bis=xk(it);
        fa_bis=fxk;
    elseif sign(fxk)*sign(fb_bis)>0
        b_bis=xk(it);
        fb_bis=fxk;
    end
end
H0=xk(it);



fBP1 = @(x) An * integral(@(x) Ainv(x), 0, x(1),'ArrayValued',true);
fCP1 = @(x) 0.5 * ( ( An / A(x) )^2 - 1 );
fDP1 = @(x) ((Vb-V0w)/(Vb-Vw(x)))^k;


fRHOaP1 = @(x) Rhoa01*(fDP1(x))^(1/k);
fMTP1 = @(x) Mb*N + Ms + SloshingMass + N*Vw(x)*rho_w + N*Rhoa01*(Vb-V0w); %Qua ho della roba inutile, si semplifica tutto

% System of ODEs
% [x(1),x(2)] = [H, Un]
fP1 = @(t,x) [(An/A(x))*x(2); 1/(fBP1(x)-rho_w*An*x(1)^2/fMTP1(x))*(Patm/rho_w-fDP1(x)*Pa01/rho_w-x(1)*drag(x(2))/fMTP1(x)-x(2)^2*(fCP1(x)+x(1)*rho_w*An/fMTP1(x)*(1-An/A(x))))];

dtP1 = 1e-5;

tspan = [0,4];

x0 = [H0;0];

options=odeset('RelTol', 1e-6);

[tP1,solP1] = ode23tb(fP1,tspan,x0,options);

HP1 = solP1(:,1);

tOutP1 = 1;
while HP1(tOutP1) > 0.00001
    tOutP1 = tOutP1 + 1;
end

tP1_bis = 0:dtP1:tP1(tOutP1);
solP1 = interp1(tP1,solP1,tP1_bis);
tP1 = tP1_bis;
tOutP1 = length(tP1);

HP1 = solP1(:,1);
UnP1 = solP1(:,2);


ThrustP1       = zeros(tOutP1,1);
MassFlowP1     = zeros(tOutP1,1);
PaP1           = zeros(tOutP1,1);
RHOaP1         = zeros(tOutP1,1);
TaP1           = zeros(tOutP1,1);
RocketMTotalP1 = zeros(tOutP1,1);
RocketAccP1    = zeros(tOutP1,1);
RocketVP1      = zeros(tOutP1,1);
RocketYP1      = zeros(tOutP1,1);
dUndtP1        = zeros(tOutP1,1);

RocketVP1(1)   = RocketV01;
RocketYP1(1)   = RocketY01;
RocketAccP1(1) = RocketAcc01;
%
%capiamo come funziona gradient
dUndtP1        = gradient(UnP1, dtP1);

for i = 1:tOutP1
    %va calcolata *N la thrust?
    ThrustP1(i)       = rho_w * An * UnP1(i)^2;
    MassFlowP1(i)     = rho_w * An * UnP1(i);
    RHOaP1(i)         = Rhoa01 * fDP1(HP1(i))^(1 / k);
    PaP1(i)           = Pa01 * ( (Vb - V0w) / (Vb - Vw(HP1(i))))^k;
    TaP1(i)           = PaP1(i) / ( R * RHOaP1(i) );
    RocketMTotalP1(i) = fMTP1(HP1(i));
    fIntP1 = -rho_w*An*(HP1(i)*dUndtP1(i)+(An/A(HP1(i)))*UnP1(i)^2);
    
    if i >= 2
        RocketAccP1(i)    = (N*ThrustP1(i) + N*fIntP1 + drag(RocketVP1(i-1)))/RocketMTotalP1(i)-g;
        RocketVP1(i)    = RocketVP1(i-1) + RocketAccP1(i) * (dtP1);
        RocketYP1(i)    = RocketYP1(i-1) + RocketVP1(i)   * (dtP1);
    end
end
Pa02           = PaP1(end);
RHOa02         = RHOaP1(end);
Ta02           = TaP1(end);
RocketMTot02   = RocketMTotalP1(end);
RocketThrust02 = ThrustP1(end);
RocketAcc02    = RocketAccP1(end);
RocketV02      = RocketVP1(end);
RocketY02      = RocketYP1(end);

%% Phase 2: Air Blowdown Phase

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Regime 1 --> Chocked Flow


%System of ODEs
% RHOt = x(1)

fP2_Chock = @(t,x) -An/Vb*((k+1)/2)^(1/(1-k)-0.5)*x(1)*(R*k*Ta02)^(0.5)*(x(1)/RHOa02)^((k-1)/2);

dtP2_Chock = 1e-5;

tspan = [0,2];
X0P2_Chock = RHOa02;
options = odeset('RelTol',1.e-6);

[tP2_Chock,solP2_Chock] = ode45(fP2_Chock,tspan,X0P2_Chock,options);

nTSP2_Chock  = length(tP2_Chock);
RHOaP2_Chock = solP2_Chock;
PaP2_Chockaux = zeros(nTSP2_Chock,1);

for i = 1:nTSP2_Chock
    PaP2_Chockaux(i) = Pa02*(RHOaP2_Chock(i)/RHOa02)^k;
end

tOutP2_Chock = 1;
while PaP2_Chockaux(tOutP2_Chock) >= Patm * ((k+1)/2)^(k/(k-1))
    tOutP2_Chock = tOutP2_Chock + 1;
end

tP2_Chock_bis = 0:dtP2_Chock:tP2_Chock(tOutP2_Chock);
solP2_Chock = interp1(tP2_Chock,solP2_Chock,tP2_Chock_bis);
tP2_Chock = tP2_Chock_bis;
tOutP2_Chock = length(tP2_Chock);

nTSP2_Chock  = length(tP2_Chock);
RHOaP2_Chock = solP2_Chock;

if tOutP2_Chock ~= 1
    
    ThrustP2_Chock       = zeros(tOutP2_Chock,1);
    MassFlowP2_Chock     = zeros(tOutP2_Chock,1);
    PaP2_Chock           = zeros(tOutP2_Chock,1);
    TaP2_Chock           = zeros(tOutP2_Chock,1);
    PnP2_Chock           = zeros(tOutP2_Chock,1);
    RHOnP2_Chock         = zeros(tOutP2_Chock,1);
    TnP2_Chock           = zeros(tOutP2_Chock,1);
    UnP2_Chock           = zeros(tOutP2_Chock,1);
    MnP2_Chock           = zeros(tOutP2_Chock,1);
    RocketMTotalP2_Chock = zeros(tOutP2_Chock,1);
    RocketAccP2_Chock    = zeros(tOutP2_Chock,1);
    RocketVP2_Chock      = zeros(tOutP2_Chock,1);
    RocketYP2_Chock      = zeros(tOutP2_Chock,1);
    
    RocketVP2_Chock(1)   = RocketV02;
    RocketYP2_Chock(1)   = RocketY02;
    
    %
    RocketVP2_Chock(1)   = RocketV02;
    RocketYP2_Chock(1)   = RocketY02;
    RocketAccP2_Chock(1) = RocketAcc02;
    
    for i = 1:tOutP2_Chock
        PaP2_Chock(i)           = Pa02 * ( RHOaP2_Chock(i) / RHOa02 )^(k);
        TaP2_Chock(i)           = Ta02 * (RHOaP2_Chock(i) / RHOa02) ^ (k - 1.0);
        RHOnP2_Chock(i)         = RHOaP2_Chock(i) * ( (k + 1.0) / 2.0 )^(1.0 / (1.0 - k));
        PnP2_Chock(i)           = PaP2_Chock(i) * ( RHOnP2_Chock(i) / RHOaP2_Chock(i) ) ^ (k);
        TnP2_Chock(i)           = TaP2_Chock(i) * ( RHOnP2_Chock(i) / RHOaP2_Chock(i) ) ^ (k - 1.0);
        UnP2_Chock(i)           = - ( k * R * TnP2_Chock(i) )^(0.5);
        MnP2_Chock(i)           = 1;
        MassFlowP2_Chock(i)     = RHOnP2_Chock(i) * An * UnP2_Chock(i);
        ThrustP2_Chock(i)       = 2 * PaP2_Chock(i) * An * ( 2.0 / (k+1.0) )^(1.0 / (k-1.0)) - Patm * An;
        
        if i == 1
            RocketMTotalP2_Chock(i) = RocketMTot02;
        else
            RocketMTotalP2_Chock(i) = RocketMTotalP2_Chock(i-1) + N*MassFlowP2_Chock(i) * (tP2_Chock(i) - tP2_Chock(i-1));
            
            RocketAccP2_Chock(i)    = ( N*ThrustP2_Chock(i) + drag(RocketVP2_Chock(i-1)) ) / RocketMTotalP2_Chock(i) - g;
            if i >= 2
                RocketVP2_Chock(i)  = RocketVP2_Chock(i-1) + RocketAccP2_Chock(i) * ( tP2_Chock(i) - tP2_Chock(i-1) );
                RocketYP2_Chock(i)  = RocketYP2_Chock(i-1) + RocketVP2_Chock(i)   * ( tP2_Chock(i) - tP2_Chock(i-1) );
            end
        end
    end
    RocketMTot03   = RocketMTotalP2_Chock(end);
    MassFlow03     = MassFlowP2_Chock(end);
    Pa03           = PaP2_Chock(end);
    Ta03           = TaP2_Chock(end);
    RHOa03         = RHOaP2_Chock(end);
    RocketThrust03 = ThrustP2_Chock(end);
    RocketAcc03    = RocketAccP2_Chock(end);
    RocketV03      = RocketVP2_Chock(end);
    RocketY03      = RocketYP2_Chock(end);
    
    MassAirExpelledP2_Chock = - N*trapz(tP2_Chock(1:tOutP2_Chock),MassFlowP2_Chock);

else
    
    Pa03                      = Pa02;
    MassFlow03                = 0.0001;
    RocketV03                 = RocketV02;
    RocketY03                 = RocketY02;
    RocketAcc03               = RocketAcc02;
    RocketMTot03              = RocketMTot02;
    MassAirExpelledP2_Chock = 0;
    tOutP2_Chock            = 1;
    MassFlowP2_Chock        = zeros(tOutP2_Chock,1);
    RocketMTotalP2_Chock    = zeros(tOutP2_Chock,1);
    RocketMTotalP2_Chock(1) = RocketMTot02;
    RocketAccP2_Chock       = zeros(tOutP2_Chock,1);
    RocketVP2_Chock         = zeros(tOutP2_Chock,1);
    RocketYP2_Chock         = zeros(tOutP2_Chock,1);
    ThrustP2_Chock          = zeros(tOutP2_Chock,1);
    RocketAccP2_Chock(1)    = RocketAccP1(end);
    RocketVP2_Chock(1)      = RocketVP1(end);
    RocketYP2_Chock(1)      = RocketYP1(end);
    ThrustP2_Chock(1)       = ThrustP1(end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Regime 2 --> Subsonic Flow

PratioTerm = @(P) (1-(Patm/P)^((k-1)/(k)));

eta = @(P) sqrt(2*k/(k-1)*(Patm/P)^(2/k)*PratioTerm(P));

pSqetaSq = @(P) P^(2)*2*k/(k-1)*(Patm/P)^(2/k)*PratioTerm(P);



% System of ODEs
% P       = x(1)
% Mdot    = x(2)

fP2_Subsonic = @(t,x) [-2*An^(2)*k^(2)/(Vb*abs(x(2))*(k-1))*(x(1)^(2)*(Patm/x(1))^((k+1)/(k))+Patm^(2)*(Patm/x(1))^(-4*(k-1)/(k))*PratioTerm(x(1)))*PratioTerm(x(1));
    ((x(1) > 1.00001*Patm | x(1) < 0.999999*Patm) & x(2) ~= 0)*(-abs(x(2))/(2*x(1))*((((k-1)/k)*(2+1/PratioTerm(x(1))*(Patm/x(1))^((k-1)/(k))))*(  -2*An^(2)*k^(2)/(Vb*abs(x(2))*(k-1))*(x(1)^(2)*(Patm/x(1))^((k+1)/(k))+Patm^(2)*(Patm/x(1))^(-4*(k-1)/(k))*PratioTerm(x(1)))*PratioTerm(x(1))  )-An^(2)*pSqetaSq(x(1))/(Vb*abs(x(2))))) + (x(2) == 0)*(1/(2*x(1))*(An^(2)*pSqetaSq(x(1))/Vb))];


% Setting options to solve ODEs system
t0P2_Subsonic   = 0;
X0P2_Subsonic   = [Pa03 MassFlow03];
tMaxP2_Subsonic = 10;
dtP2_Subsonic   = 1e-5;
tP2_Subsonic    = t0P2_Subsonic:dtP2_Subsonic:tMaxP2_Subsonic;
tspan = [0 10];
options=odeset('RelTol',1.e-6);

% Solving ODE
[tP2_Subsonic,solP2_Subsonic] = ode45(fP2_Subsonic,tspan,X0P2_Subsonic,options);

PaP2_Subsonicaux       = solP2_Subsonic(:,1);
MassFlowP2_Subsonicaux = solP2_Subsonic(:,2);

tOutP2_Subsonic = 1;
while PaP2_Subsonicaux(tOutP2_Subsonic) > 1.000001*Patm   % this is just a counter
    tOutP2_Subsonic = tOutP2_Subsonic + 1;
end

tP2_Subsonic_bis = 0:dtP2_Subsonic:tP2_Subsonic(tOutP2_Subsonic);
solP2_Subsonic = interp1(tP2_Subsonic,solP2_Subsonic,tP2_Subsonic_bis);
tP2_Subsonic = tP2_Subsonic_bis;
tOutP2_Subsonic = length(tP2_Subsonic);

PaP2_Subsonicaux       = solP2_Subsonic(:,1);
MassFlowP2_Subsonicaux = solP2_Subsonic(:,2);

PaP2_Subsonic           = PaP2_Subsonicaux(1:tOutP2_Subsonic);
MassFlowP2_Subsonic     = MassFlowP2_Subsonicaux(1:tOutP2_Subsonic);

TaP2_Subsonic           = zeros(tOutP2_Subsonic,1);
RHOaP2_Subsonic         = zeros(tOutP2_Subsonic,1);
UnP2_Subsonic           = zeros(tOutP2_Subsonic,1);
PnP2_Subsonic           =  ones(tOutP2_Subsonic,1) * Patm;
TnP2_Subsonic           = zeros(tOutP2_Subsonic,1);
RHOnP2_Subsonic         = zeros(tOutP2_Subsonic,1);
MnP2_Subsonic           = zeros(tOutP2_Subsonic,1);
ThrustP2_Subsonic       = zeros(tOutP2_Subsonic,1);
RocketMTotalP2_Subsonic = zeros(tOutP2_Subsonic,1);
RocketAccP2_Subsonic    = zeros(tOutP2_Subsonic,1);
RocketVP2_Subsonic      = zeros(tOutP2_Subsonic,1);
RocketYP2_Subsonic      = zeros(tOutP2_Subsonic,1);


RocketVP2_Subsonic(1)   = RocketV03;
RocketYP2_Subsonic(1)   = RocketY03;
RocketAccP2_Subsonic(1) = RocketAcc03;


for i = 1:tOutP2_Subsonic
    TaP2_Subsonic(i)           = 1/R*(An*PaP2_Subsonic(i)/(-MassFlowP2_Subsonic(i)))^(2)*eta(PaP2_Subsonic(i))^(2);
    RHOaP2_Subsonic(i)         = PaP2_Subsonic(i)/(R*TaP2_Subsonic(i));
    TnP2_Subsonic(i)           = TaP2_Subsonic(i)*(Patm/PaP2_Subsonic(i))^((k-1)/k);
    RHOnP2_Subsonic(i)         = Patm/(R*TnP2_Subsonic(i));
    UnP2_Subsonic(i)           = MassFlowP2_Subsonic(i)/(An*RHOnP2_Subsonic(i));
    MnP2_Subsonic(i)           = -UnP2_Subsonic(i)/((k*R*TnP2_Subsonic(i))^0.5);
    ThrustP2_Subsonic(i)       = MassFlowP2_Subsonic(i)*UnP2_Subsonic(i);
    
    if i == 1
        RocketMTotalP2_Subsonic(i) = RocketMTot03;
    else
        RocketMTotalP2_Subsonic(i) = RocketMTotalP2_Subsonic(i-1) + N*MassFlowP2_Subsonic(i) * (tP2_Subsonic(i) - tP2_Subsonic(i-1));
    end
    
    if i >= 2
        RocketAccP2_Subsonic(i) = ( N*ThrustP2_Subsonic(i) + drag(RocketVP2_Subsonic(i-1)) ) / RocketMTotalP2_Subsonic(i) - g;  %%%%% CHECK %%%%%
        RocketVP2_Subsonic(i) = RocketVP2_Subsonic(i-1) + RocketAccP2_Subsonic(i) * ( tP2_Subsonic(i) - tP2_Subsonic(i-1) );
        RocketYP2_Subsonic(i) = RocketYP2_Subsonic(i-1) + RocketVP2_Subsonic(i)   * ( tP2_Subsonic(i) - tP2_Subsonic(i-1) );
    end
    
end
RocketMTot04   = RocketMTotalP2_Subsonic(end);
MassFlow04     = MassFlowP2_Subsonic(end);
Pa04           = PaP2_Subsonic(end);
Ta04           = TaP2_Subsonic(end);
RHOa04         = RHOaP2_Subsonic(end);
RocketThrust04 = ThrustP2_Subsonic(end);
RocketAcc04    = RocketAccP2_Subsonic(end);
RocketV04      = RocketVP2_Subsonic(end);
RocketY04      = RocketYP2_Subsonic(end);

MassAirExpelledP2_Subsonic = -N*trapz(tP2_Subsonic(1:tOutP2_Subsonic), MassFlowP2_Subsonic);
TotalMassAirExpelled      =  MassAirExpelledP2_Chock + MassAirExpelledP2_Subsonic;
MassAirLeftInTank         =  Vb*RHOa04;

%% Phase 3: Free flight - Maximum apogee is reached while thrust phase has ended

t0P3      = 0;
tMaxP3    = 20;
dtP3      = 1e-5;
tP3aux    = t0P3:dtP3:tMaxP3;
nTSP3aux  = length(tP3aux);

RocketMTotalP3aux = zeros(nTSP3aux,1);
RocketAccP3aux    = zeros(nTSP3aux,1);
RocketVP3aux      = zeros(nTSP3aux,1);
RocketYP3aux      = zeros(nTSP3aux,1);

RocketVP3aux(1)   = RocketV04;
RocketYP3aux(1)   = RocketY04;

for i = 1:nTSP3aux
    
    RocketMTotalP3aux(i) = RocketMTot04;
    if i == 1
        RocketAccP3aux(i) = ( drag(RocketVP3aux(i))   ) / RocketMTotalP3aux(i) - g;
    else
        RocketAccP3aux(i) = ( drag(RocketVP3aux(i-1)) ) / RocketMTotalP3aux(i) - g;
    end
    
    if i >= 2
        RocketVP3aux(i) = RocketVP3aux(i-1) + RocketAccP3aux(i) * ( tP3aux(i) - tP3aux(i-1) );
        RocketYP3aux(i) = RocketYP3aux(i-1) + RocketVP3aux(i)   * ( tP3aux(i) - tP3aux(i-1) );
    end
    
end

[aux, nTSP3] = max(RocketYP3aux);
nTSP3 = nTSP3 + 1;


tP3            = tP3aux(1:nTSP3);
RocketMTotalP3 = RocketMTotalP3aux(1:nTSP3);
RocketAccP3    = RocketAccP3aux(1:nTSP3);
RocketVP3      = RocketVP3aux(1:nTSP3);
RocketYP3      = RocketYP3aux(1:nTSP3);


%%  Combining all phases to have complete flight history for each quantity of interest

nTSTotal  = tOutP0 + tOutP1 + tOutP2_Chock + tOutP2_Subsonic + nTSP3;


tTotal    = [ tP0(1:tOutP0), tP0(tOutP0) + tP1(1:tOutP1),...
    tP0(tOutP0) + tP1(tOutP1) + tP2_Chock(1:tOutP2_Chock),...
    tP0(tOutP0) + tP1(tOutP1) + tP2_Chock(tOutP2_Chock) + tP2_Subsonic(1:tOutP2_Subsonic),...
    tP0(tOutP0) + tP1(tOutP1) + tP2_Chock(tOutP2_Chock) + tP2_Subsonic(tOutP2_Subsonic) + tP3 ];


%
%   Initializing remaining total variables
RocketThrustTotal  = zeros(nTSTotal,1);
RocketMassTotal    = zeros(nTSTotal,1);
RocketAccTotal     = zeros(nTSTotal,1);
RocketVTotal       = zeros(nTSTotal,1);
RocketYTotal       = zeros(nTSTotal,1);
WaterMassFlowTotal = zeros(nTSTotal,1);
AirMassFlowTotal   = zeros(nTSTotal,1);
PaTotal            = zeros(nTSTotal,1);
TaTotal            = zeros(nTSTotal,1);
RHOaTotal          = zeros(nTSTotal,1);
PanTotal           = zeros(nTSTotal,1);
TanTotal           = zeros(nTSTotal,1);
RHOanTotal         = zeros(nTSTotal,1);
PwnTotal           = zeros(nTSTotal,1);
TwnTotal           = zeros(nTSTotal,1);
RHOwnTotal         = zeros(nTSTotal,1);
ManTotal           = zeros(nTSTotal,1);
UanTotal           = zeros(nTSTotal,1);
UwnTotal           = zeros(nTSTotal,1);

%   Air and water mass flows

WaterMassFlowTotal(tOutP0+1:tOutP0+tOutP1)                                                           = MassFlowP1(1:tOutP1);
AirMassFlowTotal(tOutP0+tOutP1+1:tOutP0+tOutP1+tOutP2_Chock)                                   = MassFlowP2_Chock(1:tOutP2_Chock);
AirMassFlowTotal(tOutP0+tOutP1+tOutP2_Chock+1:tOutP0+tOutP1+tOutP2_Chock+tOutP2_Subsonic)     = MassFlowP2_Subsonic(1:tOutP2_Subsonic);

%   Rocket mass

RocketMassTotal(1:tOutP0)                                                                                    = M0;
RocketMassTotal(tOutP0+1:tOutP0+tOutP1)                                                              = RocketMTotalP1(1:tOutP1);
RocketMassTotal(tOutP0+tOutP1+1:tOutP0+tOutP1+tOutP2_Chock)                                    = RocketMTotalP2_Chock(1:tOutP2_Chock);
RocketMassTotal(tOutP0+tOutP1+tOutP2_Chock+1:tOutP0+tOutP1+tOutP2_Chock+tOutP2_Subsonic)      = RocketMTotalP2_Subsonic(1:tOutP2_Subsonic);
RocketMassTotal(tOutP0+tOutP1+tOutP2_Chock+tOutP2_Subsonic+1:nTSTotal)                                  = RocketMTot04;

%   Rocket Thrust

RocketThrustTotal(1:tOutP0)                                                                                   = ThrustP0(1:tOutP0);
RocketThrustTotal(tOutP0+1:tOutP0+tOutP1)                                                             = ThrustP1(1:tOutP1);
RocketThrustTotal(tOutP0+tOutP1+1:tOutP0+tOutP1+tOutP2_Chock)                                   = ThrustP2_Chock(1:tOutP2_Chock);
RocketThrustTotal(tOutP0+tOutP1+tOutP2_Chock+1:tOutP0+tOutP1+tOutP2_Chock+tOutP2_Subsonic)     = ThrustP2_Subsonic(1:tOutP2_Subsonic);
RocketThrustTotal(tOutP0+tOutP1+tOutP2_Chock+tOutP2_Subsonic+1:nTSTotal)                                 = 0;

%   Rocket Acceleration

RocketAccTotal(1:tOutP0)                                                                                     = RocketAccP0(1:tOutP0);
RocketAccTotal(tOutP0+1:tOutP0+tOutP1)                                                               = RocketAccP1(1:tOutP1);
RocketAccTotal(tOutP0+tOutP1+1:tOutP0+tOutP1+tOutP2_Chock)                                     = RocketAccP2_Chock(1:tOutP2_Chock);
RocketAccTotal(tOutP0+tOutP1+tOutP2_Chock+1:tOutP0+tOutP1+tOutP2_Chock+tOutP2_Subsonic)       = RocketAccP2_Subsonic(1:tOutP2_Subsonic);
RocketAccTotal(tOutP0+tOutP1+tOutP2_Chock+tOutP2_Subsonic+1:nTSTotal)                                   = RocketAccP3(1:nTSP3);


%   Rocket Velocity

RocketVTotal(1:tOutP0)                                                                                       = RocketVP0(1:tOutP0);
RocketVTotal(tOutP0+1:tOutP0+tOutP1)                                                                 = RocketVP1(1:tOutP1);
RocketVTotal(tOutP0+tOutP1+1:tOutP0+tOutP1+tOutP2_Chock)                                       = RocketVP2_Chock(1:tOutP2_Chock);
RocketVTotal(tOutP0+tOutP1+tOutP2_Chock+1:tOutP0+tOutP1+tOutP2_Chock+tOutP2_Subsonic)         = RocketVP2_Subsonic(1:tOutP2_Subsonic);
RocketVTotal(tOutP0+tOutP1+tOutP2_Chock+tOutP2_Subsonic+1:nTSTotal)                                     = RocketVP3(1:nTSP3);

%   Rocket Height

RocketYTotal(1:tOutP0)                                                                                       = RocketYP0(1:tOutP0);
RocketYTotal(tOutP0+1:tOutP0+tOutP1)                                                                 = RocketYP1(1:tOutP1);
RocketYTotal(tOutP0+tOutP1+1:tOutP0+tOutP1+tOutP2_Chock)                                       = RocketYP2_Chock(1:tOutP2_Chock);
RocketYTotal(tOutP0+tOutP1+tOutP2_Chock+1:tOutP0+tOutP1+tOutP2_Chock+tOutP2_Subsonic)         = RocketYP2_Subsonic(1:tOutP2_Subsonic);
RocketYTotal(tOutP0+tOutP1+tOutP2_Chock+tOutP2_Subsonic+1:nTSTotal)                                     = RocketYP3(1:nTSP3);


%Printing maximum altitude (major quantity of interest
fprintf('Max altitude = %f\n',RocketYP3(end))






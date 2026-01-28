% ===============================================================
% Boeing 747 Landing Performance

 function landingPerformanceGUI

 %---------- aircraft constants---------

S = 525;                                     % wing area [m^2]

T_max = 4 * 251e3;                           % total rated thrust [N]

%---------- default values---------

default.mass = 260360;                       % kg (max landing weight)
default.alt = 112;                           % m (Frankfurt)
default.CLmax = 1.78;                        % landing-flap
default.gammaDeg = 3;                        % deg (glide slope)
default.wind =-1.6777;                       % m/s (+ headwind)
default.slopeDeg =-0.30;                     % deg (down-slope)
default.mu = 0.2;                            % dry runway
default.revOn = false;
default.revFrac = 0.50;

 %---------- GUI layout---------

fig = figure('Name','747 Landing Performance','NumberTitle','off',...
'Position',[100 100 600 500]);
newLbl = @(t,y) uicontrol(fig,'Style','text','Units','normalized',...
'Position',[0.05 y 0.40 0.05],'String',t,...
'HorizontalAlignment','left','FontWeight','bold');
newEdt = @(v,y) uicontrol(fig,'Style','edit','Units','normalized',...
'Position',[0.45 y 0.20 0.05],'String',num2str(v));

newLbl('Aircraft Mass (kg):' ,0.88); massEdt = newEdt(default.mass ,0.88);
newLbl('Landing Altitude (m):' ,0.80); altEdt = newEdt(default.alt ,0.80);
newLbl('Max Lift Coeff C_L_{max}:' ,0.72); clEdt = newEdt(default.CLmax,0.72);
newLbl('Approach Angle (deg):' ,0.64); angEdt = newEdt(default.gammaDeg,0.64);
newLbl('Wind (m/s, + headwind):' ,0.56); wndEdt = newEdt(default.wind ,0.56);
newLbl('Runway Slope (deg):' ,0.48); slpEdt = newEdt(default.slopeDeg,0.48);
newLbl('Braking Coefficient mu:' ,0.40); muEdt = newEdt(default.mu ,0.40);

revChk = uicontrol(fig,'Style','checkbox','Units','normalized',...
'Position',[0.05 0.32 0.40 0.05],'String','Use Reverse Thrust',...
'Value',default.revOn,'FontWeight','bold',...
'HorizontalAlignment','left','Callback',@tglRev);
newLbl('Reverse Thrust Fraction:' ,0.25);
rvEdt = newEdt(default.revFrac,0.25);
if ~default.revOn, set(rvEdt,'Enable','off'); 
end

uicontrol(fig,'Style','pushbutton','Units','normalized',...
'Position',[0.75 0.80 0.18 0.08],'String','Compute',...
'FontSize',10,'FontWeight','bold','Callback',@computeLanding);

 outTxt = uicontrol(fig,'Style','text','Units','normalized',...
'Position',[0.05 0.02 0.90 0.20],'HorizontalAlignment','left',...
'String','Results appearhere.');

  %----------callbacks---------

function tglRev(src,~), set(rvEdt,'Enable',tern(src.Value,'on','off')); 
end

 function s= tern(c,a,b); if c,s = a; else, s =b; end, end

%--------------------------------------------------------------
% main COMPUTE callback
%--------------------------------------------------------------

 function computeLanding(~,~)

%-------user parameters------
m = str2double(get(massEdt,'String'));
alt = str2double(get(altEdt ,'String'));
CLmax = str2double(get(clEdt ,'String'));
gammaD = str2double(get(angEdt ,'String'));
wind = str2double(get(wndEdt ,'String'));
slopeD = str2double(get(slpEdt ,'String'));
mu = str2double(get(muEdt ,'String'));
revOn = get(revChk,'Value');
revF = str2double(get(rvEdt ,'String'));

rho = ISA_rho(alt);
W = m*9.81;
g = 9.81;

 %-------key speeds----------

Vs = sqrt(2*W/(rho*S*CLmax));
Vapp = 1.3*Vs;
VtdA = 1.15*Vs;

 %-------approach + flare----

th = deg2rad(gammaD);
Hobs = 15.24;
s_app = Hobs/tan(th);

% simple flare model

Vbar = 0.5*(Vapp + VtdA);
t_flr = min(8, (Vapp*sin(th))/(0.2*g));
s_flare = Vbar*t_flr;
s_air = s_app + s_flare;

 %-------ground roll---------

VtdG = max(0, VtdA- wind);
t_free = 2;
s_free = VtdG*t_free;
T_rev = revF*0.5*T_max*(revOn && revF>0);
V = VtdG;s_g = 0; Ebr = 0;
dt = 0.1;
while V >0
Vair = max(0,V + wind);
L = 0.20*W*(Vair/VtdA)^2;                       % residuallift
CLi = L/(0.5*rho*S*max(Vair^2,1e-6));
CDi = 0.119+ 0.0435*CLi^2;
Fd = 0.5*rho*S*CDi*Vair^2;
Fb = mu*(W- L);
a_tot = (Fb+ Fd + T_rev +W*sin(deg2rad(slopeD)))/m;
a_tot = max(a_tot,0.1);
Vn = V-a_tot*dt; 
if Vn <0, dt = V/a_tot; Vn= 0; 
end
s_step = V*dt- 0.5*a_tot*dt^2;
s_g = s_g+ s_step; Ebr = Ebr+ Fb*s_step; V =Vn;
end
s_ground = s_free + s_g;
s_total = s_air + s_ground;

set(outTxt,'String',{ ...
sprintf('Stall SpeedV_s : %6.1f m/s (%.1f kt)',Vs ,Vs*1.944)
sprintf('Approach SpeedV_app : %6.1f m/s(%.1f kt)',Vapp,Vapp*1.944)
sprintf('Touchdown SpeedV_td : %6.1f m/s(%.1f kt)',VtdA,VtdA*1.944)
sprintf('Approach Dist(50 ft->flare): %6.0fm',s_app)
sprintf('Flare Distance: %6.0f m',s_flare)
sprintf('Ground-Roll Distance : %6.0f m',s_ground)
sprintf('TOTAL LandingDistance : %6.0f m',s_total )
sprintf('Brake EnergyAbsorbed : %6.1f MJ',Ebr/1e6)
});

%---------sweep plots---------------------------

makePlots(mu,revOn,revF, m, wind, alt, CLmax,...
t_free,VtdA,VtdG,slopeD,gammaD);
 end

%----------------------------------------------------------------------

%--------------------------------------------------------------
% PLOTTING+ sensitivity analysis
%--------------------------------------------------------------

 function makePlots(mu0, revOn, revF, m,wind, alt, CLmax,...
t_free, VtdA_base, VtdG_base, slopeD, gammaDeg)      %#ok<AGROW>

%---helperthat reproduces physicsfrom computeLanding---------

function s_tot= landDist(mass, mu,revFrac, slopeDeg)
rho = ISA_rho(alt); g = 9.81; W= mass*g;
Vs = sqrt(2*W/(rho*S*CLmax) );
Vapp = 1.3*Vs; VtdA = 1.15*Vs;
VtdG = max(0,VtdA- wind);

% airborne segment

th = deg2rad(gammaDeg);
Hobs = 15.24;
s_app = Hobs/ tan(th);
t_flr = min(8, (Vapp*sin(th))/(0.2*g));
s_flr = 0.5*(Vapp + VtdA)*t_flr;
s_air = s_app+ s_flr;

% ground roll

T_rev = revFrac * 0.5 * T_max* (revOn && revFrac> 0);
V = VtdG;s_g = 0; dt = 0.1;
while V >0
Vair = max(0, V + wind);
L = 0.20*W*(Vair/VtdA)^2;
CLi = L/ (0.5*rho*S*max(Vair^2,1e-6));
CDi = 0.119 + 0.0435*CLi^2;
Fd = 0.5*rho*S*CDi*Vair^2;
Fb = mu*(W- L);
at = (Fb+ Fd + T_rev + W*sin(deg2rad(slopeDeg)))/mass;
at = max(at,0.1);
Vn = V- at*dt; 
if Vn < 0, dt = V/at; Vn= 0; 
end
s_g = s_g+ V*dt- 0.5*at*dt^2;
V = Vn;
end
s_tot = s_air+ t_free*VtdG + s_g;
end

%---------------------------------------------------------------

 %--------sweep 1: runway frictioncoefficient-----------------

muVec = linspace(0.20, 0.80, 13);
dMu = arrayfun(@(mu) landDist(m, mu,revF, slopeD), muVec);

%--------sweep 2: runway slope(deg)--------------------------

slopeVec = linspace(-1.0, 1.5, 26);                                     % downhill(-) to uphill(+)
dSlope = arrayfun(@(sl) landDist(m, mu0,revF, sl), slopeVec);

 %-----------------------plotting------------------------------

figure('Name','747 Landing-Distance Sensitivities','Color','white');

subplot(1,2,1)
plot(muVec, dMu/1000, 'LineWidth',1.4);
grid on
xlabel('Runway friction coefficient \\mu');
ylabel('Landing distance [km]');
title('Sensitivity to\\mu');
subplot(1,2,2)
plot(slopeVec, dSlope/1000, 'LineWidth',1.4);
grid on
xlabel('Runway slope[deg] (+ uphill)');
ylabel('Landing distance [km]');
title('Sensitivity torunway slope');
 end

%--------------------------------------------------------------
 end           % GUIwrapper
%----------------------------------------------------------------

%----------ISA density Feedback---------

function rho= ISA_rho(h)
rho0 = 1.225; T0 = 288.15; P0= 101325; L = 0.0065; R = 287.05; g= 9.81;
h = max(h,0);
if h <= 11000
T = T0-L*h; P = P0*(T/T0)^(g/(R*L)); rho = P/(R*T);
else
T11 = T0- L*11000; P11 = P0*(T11/T0)^(g/(R*L));
P = P11*exp(-g*(h-11000)/(R*T11)); rho =P/(R*T11);
end
end

function cruisePerformanceGUI

%======================================================================
% Boeing 747-400- Cruise Performance
%---------------------------------------------------------------------


%----------------- AIRCRAFT / ENGINE CONSTANTS----------------

S = 541.2;                                                            % wing planform area [m^2]
SFC_cruise = 0.0575119944;                                            % kg/(N h)
g = 9.80665;                                                          % gravity [m/s^2]

%----------------- DEFAULT GUI VALUES-------------------------

def.mass0 = 387282;                                                   % kg (= weight after climb in mission)
def.altCru = 10668;                                                   % m (FL350)
def.mach = 0.85;                                                      % cruise Mach number
def.timeHr = 5.0;                                                     % cruise duration [h]


%----------------- BUILD GUI----------------------------------

f = figure('Name','747 Cruise Perf','NumberTitle','off',...
'Position',[120 120 650 430],'Resize','off');

makeLabel('Mass at start of cruise [kg]', 350); mEdit = makeEdit(def.mass0,350);
makeLabel('Cruise altitude [m]', 310); aEdit = makeEdit(def.altCru,310);
makeLabel('Cruise Mach number ', 270); MEdit = makeEdit(def.mach,270);
makeLabel('Cruise duration [h]', 230); tEdit = makeEdit(def.timeHr,230);


uicontrol(f,'Style','push','Position',[440 300 170 45],'FontWeight','bold',...
'String','Compute Cruise','Callback',@runCruise);


 outTxt = uicontrol(f,'Style','text','Position',[25 25 600 200],...
'Horizontal','left','FontSize',10,'String',' ');

%==============================================================
% FEEDBACK- compute cruise
%==============================================================

function runCruise(~,~)
try
m = str2double(get(mEdit,'String'));                                         % kg
h = str2double(get(aEdit,'String'));                                         % m
M = str2double(get(MEdit,'String'));                                         % Mach
tHr = str2double(get(tEdit,'String'));                                       % h
if any(isnan([m,h,M,tHr])) || m<=0 || h<=0 || M<=0 || tHr<=0
error('All inputs must be positive numbers.');
end

 %---atmosphere at cruise altitude

[rho, a] =ISA_atmosphere(h);
V = M *a;                                                                    % TAS [m/s]
dt = 60;                                                                     % integration timestep [s] (1 min)
N = ceil(tHr*3600 / dt);                                                     % numberof steps
fuelTotal =0;

%---integration loop--------------------------------

for k = 1:N
W = m* g;                                                                    % weight [N]
CL = W/ (0.5*rho*V^2*S);                                                     % lift coefficient
CD = dragCoefficient(CL,'cruise');
D = 0.5*rho*V^2*S*CD;                                                        % drag (= thrust required )
if CL >1.5                                                                   % simple sanity:747 cant sustain CL>1.5
error(['At %.0f mand Mach %.2f the required lift coefficient(%.2f) ' ...
'is too high-choose lower altitude, shorter time, orlower mass.'], ...
h,M,CL);
end

                                                                             % Fuel flow (kg) this minute =T_req * SFC * dt/3600
fuel = D* SFC_cruise * (dt/3600);
                                                                             % Preventnegative mass if entered huge time

if fuel >m*0.01                                                              % more than1 % of remaining mass per minute? unlikely
warning('Fuel burn very high; stopping integration early.');
break;
end
fuelTotal =fuelTotal + fuel;
m = m- fuel;                                                                 % update mass
end
tAct = N*dt/3600;                                                            % actual integratedtime [h]

%--------------summary text-------------------------

txt = sprintf('CRUISE SUMMARY-Boeing 747-400\n');
txt = [txt sprintf('Start mass :%9.0f kg\n', str2double(get(mEdit,'String')))];
txt = [txt sprintf('Cruise altitude: %7.0f m (%.0f ft)\n', h, h/0.3048)];
txt = [txt sprintf('Mach number: %.2f (TAS = %.0fm/s)\n', M,V)];
txt = [txt sprintf('Requested time: %.2f h\n', tHr)];
txt = [txt sprintf('Integrated time: %.2f h\n', tAct)];
txt = [txt sprintf('--------------------------------------------\n')];
txt = [txt sprintf('Fuel burned: %7.0f kg\n', fuelTotal)];
txt = [txt sprintf('Final mass :%9.0f kg\n', m)];
set(outTxt,'String',txt); drawnow;
catch ME
errordlg(ME.message,'Cruise Error'); rethrow(ME);
end
end

%==============================================================
% GUI HELPERS
%==============================================================

function makeLabel(str,y)
     uicontrol(f,'Style','text','Position',[30 y 200 25],'String',str,'Horizontal','left');
end
function e= makeEdit(val,y)
e = uicontrol(f,'Style','edit','Position',[250 y 120 25],...
'String',num2str(val));
end
end 
%=========================== END GUI FUNCTION ==================

%======================================================================
% LOCAL FUNCTIONS
%======================================================================

function [rho,a] = ISA_atmosphere(h)
% 1976 ISA (troposphere + isothermal lower stratosphere)
T0=288.15; p0=101325; L=0.0065; R=287.058; g=9.80665;
if h<11000
T=T0-L*h; p=p0*(T/T0)^(g/(L*R));
else
T=216.65; p11=p0*(T/T0)^(g/(L*R));
p=p11*exp(-g*(h-11000)/(R*T));
end
rho=p/(R*T); a=sqrt(1.4*R*T);
end
 function CD = dragCoefficient(CL,phase)
if strcmpi(phase,'cruise'), CD0 = 0.025; else, CD0 = 0.045; end
K = 0.0435; CD = CD0 + K*CL^2;
 end
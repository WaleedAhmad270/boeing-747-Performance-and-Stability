function ClimbPerformanceGui
%======================================================================
% Boeing 747-400 - Climb Performance (GUI)
%======================================================================

%----------------- CONSTANTS ----------------
S = 541.2;                                      % wing area [m^2]
T_max_SL = 4 * 208390.66251372406;              % max thrust SL [N]
T_climb  = 0.8 * T_max_SL;                      % climb thrust [N]
SFC_climb = 0.0508770819;                       % (your assumed SFC model)
g = 9.80665;

%----------------- DEFAULTS ------------------
d.mass0 = 396505.6289676431;                    % kg
d.alt0  = 173;                                  % m
d.altC  = 10668;                                % m (FL350)

%----------------- GUI -----------------------
f = figure('Name','747 ClimbPerf','NumberTitle','off', ...
    'Position',[200 200 520 380], 'Resize','off');

makeLabel('Initial mass [kg]', 330);  mEdit = makeEdit(d.mass0,330);
makeLabel('Cruise altitude [m]',290); cEdit = makeEdit(d.altC,290);
makeLabel('Airport elevation [m]',250); aEdit = makeEdit(d.alt0,250);

uicontrol(f,'Style','pushbutton','Position',[330 285 170 45], ...
    'FontWeight','bold','String','Compute Climb','Callback',@runClimb);

outTxt = uicontrol(f,'Style','edit','Position',[25 25 470 205], ...
    'HorizontalAlignment','left','Max',10,'Min',1, ...   % multiline
    'FontSize',10,'String','');

%==============================================================
% CALLBACK
%==============================================================
function runClimb(~,~)
    try
        m  = str2double(get(mEdit,'String'));   % current mass [kg]
        h0 = str2double(get(aEdit,'String'));   % airport elevation [m]
        hC = str2double(get(cEdit,'String'));   % cruise altitude [m]

        if isnan(m) || isnan(h0) || isnan(hC) || m<=0 || hC<=h0
            error('Check numeric inputs: mass > 0, cruise altitude > airport elevation.');
        end

        %---- segment definitions ---
        segAlt   = [3048; 7620];                % 10k ft, 25k ft (m)
        segSpeed = [131.7; 228.2];              % 250 kt, 340 kt (m/s)
        segLabel = {'0-10 kft (250 kt)'; '10-25 kft (340 kt)'};

        h = h0;
        fuelTot = 0;
        segInfo = {}; % {label,fuel,time}

        %---- low-altitude segments ---
        for s = 1:numel(segAlt)
            if h < segAlt(s) && hC > h
                [fuel,dt,m,h] = climbSeg(m, h, min(hC,segAlt(s)), segSpeed(s));
                fuelTot = fuelTot + fuel;
                segInfo(end+1,:) = {segLabel{s}, fuel, dt}; %#ok<AGROW>
            end
        end

        %---- high-altitude (Mach) segment ---
        if hC > h
            [~,a] = ISA_atmosphere((h+hC)/2);
            V = 0.84 * a; % Mach 0.84
            [fuel,dt,m,h] = climbSeg(m, h, hC, V);
            fuelTot = fuelTot + fuel;
            segInfo(end+1,:) = {'25 kft - Cruise (M 0.84)', fuel, dt}; %#ok<AGROW>
        end

        %---- summary text ---
        txt = sprintf('CLIMB SUMMARY - Boeing 747-400\n');
        txt = [txt sprintf('Initial mass      : %9.0f kg\n', m + fuelTot)];
        txt = [txt sprintf('Airport elevation : %9.0f m\n', h0)];
        txt = [txt sprintf('Cruise altitude   : %9.0f m\n\n', hC)];

        for k = 1:size(segInfo,1)
            txt = [txt sprintf('%-24s  Fuel %7.1f kg   Time %7.1f s\n', ...
                segInfo{k,1}, segInfo{k,2}, segInfo{k,3})]; %#ok<AGROW>
        end

        txt = [txt sprintf('\n--------------------------------------------\n')];
        txt = [txt sprintf('Total climb fuel  : %7.1f kg\n', fuelTot)];
        txt = [txt sprintf('Mass after climb  : %9.1f kg\n', m)];

        set(outTxt,'String',txt);

    catch ME
        errordlg(ME.message,'Climb Error');
        rethrow(ME);
    end
end

%==============================================================
% ONE CLIMB SEGMENT
%==============================================================
function [fuel,timeNew,mNew,h2out] = climbSeg(mCur,h1,h2in,V)
    % Inputs: mass [kg], start alt h1 [m], end alt h2in [m], TAS V [m/s]
    % Returns: fuel [kg], time [s], new mass, new altitude

    [rho,~] = ISA_atmosphere((h1+h2in)/2);
    W = mCur * g;

    CL = W / (0.5*rho*V^2*S);
    CD = dragCoefficient(CL,'climb');
    D  = 0.5*rho*V^2*S*CD;

    excess = T_climb - D;
    if excess <= 0
        error('Insufficient climb thrust near %.0f m: Drag %.0f N > Thrust %.0f N.', h1, D, T_climb);
    end

    gamma = asin(excess / W);
    ROC = V * sin(gamma);                   % rate of climb [m/s]
    timeNew = (h2in - h1) / ROC;            % time [s]

    fuel = SFC_climb * T_climb * (timeNew/3600); % your model (kg)
    mNew = mCur - fuel;
    h2out = h2in;
end

%==============================================================
% SMALL GUI HELPERS
%==============================================================
function makeLabel(str,y)
    uicontrol(f,'Style','text','Position',[25 y 200 20], ...
        'String',str,'HorizontalAlignment','left');
end

function e = makeEdit(val,y)
    e = uicontrol(f,'Style','edit','Position',[220 y 100 25], ...
        'String',num2str(val));
end

end % end main GUI function

%======================================================================
% LOCAL FUNCTIONS
%======================================================================
function [rho,a] = ISA_atmosphere(h)
% 1976 ISA (troposphere + simple stratosphere)
T0=288.15; p0=101325; L=0.0065; R=287.058; g=9.80665;
if h < 11000
    T = T0 - L*h;
    p = p0*(T/T0)^(g/(L*R));
else
    T = 216.65;
    p11 = p0*(T/T0)^(g/(L*R));
    p = p11*exp(-g*(h-11000)/(R*T));
end
rho = p/(R*T);
a   = sqrt(1.4*R*T);
end

function CD = dragCoefficient(CL,phase)
if strcmpi(phase,'climb')
    CD0 = 0.045;
else
    CD0 = 0.025;
end
K = 0.0435;
CD = CD0 + K*CL^2;
end

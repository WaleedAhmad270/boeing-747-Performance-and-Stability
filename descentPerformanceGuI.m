function descentPerformanceGUI
%======================================================================
% Boeing 747-400 - Descent 
%======================================================================

%----------------- CONSTANTS (aircraft / engine) ----------------------
S = 541.2;                       % wing area [m^2]
SFC_idle = 0.057511;             % kg/(N h)  
g = 9.80665;                     % gravity [m/s^2]
TmaxSL_each = 208390.66251372406; % thrust per engine at SL [N]
T_idle = 0.061 * 4 * TmaxSL_each; % idle thrust total [N] 

%----------------- DEFAULT INPUTS -------------------------------------
def.mass0 = 299240;              % kg
def.alt1  = 11887.2;             % m (FL390)
def.altAPT = 112;                % m
def.vs = 8.89;                   % m/s (NOTE: this is ~1750 ft/min)

%----------------- BUILD GUI ------------------------------------------
f = figure('Name','747 Descent Perf','NumberTitle','off', ...
    'Position',[140 140 640 420],'Resize','off');

makeLabel(f,'Mass at start of descent [kg]', 350);
mEdit  = makeEdit(f,def.mass0,350);

makeLabel(f,'Start altitude [m]', 310);
h1Edit = makeEdit(f,def.alt1,310);

makeLabel(f,'Airport elevation [m]', 270);
h2Edit = makeEdit(f,def.altAPT,270);

makeLabel(f,'Target descent rate (positive down) [m/s]', 230);
vsEdit = makeEdit(f,def.vs,230);

uicontrol(f,'Style','push','Position',[435 300 165 45],'FontWeight','bold', ...
    'String','Compute Descent','Callback',@runDescent);

outTxt = uicontrol(f,'Style','text','Position',[25 25 590 200], ...
    'Horizontal','left','FontSize',10,'String',' ');

%======================================================================
% CALLBACK
%======================================================================

function runDescent(~,~)
    try
        m0   = str2double(get(mEdit,'String'));
        h1   = str2double(get(h1Edit,'String'));
        hAPT = str2double(get(h2Edit,'String'));
        vz   = str2double(get(vsEdit,'String')); % positive downward

        if any(isnan([m0,h1,hAPT,vz])) || m0<=0 || h1<=hAPT || vz<=0
            error('Inputs must be numeric: mass>0, start alt > airport alt, descent rate > 0.');
        end

        fuelTot = 0;
        m = m0;
        segInfo = {}; % {label, fuel, time}

        % Segment 1: start -> 10,000 ft at M0.84

        h10 = 3048; % 10,000 ft in m
        if h1 > h10
            [fuel, dt] = descentSeg(m, h1, max(h10,hAPT), 0.84, vz);
            segInfo(end+1,:) = {'FL390-10k ft (M0.84)', fuel, dt};
            fuelTot = fuelTot + fuel;
            m = m - fuel;
        end

        % Segment 2: 10,000 ft -> airport at 250 kt (TAS)
        if hAPT < h10
            V250 = 250 * 0.514444; % 250 knots -> m/s
            [fuel, dt] = descentSeg(m, max(h10,hAPT), hAPT, V250, vz);
            segInfo(end+1,:) = {'10k ft-APT (250 kt)', fuel, dt};
            fuelTot = fuelTot + fuel;
            m = m - fuel;
        end

        if isempty(segInfo)
            
            segInfo = {'No descent segment computed', 0, 0};
        end

        tMin = sum(cell2mat(segInfo(:,3))) / 60;

        % Summary

        txt = sprintf('DESCENT SUMMARY - Boeing 747-400\n');
        txt = [txt sprintf('Start mass        : %9.0f kg\n', m0)];
        txt = [txt sprintf('Start altitude    : %7.0f m (%.0f ft)\n', h1, h1/0.3048)];
        txt = [txt sprintf('Airport elevation : %7.0f m (%.0f ft)\n', hAPT, hAPT/0.3048)];
        txt = [txt sprintf('Target VS         : %4.1f m/s (%.0f ft/min)\n', vz, vz*196.85)];

        for k = 1:size(segInfo,1)
            txt = [txt sprintf(' %-23s : Fuel %6.0f kg | Time %5.0f s\n', ...
                segInfo{k,1}, segInfo{k,2}, segInfo{k,3})];
        end

        txt = [txt sprintf('--------------------------------------------\n')];
        txt = [txt sprintf('Total descent fuel: %7.0f kg\n', fuelTot)];
        txt = [txt sprintf('Landing mass      : %9.0f kg\n', m)];
        txt = [txt sprintf('Total descent time: %5.1f min\n', tMin)];

        set(outTxt,'String',txt); drawnow;

    catch ME
        errordlg(ME.message,'Descent Error');
        rethrow(ME);
    end
end

%======================================================================
% ONE DESCENT SEGMENT (idle thrust)
%======================================================================

function [fuel,timeSec] = descentSeg(mStart, hTop, hBot, MachOrTas, vz)
    % MachOrTas: if <= 1.5 treated as Mach; else treated as TAS [m/s]

    hMid = 0.5*(hTop + hBot);
    [rhoMid, aMid] = ISA_atmosphere(hMid);

    if MachOrTas <= 1.5
        V = MachOrTas * aMid;   % Mach
    else
        V = MachOrTas;          % TAS
    end

    W  = mStart * g;
    CL = W / (0.5*rhoMid*V^2*S);
    CD = dragCoefficient(CL,'cruise');
    D  = 0.5*rhoMid*V^2*S*CD; %#ok<NASGU> 

    dh = hTop - hBot;           % altitude loss [m]
    timeSec = dh / vz;          % seconds
    fuel = T_idle * SFC_idle * (timeSec/3600); % kg
end

%======================================================================
% ISA ATMOSPHERE (rho, a, T)
%======================================================================

function [rho,a,T] = ISA_atmosphere(h)
    T0=288.15; p0=101325; L=0.0065; R=287.058; g0=9.80665;
    if h < 11000
        T = T0 - L*h;
        p = p0*(T/T0)^(g0/(L*R));
    else
        T = 216.65;
        p11 = p0*(T/T0)^(g0/(L*R));
        p = p11*exp(-g0*(h-11000)/(R*T));
    end
    rho = p/(R*T);
    a = sqrt(1.4*R*T);
end

function CD = dragCoefficient(CL,phase)
    if strcmpi(phase,'cruise')
        CD0 = 0.025;
    else
        CD0 = 0.045;
    end
    K = 0.0435;
    CD = CD0 + K*CL^2;
end

%======================================================================
% GUI HELPERS
%======================================================================

function makeLabel(figHandle, str, y)
    uicontrol(figHandle,'Style','text','Position',[25 y 300 22], ...
        'HorizontalAlignment','left','String',str);
end

function hEdit = makeEdit(figHandle, defaultVal, y)
    hEdit = uicontrol(figHandle,'Style','edit','Position',[330 y 80 26], ...
        'String',num2str(defaultVal));
end

end

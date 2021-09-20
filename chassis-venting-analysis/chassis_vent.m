%
%   chassisVentingAnalysis.m
%   Calvin Boyle - 2021
%   Carnegie Mellon University
%

%   import Falcon 9 fairing venting data
%   data ripped from pg. 31 (Figure 4-10) https://www.mach5lowdown.com/wp-content/uploads/2020/04/falcon_9_users_guide_rev_2.0-1.pdf
F9_DATA = readtable('f9_vent_rate.csv');
f9_pressure = 100000; %Pa

%   MoonRanger chassis initial conditions
ch_pressure = 100000; %Pa
CH_VOLUME = 0.022; %m^3
CH_TEMP = 300; %K
ch_mass = idealGasMass(ch_pressure, CH_VOLUME, CH_TEMP);
global VENT_AREA
VENT_AREA = 0.000507; %m^2

%   analysis time-stepping
global T_STEP
T_STEP = .1; %sec
t = 1:T_STEP:150;

%   initialize arrays for data storage
f9_pressure_log = zeros(1, length(t));
global f9_vent_log
f9_vent_log = zeros(1, length(t));
ch_pressure_log = zeros(1, length(t));
delta_pressure_log = zeros(1, length(t));


%main loop
for c = 1:length(t)
    f9_pressure = f9_pressure + pressureStep(t(c), F9_DATA, c);
    f9_pressure_log(c) = pa2psi(f9_pressure);
    ch_mass = ch_mass - feltMassFlow(ch_pressure - f9_pressure);
    ch_pressure = idealGasPressure(ch_mass, CH_VOLUME, CH_TEMP);
    ch_pressure_log(c) = pa2psi(ch_pressure);
    delta_pressure_log(c) = ch_pressure_log(c) - f9_pressure_log(c);
end

%plot and formatting
    tiledlayout(2,2)
    nexttile
    plot(t, ch_pressure_log);
    hold on
    plot(t, f9_pressure_log);
    grid on
    set(gca,'Xtick',0:10:150)
    set(gca,'XtickLabel',0:20:150)
    title('Pressure environments During Ascent')
    xlabel('T+ (s)')
    ylabel('Absolute Pressure (psi)')
    legend({'MoonRanger Internal', 'F9 Fairing'}, 'Location', 'southwest')
    
    nexttile
    plot(t, f9_vent_log);
    grid on
    set(gca,'Xtick',0:10:150)
    set(gca,'XtickLabel',0:20:150)
    title('Falcon 9 Fairing Vent Rate')
    xlabel('T+ (s)')
    ylabel('Vent Rate (psi/s)')
    
    nexttile([1 2])
    plot(t, delta_pressure_log);
    grid on
    set(gca,'Xtick',0:10:150)
    set(gca,'XtickLabel',0:10:150)
    title('MoonRanger Internal Pressurization')
    xlabel('T+ (s)')
    ylabel('Absolute Pressure (psi)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%calculation of ideal gas pressure given mass and volume
%P (Pa)
%m (kg)
%V (m^3)
%T (K)
function P = idealGasPressure(m, V, T)
    R_AIR = 0.287058; %J/g-K
    P = (m * 1000 * R_AIR * T) / V;
end


%calculation of ideal gas mass given pressure and volume
%P (Pa)
%m (kg)
%V (m^3)
%T (K)
function m = idealGasMass(P, V, T)
    R_AIR = 0.287058; %J/g-K
    m = (P * V) / (R_AIR * T) / 1000;
end


%linear interpolation of felt mass flow
%m_del (kg)
%p_delta (Pa)
function m_del = feltMassFlow(p_delta)
    global T_STEP
    global VENT_AREA
    m_del = p_delta * 0.000260369 * VENT_AREA * T_STEP;
end


%linear interpolation of fairing vent rate
function v_rate = pressureStep(time, T, c)
    global T_STEP
    global f9_vent_log

    for row = 1:height(T)
       if time >= T{row, 1} && time <=T{row+1, 1}
           break
       end
    end
    
    lower_time = T{row, 1};
    lower_vent = T{row, 2};
    upper_time = T{row+1, 1};
    upper_vent = T{row+1, 2};
    
    tp = (time - lower_time) / (upper_time - lower_time);
    f9_vent_log(c) = (tp * (upper_vent - lower_vent) + lower_vent);
    v_rate =  f9_vent_log(c) * 6894.76 * T_STEP;
end


%converts pascals to psi
function psi = pa2psi(p)
    psi = p / 6894.76;
end
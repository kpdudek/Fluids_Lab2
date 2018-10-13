function Fluids_Lab2
%%%     BEGINNING OF ANALYSIS     %%%
[m,b] = calibration_fit();

Cross_Section_Velocities(m,b)


%%%   Table 2 of Packet: Static Pressure along the test section   %%%
volts_along_tunnel = [1.258,1.261,1.267,1.275,1.284,1.290,1.390,1.463,1.476,1.478];



function [m,b] = calibration_fit()
%%%   Values Taken for Calibration   %%%
v0 = 1.241; %Voltage at zero pressure | for calibration
h0 = .203;%in_m(.203); %Units in inches | height at zero pressure | for calibration
v1 = 1.452; %Voltage at calibration pressure
h1 = .519;%in_m(.519); %Units in inches | height at calibration pressure

v = [v0,v1];
h = [h0,h1];
p1 = ((h(2)-h(1))*.036 * 6894.75729);%*1.204*9.8);
press = [0,p1];

cal = fitlm(press,v,'linear');
b=table2array(cal.Coefficients(1,'Estimate'));
m=table2array(cal.Coefficients(2,'Estimate'));

figure('Visible','on','Name','Calibration Fit')
plot(press,v,'.','MarkerSize',20)
hold on
linfit = press.*m+b;
plot(press,linfit,'r')

legend('Voltage','Calibration Fit','Location','Northwest')
xlabel('Pressure (Pa)')
ylabel('Voltage (V)')
p_title = sprintf('Calibration Fit: V = \\DeltaP*%.4f + %.3f',m,b);
title(p_title)

function Cross_Section_Velocities(m,b)
v0 = 1.241;
%%%   Table 1 of Packet: Pressure and Velocity Profiles   %%%
port8_pose = [66,71,76,81,86,91,96]; %Units in mm
port8_volt = [1.452,1.451,1.451,1.452,1.453,1.454,1.460]; %Units in Volts

port1_pose = [8,18,28,38,48,58,68,78,88,98,108,118,128,138,148,158]; %Units in mm
port1_volt = [1.252,1.253,1.252,1.256,1.253,1.254,1.253,1.254,1.254,1.254,1.255,1.253,1.255,1.253,1.254,1.253];

port5_pose = [46,56,66,76,86,96,106,116]; %Units in mm
port5_volt = [1.279,1.279,1.279,1.280,1.280,1.280,1.280,1.281];

port10_pose = [30,40,50,60,70,80,90,100,110,120,130]; %Units in mm
port10_volt = [1.232,1.233,1.234,1.310,1.467,1.470,1.471,1.405,1.254,1.227,1.227];

%%% Port8
vel8 = zeros(1,length(port8_volt));
for i = 1:length(port8_volt)
    if port8_volt(i) > 0
        vel8(i) = sqrt(2*((port8_volt(i)-b)/m)/1.204);
    else
        vel8(i) = NaN;
    end
end

%%% Port1
vel1 = zeros(1,length(port1_volt));
for i = 1:length(port1_volt)
    if port1_volt(i) > v0
        vel1(i) = sqrt(2*((port1_volt(i)-b)/m)/1.204);
    else
        vel1(i) = NaN;
    end
end

%%% Port5
vel5 = zeros(1,length(port5_volt));
for i = 1:length(port5_volt)
    if port5_volt(i) > v0
        vel5(i) = sqrt(2*((port5_volt(i)-b)/m)/1.204);
    else
        vel5(i) = NaN;
    end
end

%%% Port10
vel10 = zeros(1,length(port10_volt));
for i = 1:length(port10_volt)
    if port10_volt(i) > v0
        vel10(i) = sqrt(2*((port10_volt(i)-b)/m)/1.204);
    else
        vel10(i) = NaN;
    end
end

figure('Visible','on','Name','Cross Section Velocities')
plot(vel8,port8_pose,'+r',vel1,port1_pose,'og',vel5,port5_pose,'*b',vel10,port10_pose,'xk')
legend('Point 8','Point 1','Point 5','Point 10')
xlabel('Velocity (m/s)')
ylabel('Vertical Position (mm)')
title('Vertical Position vs Velocity at 4 Cross sections')

function out = in_m(in)
out = in.*.0254;
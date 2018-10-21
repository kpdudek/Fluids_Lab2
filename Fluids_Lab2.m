function Fluids_Lab2
%%%     BEGINNING OF ANALYSIS     %%%
[m,b] = calibration_fit();
static = static_pressures(m,b);
[stag8,stag1,stag5,stag10] = stagnation_pressure(m,b);
[vel8,vel1,vel5,vel10] = Cross_Section_Velocities(1.204,stag8,stag1,stag5,stag10,static);
[ave_vel8,ave_vel1,ave_vel5,ave_vel10] = average_velocites(vel8,vel1,vel5,vel10);
table_calcs(static,ave_vel8,ave_vel1,ave_vel5,ave_vel10)



function [m,b] = calibration_fit()
%%%   Values Taken for Calibration   %%%
v0 = 1.241; %Voltage at zero pressure | for calibration
h0 = in_m(.203); %Units in inches | height at zero pressure | for calibration
v1 = 1.452; %Voltage at calibration pressure
h1 = in_m(.519); %Units in inches | height at calibration pressure

v = [v0,v1];
h = [h0,h1];
p1 = (h(2)-h(1))*1000*9.8*2; %.036 * 6894.75729);%*1.204*9.8);
press = [0,p1];

cal = fitlm(press,v,'linear');
b=table2array(cal.Coefficients(1,'Estimate'));
m=table2array(cal.Coefficients(2,'Estimate'));

figure('Visible','on','Name','Calibration Fit')
plot(press,v,'.','MarkerSize',25)
hold on
linfit = press.*m+b;
plot(press,linfit,'r')

legend('Voltage','Calibration Fit','Location','Northwest')
xlabel('Differential Pressure (Pa)')
ylabel('Voltage (V)')
p_title = sprintf('Calibration Fit: V = %.4fP + %.3f',m,b);
title(p_title)

function static = static_pressures(m,b)
%%%   Table 2 of Packet: Static Pressure along the test section   %%%
volts_along_tunnel = [1.258,1.261,1.267,1.275,1.284,1.290,1.390,1.463,1.476,1.478];

static = zeros(1,length(volts_along_tunnel));
for i = 1:length(volts_along_tunnel)
    static(i) = (volts_along_tunnel(i)-b)/m;
end

function [stag8,stag1,stag5,stag10] = stagnation_pressure(m,b)
port8_volt = [1.452,1.451,1.451,1.452,1.453,1.454,1.460]; %Units in Volts
port1_volt = [1.252,1.253,1.252,1.256,1.253,1.254,1.253,1.254,1.254,1.254,1.255,1.253,1.255,1.253,1.254,1.253];
port5_volt = [1.279,1.279,1.279,1.280,1.280,1.280,1.280,1.281];
port10_volt = [1.232,1.233,1.234,1.310,1.467,1.470,1.471,1.405,1.254,1.227,1.227];
v0 = 1.241;

vel8 = zeros(1,length(port8_volt));
stag8 = zeros(1,length(port8_volt));
for i = 1:length(port8_volt)
    if port8_volt(i) > v0
        stag8(i) = (port8_volt(i)-b)/m;
    else
        stag8(i) = NaN;
    end
end

%%% Port1
vel1 = zeros(1,length(port1_volt));
stag1 = zeros(1,length(port1_volt));
for i = 1:length(port1_volt)
    if port1_volt(i) > v0
        stag1(i) = (port1_volt(i)-b)/m;
    else
        stag1(i) = NaN;
    end
end

%%% Port5
vel5 = zeros(1,length(port5_volt));
stag5 = zeros(1,length(port5_volt));
for i = 1:length(port5_volt)
    if port5_volt(i) > v0
        stag5(i) = (port5_volt(i)-b)/m;
    else
        stag5(i) = NaN;
    end
end

%%% Port10
vel10 = zeros(1,length(port10_volt));
stag10 = zeros(1,length(port10_volt));
for i = 1:length(port10_volt)
    if port10_volt(i) > v0
        stag10(i) = (port10_volt(i)-b)/m;
    else
        stag10(i) = NaN;
    end
end

function [vel8,vel1,vel5,vel10] = Cross_Section_Velocities(p,stag8,stag1,stag5,stag10,static)
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

static = zeros(1,10);
%%% Port8
vel8 = zeros(1,length(port8_volt));
for i = 1:length(port8_volt)
    if port8_volt(i) > v0
        vel8(i) = sqrt(2*(stag8(i)-static(8))/p);
    else
        vel8(i) = NaN;
    end
end

%%% Port1
vel1 = zeros(1,length(port1_volt));
for i = 1:length(port1_volt)
    if port1_volt(i) > v0
        vel1(i) = sqrt(2*(stag1(i)-static(1))/p);
    else
        vel1(i) = NaN;
    end
end

%%% Port5
vel5 = zeros(1,length(port5_volt));
for i = 1:length(port5_volt)
    if port5_volt(i) > v0
        vel5(i) = sqrt(2*(stag5(i)-static(5))/p);
    else
        vel5(i) = NaN;
    end
end

%%% Port10
vel10 = zeros(1,length(port10_volt));
for i = 1:length(port10_volt)
    if port10_volt(i) > v0
        vel10(i) = sqrt(2*(stag10(i)-static(10))/p);
    else
        vel10(i) = NaN;
    end
end

figure('Visible','on','Name','Cross Section Velocities')
plot(vel8,port8_pose,'+r',vel1,port1_pose,'og',vel5,port5_pose,'*b',vel10,port10_pose,'xk')
xlim([0,20])
ylim([0,170])
legend('Point 8','Point 1','Point 5','Point 10')
xlabel('Velocity (m/s)')
ylabel('Vertical Position (mm)')
title('Vertical Position vs Velocity at 4 Cross sections')

function [ave_vel8,ave_vel1,ave_vel5,ave_vel10] = average_velocites(vel8,vel1,vel5,vel10)
h8 = .005;
h1 = .010;
h5 = .010;
h10 = .010;
disp(vel10)

sum_vels = [];
for i =1:length(vel8)-1
    sum_vels(end+1) = vel8(i)+vel8(i+1);
end
ave_vel8 = sum(sum_vels)/(2*(length(vel8)-1));

sum_vels = [];
for i =1:length(vel1)-1
    sum_vels(end+1) = vel1(i)+vel1(i+1);
end
ave_vel1 = sum(sum_vels)/(2*(length(vel1)-1));

sum_vels = [];
for i =1:length(vel5)-1
    sum_vels(end+1) = vel5(i)+vel5(i+1);
end
ave_vel5 = sum(sum_vels)/(2*(length(vel5)-1));

sum_vels = [];
for i =4:8
    sum_vels(end+1) = vel10(i)+vel10(i+1);
end
ave_vel10 = sum(sum_vels)/(2*5);

%disp([ave_vel1,ave_vel5,ave_vel8,ave_vel10])

function table_calcs(static,ave_vel8,ave_vel1,ave_vel5,ave_vel10)
areas = [33555,29643,25655,21897,17986,14151,10240,8323,13753,24830] .* 10^-6;
perims = [];
for i = 1:10
    l = [329,278,226,177,126,76,25,0,25,76];
    if i <= 8
        perims(end+1) = (2*203) + (2*(41+(2*l(i)*.1889)));
    else
        perims(end+1) = (2*203) + (2*(41+(2*l(i)*.535)));
    end
end
perims = perims*10^-3;
    
ave_vels = (areas(1)./areas).*ave_vel1;
ave_vels(5) = ave_vel5;
ave_vels(8) = ave_vel8;
ave_vels(10) = ave_vel10;

hd = (4.*areas)./perims;
re = (1.204.*ave_vels.*hd)./(18.13*10^-6);
%disp(re)

flow = areas.*ave_vels;
%disp(flow)

stag = static + (.5 * 1.204 .* ave_vels.^2);
%disp(stag)

for i = 1:10
    fprintf('At port %d:\nV_ave = % .3f\nRe = %.3f\nP_s = %.3f\nQ = %.3f\nP_t = %.3f\n\n',i,ave_vels(i),re(i),...
        static(i),flow(i),stag(i));
end

function out = in_m(in)
out = in.*.0254;

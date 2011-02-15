%FDTD1D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize MATLAB
close all; clc;
clear all; 

%Constants
c0 = 299792458; %m/s
e0 = 8.854187817*10^-12; %F/m
u0 = 1.256637061*10^-6; %H/m


%Physical Environment
dz = 1.4286*10^-8; %meters
dt = 4.7652*10^-17; %secs

%Simulated Environment
Nz = 180;
STEPS = 1000;

%Material Vectors
ER = ones([1 Nz]);
UR = ones([1 Nz]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FDTD Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Update Coefficients
mER = (c0*dt/dz)./ER;
mHR = (c0*dt/dz)./UR;

% Initialize Feilds
Ey = zeros([1 Nz]);
Hx = zeros([1 Nz]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:STEPS

  % Calculate H
  for nz = 1:Nz-1
    Hx(nz) = Hx(nz) + mHR(nz)*(Ey(nz+1)-Ey(nz));
  end
  Hx(Nz) = Hx(Nz) + mHR(Nz)*(0 - Ey(Nz));

  % Calculate E  
  Ey(1) = Ey(1) + mER(1)*(Hx(1) - 0);
  for nz = 2:Nz
    Ey(nz) = Ey(nz) + mER(nz)*(Hx(nz)-Hx(nz-1)); 
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure;
SetFigure(fig, 'HW#3-P2', [680 274 965 826]);

%Plot Magnetic Field
subplot(211)
h = plot(Hx, '-r', 'LineWidth', 2);
title('Magnetic Field');
h = get(h, 'Parent');
set(h, 'Fontsize', 14);
xlabel('z');
ylabel('Hx', 'Rotation', 0);
set(gca,'YTickLabel',{'1','0.5','0', '-0.5', '-1'})

%Plot Electric Field
subplot(212)
h = plot(Ey, '-b', 'LineWidth', 2);
title('Electric Field');
h = get(h, 'Parent');
set(h, 'Fontsize', 14);
xlabel('z');
ylabel('Ey', 'Rotation', 0);
set(gca,'YTickLabel',{'1','0.5','0', '-0.5', '-1'})





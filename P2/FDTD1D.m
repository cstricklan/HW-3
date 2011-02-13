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
mER = (c0*dt)./ER;
mHR = (c0*dt)./UR;

% Initialize Feilds
Ey = zeros([1 Nz]);
Hx = zeros([1 Nz]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:STEPS
  for nz = 1:Nz
    Hx1 = nz-1;
    Ey1 = nz+1;
    
    %Handle Boundary
    if(Hx1==0)
      Hx1=Nz;
    end
    
    if(Ey1>Nz)
      Ey1 = 1;
    end
    
    
    % Calculate H
    Hx(nz) = Hx(nz) + mHR(nz)*(Ey(Ey1)-Ey(nz))/dz;

    % Calculate E
    Ey(nz) = Ey(nz) + mER(nz)*(Hx(nz)-Hx(Hx1))/dz; 
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure;
SetFigure(fig, 'HW#3-P2', [680 274 965 826]);

%Plot Magnetic Field
subplot(211);ylabel('ytest');title('title1');
h = imagesc(Hx);
colorbar;
title('Magnetic Field');
h = get(h, 'Parent');
set(h, 'Fontsize', 14);
xlabel('z');
ylabel('', 'Rotation', 0);
set(gca,'YTickLabel',{'','',''})

%Plot Electric Field
subplot(212);ylabel('ytest2');title('title2');
h = imagesc(Ey);
colorbar;
title('Electric Field');
h = get(h, 'Parent');
set(h, 'Fontsize', 14);
xlabel('z');
ylabel('', 'Rotation', 0);
set(gca,'YTickLabel',{'','',''})





%% MATLAB code to implement numerical solutions to the two dimensional convection equation using Finite Volume Method
%% Equation: grad.(rho*<V>*phi) = 0;
%% Interpolation method: Deferred correction
%% Version: 1:0 | Date: 28th March 2020
%% Author : R Surya Narayan
%% PREAMBLE:
clear all;clc;close all;
%% Variable setup and declaration
%% Boundary and Mesh parameters
length = 1; %length along the positive x-directon of the flow domain 
breadth = 1;%length along the positive y-direction of the flow domain
mesh_x = 50; %number of cells discretized along the x-direction
mesh_y = 50; %number of cells discretized along the y-direction
dx = length/mesh_x; %cell size along the x-direction 
dy = breadth/mesh_y; %cell size along the y-direction
iterations = 10000; %number of iterations to be carried out
beta = 1; %Blending factor determining proportion of differencing and upwinding interpolation applied
%% Fluid parameters
u = 2; %x-velocity
v = 2; %y-velocity
rho = 1; %density of the fluid
m_e = rho*u*dy; % east constant for a unit cell
m_w = rho*u*dy; % west constant for a unit cell
m_n = rho*v*dx; % north constant for a unit cell
m_s = rho*v*dx; % south constant for a unit cell
phi_w_bc = 1; %west boundary condition
phi_s_bc = 0; %south boundary condition
%% Solution variables
phi = zeros(mesh_x+2,mesh_y+2); %convection scalar 
AE = zeros(mesh_x,mesh_y); %East coefficient for a cell considering upwinding interpolation scheme
AW = zeros(mesh_x,mesh_y); %West coefficient for a cell considering upwinding interpolation scheme
AN = zeros(mesh_x,mesh_y); %North coefficient for a cell considering upwinding interpolation scheme
AS = zeros(mesh_x,mesh_y); %South coefficient for a cell considering upwinding interpolation scheme
SU = zeros(mesh_x,mesh_y); %Source term for a cell considering upwinding interpolation scheme
SP = zeros(mesh_x,mesh_y); %Source term coefficient for a cell center considering upwinding interpolation scheme
AEC = zeros(mesh_x,mesh_y); %East coefficient for a cell considering central differencing interpolation scheme
AWC = zeros(mesh_x,mesh_y); %West coefficient for a cell considering central differencing interpolation scheme
ANC = zeros(mesh_x,mesh_y); %North coefficient for a cell considering central differencing interpolation scheme
ASC = zeros(mesh_x,mesh_y); %South coefficient for a cell considering central differencing interpolation scheme
SUC = zeros(mesh_x,mesh_y); %Source term for a cell considering upwinding interpolation scheme
SPC = zeros(mesh_x,mesh_y); %Source term coefficient for a cell center considering upwinding interpolation scheme
%% Variable Initialization
%Upwinding coefficients for interior cells 
AE(:,:) = 0;
AW(:,:) = m_w;
AN(:,:) = 0;
AS(:,:) = m_s;
%central differencing coefficients for interior cells
AEC(:,:) = -m_e/2;
AWC(:,:) = m_w/2;
ANC(:,:) = -m_n/2;
ASC(:,:) = m_s/2;
%% Boundary conditions
%upwinding coefficients
%west boundary
AW(:,1) = 0;
SP(:,1) = -m_e;
SU(:,1) = m_w*phi_w_bc;
%east boundary
AE(:,mesh_x) = 0;
SP(:,mesh_x) = -(m_e-m_w);
%south boundary
AS(1,:) = 0;
SP(1,:) = -m_n;
SU(1,:) = m_s*phi_s_bc;
%north boundary
AN(mesh_y,:) = 0;
SP(mesh_y,:) = -(m_n-m_s);
%Central differencing coefficients
%west boundary
AWC(:,1) = 0;%cut the link
SPC(:,1) = -m_e;
SUC(:,1) = m_w*phi_w_bc;
%east boundary
AEC(:,mesh_x) = 0;
SPC(:,mesh_x) = -(m_e-m_w);
%south boundary
ASC(1,:) = 0;
SPC(1,:) = -m_n;
SUC(1,:) = m_s*phi_s_bc;
%north boundary
ANC(mesh_y,:) = 0;
SPC(mesh_y,:) = -(m_n-m_s);
%% corner cell coefficients
%central differencing 
%southwest corner
AWC(1,1) = 0;
ASC(1,1) = 0;
SPC(1,1) = -m_e-m_n;
SUC(1,1) = m_w*0 + m_s*phi_s_bc;
%northwest corner
ANC(mesh_y,1) = 0;
AWC(mesh_y,1) = 0;
SUC(mesh_y,1) = m_w*phi_w_bc;
SPC(mesh_y,1) = -m_e-(m_n-m_s);
%southeast corner
ASC(1,mesh_x) = 0;
AEC(1,mesh_x) = 0;
SPC(1,mesh_x) = -m_n-(m_e-m_w);
%northeast corner
ANC(mesh_y,mesh_x) = 0;
AEC(mesh_y,mesh_x) = 0;
SPC(mesh_y,mesh_x) = -(m_n-m_s)-(m_e-m_n);
%upwinding coefficients 
%southwest corner
AW(1,1) = 0;
AS(1,1) = 0;
SP(1,1) = -m_e-m_n;
SU(1,1) = m_w*0 + m_s*phi_s_bc;
%northwest corner
AN(mesh_y,1) = 0;
AW(mesh_y,1) = 0;  
SU(mesh_y,1) = m_w*phi_w_bc;
SP(mesh_y,1) = -m_e-(m_n-m_s);
%southeast corner
AS(1,mesh_x) = 0;
AE(1,mesh_x) = 0;
SP(1,mesh_x) = -m_n-(m_e-m_w);
%northeast corner
AN(mesh_y,mesh_x) = 0;
AE(mesh_y,mesh_x) = 0;
SP(mesh_y,mesh_x) = -(m_n-m_s)-(m_e-m_n);

%finalize AP and APC b  efore solving
AP = AE+AW+AS+AN-SP;
APC = AEC+AWC+ASC+ANC-SPC;
%% ITERATIVE SOLVER-Deferred correction
for k = 1:iterations
    for i = 1:mesh_y
        for j= 1:mesh_x
        phi(i+1,j+1) = (( AE(i,j)*phi(i+1,j+2) + AW(i,j)*phi(i+1,j) + AN(i,j)*phi(i+2,j+1) + AS(i,j)*phi(i,j+1) + SU(i,j))/AP(i,j))...
            -(beta/AP(i,j))*((APC(i,j)*phi(i+1,j+1) - AEC(i,j)*phi(i+1,j+2) - AWC(i,j)*phi(i+1,j) - ANC(i,j)*phi(i+2,j+1)...
            - ASC(i,j)*phi(i,j+1) - SUC(i,j)) - ( AP(i,j)*phi(i+1,j+1) - AE(i,j)*phi(i+1,j+2) - AW(i,j)*phi(i+1,j) - AN(i,j)*phi(i+2,j+1) -AS(i,j)*phi(i,j+1)- SU(i,j)));
        end
    end
end
%% POST PROCESSING
%% WRITE THE CSV FILE FOR PARAVIEW (OR) GENERATE X,Y and Phi vectors to get a surface plot of contours
csv_matrix = []; %matrix that stores coordinates and results at those coordinates
X = [0];Y=[0];
for i = 1:mesh_x
    x = dx/2 + (i-1)*dx;
    X = [X,x];
end
for j = 1:mesh_y
    y = dy/2 + (j-1)*dy;
    Y = [Y,y];
end
X = [X,X(end)+dx/2];
Y = [Y,Y(end)+dy/2];
for i = 1:mesh_x
    for j = 1:mesh_y
        x = dx/2 + (i-1)*dx;
        y = dy/2 + (j-1)*dy;
        z = 0;
        csv_matrix = [csv_matrix;x, y, z, phi(j+1,i+1)];
    end
end
phi(1,:) = phi_s_bc; %south boundary line
phi(:,1) = phi_w_bc; %west boundary line
phi(mesh_y+2,:) = phi(mesh_y+1,:);%north boundary line
phi(1,mesh_x+2) = phi(1,mesh_x+1);%east boundary line
%append the boundary values
%west boundary
for i = 1:mesh_y
    x=0;
    y = dy/2 + (i-1)*dy;
    z=0;
    csv_matrix = [csv_matrix;x,y,z,phi_w_bc];
end
%south boundary
for i= 1:mesh_x
    x = dx/2 + (i-1)*dx;
    y = 0;
    z = 0;
    csv_matrix = [csv_matrix;x,y,z,phi_s_bc];
end
%east boundary
for i = 1:mesh_y
    x=1;
    y = dy/2 + (i-1)*dy;
    z=0;
    csv_matrix = [csv_matrix;x,y,z,phi(i+1,mesh_x+1)];
end
%north boundary
for i= 1:mesh_x
    x = dx/2 + (i-1)*dx;
    y = 1;
    z = 0;
    csv_matrix = [csv_matrix;x,y,z,phi(mesh_x+1,i+1)];
end
%corner values
csv_matrix = [csv_matrix;0,0,0,phi(1,mesh_x+2)];%south east
csv_matrix = [csv_matrix;1,0,0,phi(1,1)];%south west
csv_matrix = [csv_matrix;0,1,0,phi(mesh_y+2,mesh_x+2)];%north east
csv_matrix = [csv_matrix;1,1,0,phi(mesh_y+2,1)];%north west
%write into a text file to later convert to a CSV
%writematrix(csv_matrix,'convection_onebeta_sw_modified.txt') uncomment
%this to write the matrix into a text file
%MATLAB postprocessing stuff
contourf(X,Y,phi,'ShowText','on');
colormap(jet);
title('Convection plot');





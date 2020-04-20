%% A MATLAB Code to discretize and solve numerically the two-dimensional form of the diffusion equation 
%% DISCRETIZATION PRINCIPLE: FVM - Finite Volume Method
%% EQUATION: grad.(gamma*grad(phi))+ source = 0 or grad^2(phi*gamma)+source =0
%% VERSION 2.0 | Date:27th March 2020
%% SOLVING METHOD USED FOR EQUATIONS: Successive Under Relaxation method of iterations 
%% Author: R Surya Narayan
%% PREAMBLE
clear all; close all; clc;
%% VARIABLE SETUP and DECLARATION
%% Boundary and Mesh Parameters
length = 1; %length along the positive x-directon of the flow domain 
breadth = 1;%length along the positive y-direction of the flow domain
mesh_x = 25; %number of cells discretized along the x-direction
mesh_y = 25; %number of cells discretized along the y-direction
dx = length/mesh_x; %cell size along the x-direction
dy = breadth/mesh_y; %cell size along the y-direction
iterations = 3000; %number of iterations to be carried out
%% Variables to solve for including coefficients
gamma = 1; %thermal conductivity
omega = 1.7; %the constant for the SUR method
phi = zeros(mesh_x+2,mesh_y+2); %the temperature function-set to zero at all points initially as a guess 
A_n = ones(mesh_x,mesh_y); %north coefficient
A_s = ones(mesh_x,mesh_y); %south coefficient
A_e = ones(mesh_x,mesh_y); %east coefficient
A_w = ones(mesh_x,mesh_y); %west coefficient
A_p = zeros(mesh_x,mesh_y);%coefficient of p
S_u = zeros(mesh_x,mesh_y); %Constant Source term
S_p = zeros(mesh_x,mesh_y); %Source term as a coefficient of the phi_P term
%% INITITIALIZE BOUNDARY CONDITIONS ON COEFFICIENTS
%west boundary condition
A_w(:,1) = 0; %cut the link
S_u(:,1) = 0; %use the source terms 
S_p(:,1) =-2;
%east boundary condition
A_e(:,mesh_x)=0; %cut the link
S_p(:,mesh_x)=-2; %use the source terms
for i = 1:mesh_y
    S_u(i,mesh_x) = 2*((dy/2) + (i-1)*(dy)); 
end
%north boundary condition
A_n(mesh_y,:)=0; %cut the link
S_p(mesh_y,:)=-2;
for j=1:mesh_x
    S_u(mesh_y,j) = 2*((dx/2)+(j-1)*(dx));
end
%south boundary condition
A_s(1,:)=0; 
S_u(1,:)=0;
S_p(1,:)=-2;
%corner cell boundary conditions
%south-west corner cell
A_w(1,1)=0; %cut the west link
A_s(1,1)=0; %cut the south link
S_u(1,1)=0; %use the source terms
S_p(1,1)=-4; %-2 from the south and -2 from the west 
%south-east corner cell
A_s(1,mesh_x)=0; %cut the south link
A_e(1,mesh_x)=0; %cut the east link
S_u(1,mesh_x)=2*(dy/2);
S_p(1,mesh_x)=-4;
%north-east corner cell
A_n(mesh_y,mesh_x)=0;%cut the north link
A_e(mesh_y,mesh_x)=0;%cut the east link
S_u(mesh_y,mesh_x)=2*((dx/2)+(mesh_x-1)*(dx)) + 2*((dy/2)+(mesh_y-1)*(dy));
S_p(mesh_y,mesh_x)=-4;
%north-west boundary condition
A_n(mesh_y,1)=0; %cut the north link
A_w(mesh_y,1)=0; %cut the west link
S_u(mesh_y,1)= 2*(dx/2);
S_p(mesh_y,1)=-4;
A_p = A_s + A_w + A_n + A_e-S_p;
%% ITERATIVE SOLVER
%outer k-loop to keep track of the number of iterations
for k = 1:iterations
    for i = 1:mesh_x
        for j=1:mesh_y
            phi(i+1,j+1) = phi(i+1,j+1)+(omega/A_p(i,j))*( (A_w(i,j)*phi(i+1,j)) + (A_e(i,j)*phi(i+1,j+2)) + (A_n(i,j)*phi(i+2,j+1)) + (A_s(i,j)*phi(i,j+1))+ S_u(i,j)- (A_p(i,j)*phi(i+1,j+1)));
        end
    end
end
%% INITIALIZE TRUE SOLUTION 
true_phi = zeros(mesh_x+2,mesh_y+2);
for j = 2:mesh_x+1
    for i = 2:mesh_y+1
        true_phi(i,j) = (dx/2 + (j-2)*dx)*(dy/2 + (i-2)*dy);
    end
end
%phi on the east boundary line
true_phi(mesh_x+2,mesh_y+2) = length*breadth;
true_phi(1,mesh_y+2) = dy/2;
for i = 3:mesh_x+1
    true_phi(i,mesh_y+2) = (dy/2)+(i-2)*dy;
end
%phi on the north boundary line
true_phi(mesh_x+2,1)= dx/2;
for j = 3:mesh_y+1
    true_phi(mesh_x+2,j) = (dx/2)+(j-2)*dx;
end
true_phi(mesh_x+2,mesh_y+2) = 1;
%% WRITE THE CSV FILE FOR PARAVIEW
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
%the following lines of code are for exporting the data from MATLAB to a
%text file that can be later used for converting it into a CSV file
for i = 1:mesh_x
    for j = 1:mesh_y
        x = dx/2 + (i-1)*dx;
        y = dy/2 + (j-1)*dy;
        z = 0;
        csv_matrix = [csv_matrix;x, y, z, phi(i+1,j+1)];
    end
end
phi(1,:) = true_phi(1,:); %south boundary line
phi(:,1) = true_phi(:,1); %west boundary line
phi(mesh_y+2,:) = true_phi(mesh_y+2,:); %north boundary line
phi(:,mesh_x+2) = true_phi(:,mesh_x+2); %east boundary line
err = abs(true_phi-phi);
%append the boundary values
%west boundary
for i = 1:mesh_y
    x=0;
    y = dy/2 + (j-1)*dy;
    z=0;
    csv_matrix = [csv_matrix;x,y,z,0];
end
%south boundary
for i= 1:mesh_x
    x = dx/2 + (i-1)*dx;
    y = 0;
    z = 0;
    csv_matrix = [csv_matrix;x,y,z,0];
end
%east boundary
for i = 1:mesh_y
    x=0;
    y = dy/2 + (j-1)*dy;
    z=0;
    csv_matrix = [csv_matrix;x,y,z,y];
end
%north boundary
for i= 1:mesh_x
    x = dx/2 + (i-1)*dx;
    y = 0;
    z = 0;
    csv_matrix = [csv_matrix;x,y,z,x];
end
%corner values
csv_matrix = [csv_matrix;0,0,0,0];%south east
csv_matrix = [csv_matrix;1,0,0,0];%south west
csv_matrix = [csv_matrix;0,1,0,0];%north east
csv_matrix = [csv_matrix;1,1,0,1];%north west
%writematrix(csv_matrix) uncomment for writing the CSV file
%type 'CSV_File _of_coordinates.txt'; this too
contourf(X,Y,phi,'ShowText','on');
colormap(jet);
title('Diffusion plot');
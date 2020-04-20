%A MATLAB code to numerically simulate the laminar flow of a fluid through
%a rectangular pipe numerically
clear all;clc;
%% Variable initialization
hx = 0.004; %Size of each mesh element or grid-spacing along x-direction
hy = 0.004; %size of each mesh element or grid-spacing along y-direction
PressureGrad = 1; %Pressure gradient along the z-direction
mu = 0.01; %dynamic viscosity in SI units (poise)
c = -(PressureGrad)/mu;
length = 0.4; %length of the cross-section of the pipe in meter
breadth = 0.4; %breadth of the cross-section of the pipe in meter
Ni = length/hx; %counter limit for i
Nj = breadth/hy; %counter limit for j
Gridpoints = (Ni-1)*(Nj-1); %total number of variables that we have to solve for 
w = zeros(Nj+1,Ni+1); %velocity of flow at different points
ae = zeros(Nj+1,Ni+1);
aw = zeros(Nj+1,Ni+1);
an = zeros(Nj+1,Ni+1);
as = zeros(Nj+1,Ni+1);
%% Boundary Conditions
%on the velocity
w(1,:)=0;
w(Nj+1,:)=0;
w(:,Ni+1)=0;
w(:,1)=0;
%initialize the coefficients
;
ae(:,:)= 1/((hx*hx)*2*((1/(hx^2))+(1/(hy^2))));
aw(:,:)= 1/((hx*hx)*2*((1/(hx^2))+(1/(hy^2))));
an(:,:)= 1/((hy*hy)*2*((1/(hx^2))+(1/(hy^2))));
as(:,:)= 1/((hy*hy)*2*((1/(hx^2))+(1/(hy^2))));
%apply boundary conditions
%south boundary
as(1,:) = 0; %cut the links
%north boundary
an(Nj+1,:) =0; %cut the links
%east boundary
ae(Nj+1,:) = 0; %cut the links
%west boundary
aw(1,:)=0; %cut the links
%% Solver using Gauss siedel
iterations =1000
for k = 1:iterations
    for i = 2:Ni
        for j=2:Nj
            w(i,j) = ae(i,j)*w(i,j+1) + aw(i,j)*w(i,j-1)+an(i,j)*w(i+1,j)+as(i,j)*w(i-1,j)-c/(2*((1/(hx^2))+(1/(hy^2))));
        end
    end
end
%% Post processing and generating the image for the area of cross-section
X = [0:hx:length];
Y = [0:hy:breadth];
contourf(X,Y,w,'ShowText','on');
colormap(jet);
title('Viscous flow through a square pipe');


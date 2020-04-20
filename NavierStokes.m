%% MATLAB CODE TO SOLVE THE TWO DIMENSIONAL STEADY STATE INCOMPRESSIBLE FORM OF THE NAVIER STOKES EQUATIONS
%%WITH MASS CONSERVATION COUPLING USING Semi-Implicit Method for Pressure
%%linked Coupling (SIMPLE) Algorithm 
%% Author: R Surya Narayan
%% Date: 8th April 2020 | VERSION V1.0
%% Discretization pronciple: FVM (Finite Volume Method)
%% EQUATION SOLVED: grad.(<V>) = 0; (Continuity/mass conservation) and steady incompressible Navier Stokes- grad.(rho*Vj<V>)=-grad(P)+grad.(mu*grad(Vj))
%% STAGGERED GRID Scheme used to represent different control volumes for velocity and pressure
%% Interpolation method used for velocities: First order accurate upwinding scheme
%% Problem: Flow of liquid through a driven cavity
%% PREAMBLE
clear all; clc; close all;
%% Problem setup and Variable declaration
%% Boundary and Mesh Parameters
length = 1; %length along the positive x-directon of the flow domain 
breadth = 1;%length along the positive y-direction of the flow domain
mesh_x = 60; %number of cells discretized along the x-direction
mesh_y = 60; %number of cells discretized along the y-direction
dx = length/mesh_x; %cell size along the x-direction
dy = breadth/mesh_y; %cell size along the y-direction
% fluid properties
mu = 1; %viscosity of the fluid 
rho = 1000; %density of the fluid
%iteration parameters and relaxation factors
omega_u = 0.7; %relaxation parameter for u momentum 
omega_v = 0.7; %relaxation parameter for v momentum
omega_p = 0.3; %relaxation parameter for pressure correction
outer_iterations = 1000; %number of times SIMPLE is going to be iterated through for convergence
iter_v = 10; %number of iterations for u and v momentum solvers
iter_p = 100; %number of iterations for pressure correction solver
%% Variables to solve for including coefficients
u = zeros(mesh_y+2,mesh_x+1);%u momentum
v = zeros(mesh_y+1,mesh_x+2);%v momentum
u_old = zeros(mesh_y+2,mesh_x+1); %guessed x-momentum on the cell centers
v_old = zeros(mesh_y+1,mesh_x+2); %guessed y-momentum on the cell centers
p_prime = zeros(mesh_y+2,mesh_x+2); %pressure correction
pressure = zeros(mesh_y+2,mesh_x+2);%pressure term
apu = ones(mesh_y+2,mesh_x+2);%coefficient of p-term in u momentum
apv = ones(mesh_y+2,mesh_x+2);%coefficient of p-term in v momentum
app = ones(mesh_y+2,mesh_x+2);%coefficient of p-term for pressure corrections
ae = zeros(mesh_y+2,mesh_x+2);%east coefficient for velocities
as = zeros(mesh_y+2,mesh_x+2);%south coefficient for velocities
an = zeros(mesh_y+2,mesh_x+2);%north coefficient for velocities
aw = zeros(mesh_y+2,mesh_x+2);%west coefficient for velocities
source = zeros(mesh_y+2,mesh_x+2); %source term at each stage to check mass conservation
%initialize variables for the driven cavity problem
u(mesh_y+2,:) = 5; %u velocity is fixed along the north boundary
%other four sides are the walls hence the velocity is zero, this happens to
%be our default declaration
u_old = u;
v_old = v;
%% Solver using SIMPLE-algorithm
%Main outer loop
for k = 1:outer_iterations
    %X-Momentum solver
    %initialize coefficients for the u-momentum the centers 
    for j = 2:mesh_x
        for i=2:mesh_y+1
            ae(i,j) = max(-0.5*rho*dy*(u_old(i,j)+u_old(i,j+1)),0) + mu*dy/dx;
            aw(i,j) = max(0.5*rho*dy*(u_old(i,j)+u_old(i,j-1)),0) + mu*dy/dx;
            an(i,j) = max(-0.5*rho*dx*(v_old(i,j)+v_old(i,j+1)),0) + mu*dx/dy;
            as(i,j) = max(0.5*rho*dx*(v_old(i-1,j)+v_old(i-1,j+1)),0) + mu*dx/dy;
        end
    end
    %correct the boundary values on the north and south boundary for the
    %X-momentum
    for j=2:mesh_x
        an(mesh_x+1,j) = max(-0.5*rho*dx*(v_old(mesh_y+1,j)+v_old(mesh_y+1,j+1)),0) + mu*dx/(dy/2);
        as(2,j) = max(0.5*rho*dx*(v_old(1,j)+v_old(1,j+1)),0) + mu*dx/(dy/2);
    end
    apu = ae+aw+as+an;
    %apu(15,14) = 1e30; to block out cells
    apu = apu/omega_u;
    %iterate on the x-momentum equations
    for var = 1:iter_v
        for j = 2:mesh_x
            for i = 2:mesh_y+1
                u(i,j) = (1-omega_u)*u_old(i,j) + (1/apu(i,j))*...
                    (ae(i,j)*u(i,j+1) + aw(i,j)*u(i,j-1) + an(i,j)*u(i+1,j)+ as(i,j)*u(i-1,j) + dy*(pressure(i,j)-pressure(i,j+1)));
            end
        end
    end
    %Y-momentum solver
    %write the values of coefficients
    for j = 2:mesh_x+1
        for i = 2:mesh_y
            ae(i,j) = max(-0.5*rho*dy*(u_old(i+1,j)+u_old(i,j)),0)+ mu*dy/dx;
            aw(i,j) = max(0.5*rho*dy*(u_old(i+1,j-1)+u_old(i,j-1)),0) + mu*dy/dx;
            an(i,j) = max(-0.5*rho*dx*(v_old(i,j)+v_old(i+1,j)),0) + mu*dx/dy;
            as(i,j) = max(0.5*rho*dx*(v_old(i,j)+v_old(i-1,j)),0) + mu*dx/dy;
        end
    end
    %overwrite east west boundary values
    for i = 2:mesh_y
        ae(i,mesh_y+1) = max(-0.5*rho*dy*(u_old(i+1,mesh_y+1)+u_old(i,mesh_y+1)),0) + mu*dy/(dx/2);
        aw(i,2) = max(0.5*rho*dy*(u_old(i+1,1)+u_old(i,1)),0) + mu*dy/(dx/2);
    end
    apv = ae+aw+as+an;
    %iterative solver
    apv = apv/omega_v;
    for var = 1:iter_v
        for j = 2:mesh_x+1
            for i = 2:mesh_y
                v(i,j) = v_old(i,j)*(1-omega_v) + (1/apv(i,j))*...
                    (ae(i,j)*v(i,j+1) + aw(i,j)*v(i,j-1) + as(i,j)*v(i-1,j)+ an(i,j)*v(i+1,j) + dx*(pressure(i,j)-pressure(i+1,j)));
            end
        end
    end
    %pressure correction equation
    %assign values to the coeffi+cients
    for i = 2:mesh_y+1
        for j = 2:mesh_x+1
            ae(i,j) = (rho*dy^2)/apu(i,j);
            aw(i,j) = (rho*dy^2)/apu(i,j-1);
            as(i,j) = (rho*dx^2)/apv(i-1,j);
            an(i,j) = (rho*dx^2)/apv(i,j);
        end
    end
    %pressure corrections shouldn't be applied on the boundary velocities
    %hence set the values to zero
    ae(:,mesh_x+1) = 0;
    aw(:,2) = 0;
    an(mesh_y+1,:)=0;
    as(2,:) = 0;
    app = ae+aw+as+an;
    app(2,2) = 1e30;
    %compute the source term
    source=0;
    for i = 2:mesh_y+1
        for j = 2:mesh_x+1
            source(i,j) = rho*dy*(u(i,j)-u(i,j-1)) + rho*dx*(v(i,j)-v(i-1,j));
        end
    end
    if rem(k,100)==0
        fprintf('%0.4f',total);
        fprintf('\n');
    end
    %SOR iterations for pressure correction
    for var = 1:1100
        for j = 2:mesh_x+1
            for i = 2:mesh_y+1
                p_prime(i,j) = p_prime(i,j) + (1.7/app(i,j))*...
                (ae(i,j)*p_prime(i,j+1) + aw(i,j)*p_prime(i,j-1)...
                + an(i,j)*p_prime(i+1,j) + as(i,j)*p_prime(i-1,j)...
                -source(i,j) - (p_prime(i,j)*app(i,j)));
            end
        end
    end
    %apply the pressure corrections
    for j = 2:mesh_x+1
        for i = 2:mesh_y+1
            pressure(i,j) = pressure(i,j)+omega_p*p_prime(i,j);
        end
    end
    %apply these corrections to get new velocities
    %u-momentum
    for i = 2 :mesh_y+1
        for j = 2:mesh_x
            u(i,j) = u(i,j) + (dy/apu(i,j))*(p_prime(i,j)-p_prime(i,j+1));
        end
    end
    %v-momentum
    for i = 2:mesh_y
        for j = 2:mesh_y+1
            v(i,j) = v(i,j)+(dx/apv(i,j))*(p_prime(i,j)-p_prime(i+1,j));
        end
    end
    %update velocities
    u_old = u;
    v_old = v;
    %recompute source term to check mass conservation
     for i = 2:mesh_y+1
        for j = 2:mesh_x+1
            source(i,j) = rho*dy*(u(i,j)-u(i,j-1)) + rho*dx*(v(i,j)-v(i-1,j));
        end
    end
    total = sqrt(sum(source.^2,'all'));
    %monitor mass conservation
    if rem(k,100)==0
        fprintf('%0.4f',total);
        fprintf('\n');
    end
end
%END of the NAVIER STOKES SOLVER
%% POST PROCESSING AND DISPLAY PLOTS
%quiver-plot: a type of plot where the velocities are represented with
%arrows give (x,y) components at a given coordinate or (u,v) at (x,y)
%hence interpolate staggered mesh-velocities to the corners of the main
%control volume after computing X, Y vectors as cell corner coordinates

%X-coordinates computation
X = [];
for i = 1:mesh_y+1
        x = (i-1)*dx;
        X = [X,x];
end

%Y-coordinates computation
Y = [];
for j = 1:mesh_x+1
        y = (j-1)*dy;
        Y = [Y,y];
end

%% X-velocity at grid-points computation
%interior grid points
for i = 2:mesh_y
    for j = 2:mesh_x
        U(i,j) = (u(i+1,j)+u(i,j))/2;
    end
end
%boundaries excluding corners
%north boundary
U(mesh_y+1,2:mesh_x) = u(mesh_y+2,2:mesh_x);
%south boundary
U(1,2:mesh_x) = u(1,2:mesh_x);
%east boundary 
U(2:mesh_y,mesh_x+1) = (u(2:mesh_y,mesh_x+1)+ u(3:mesh_y+1,mesh_x+1))/2;
%west boundary
U(2:mesh_y,1) = (u(2:mesh_y,1)+ u(3:mesh_y+1,1))/2;
%corner points
%south-west
U(1,1) = 0;
%south-east
U(mesh_x+1,1) = 0;
%north-west
U(mesh_y+1,1)=0;
%north-east
U(mesh_y+1,mesh_x+1) = 0;

%% Y-velocity at grid-points computation
%interior grid points
for i = 2:mesh_y
    for j = 2:mesh_x
        V(i,j) = (v(i,j+1)+v(i,j))/2;
    end
end
%boundaries excluding the corners
%north 
V(mesh_y+1,2:mesh_x) = 0.5*(v(mesh_y+1,2:mesh_x)+ v(mesh_y+1,3:mesh_x+1));
%south 
V(1,2:mesh_x) = 0.5*(v(1,2:mesh_x)+ v(1,3:mesh_x+1));
%east
V(2:mesh_y,mesh_x+1) = 0.5*(v(2:mesh_y,mesh_x+1)+v(3:mesh_y+1,mesh_x+1));
%west
V(2:mesh_y,1) = 0.5*(v(2:mesh_y,1)+v(3:mesh_y+1,1));
%corner grid-points
%southwest
V(1,1) = 0;
%southeast
V(1,mesh_x+1) = 0;
%northwest
V(mesh_y+1,1) = 0;
%northeast
V(mesh_y+1,mesh_x+1) = 0;
%quiverplot for velocities
figure(1);
quiver(X,Y,U,V,'Color','k','LineWidth',2);
%% pressure terms computation
P = ones(mesh_y+1,mesh_x+1);
%interior points
for i = 2:mesh_y
    for j=2:mesh_x
        P(i,j) = 0.25*(pressure(i,j) + pressure(i+1,j) + pressure(i,j+1) + pressure(i+1,j+1));
    end
end
%boundaries except corners
%west boundary
P(2:mesh_y,1) = 0.5*(pressure(2:mesh_y,1)+ pressure(3:mesh_y+1,1));
%east 
P(2:mesh_y,mesh_x+1) = 0.5*(pressure(2:mesh_y,mesh_x+1) + pressure(3:mesh_y+1,mesh_x+1));
%north
P(mesh_y+1,2:mesh_x) = 0.5*(pressure(mesh_y+1,2:mesh_x) + pressure(mesh_y+1,3:mesh_x+1));
%south 
P(1,2:mesh_x) = 0.5*(pressure(1,2:mesh_x) + pressure(1,3:mesh_x+1));
%corner points
%southwest
P(1,1) = pressure(2,2);
%southeast
P(1,mesh_x+1) = pressure(2,mesh_x+1);
%northeast
P(mesh_y+1,mesh_x+1) = pressure(mesh_y+1,mesh_x+1);
%northwest
P(mesh_y+1,1) = pressure(mesh_y+1,2);
figure(2);
contourf(X,Y,U,'ShowText','off','LineColor','none','LevelStep',0.01);
colormap(jet);
title('X-velocity Contours for lid-driven cavity flow');
figure(3);
contourf(X,Y,V,'ShowText','off','LineColor','none','LevelStep',0.01);
colormap(jet);
title('Y-velocity Contours for lid-driven cavity flow');
figure(4);
contourf(X,Y,P,'ShowText','off','LineColor','none','LevelStep',0.01);
colormap(jet);
title('Pressure Contours for lid-driven cavity flow');



%% MATLAB CODE TO SOLVE THE TWO DIMENSIONAL STEADY STATE INCOMPRESSIBLE FORM OF THE NAVIER STOKES EQUATIONS
%%Pressure-Velocity coupling done using SIMPLE Algorithm 
%% FLOW DOMAIN: Rectangular Channel 
%% Author: R Surya Narayan
%% Date: 2nd Feb 2020 | VERSION V1.2
%% Discretization pronciple: FVM (Finite Volume Method)
%% EQUATION SOLVED: grad.(<V>) = 0; (Continuity/mass conservation) and steady incompressible Navier Stokes- grad.(rho*Vj<V>)=-grad(P)+grad.(mu*grad(Vj))
%% STAGGERED GRID Scheme used to represent different control volumes for velocity and pressure
%% Interpolation method used for velocities: First order accurate upwinding scheme
%% PREAMBLE
clear all; clc; close all;
%% PROBLEM SETUP: set prob_type using the following legend:
% 1: fully developed channel flow
% 2: Pressure driven channel flow
% 3: Lid-driven cavity
% iblock: set to 1 (to block cells out)
prob_type=3;
iblock=0;
%% Boundary and Mesh Parameters
if (prob_type==3)
    %Case: Lid-driven cavity
    %ensure a square domain for the lid-driven cavity
    length = 1; %length along the positive x-directon of the flow domain 
    breadth = 1;%length along the positive y-direction of the flow domain
    mesh_x = 60; %number of cells discretized along the x-direction
    mesh_y = 60; %number of cells discretized along the y-direction
else
    length = 5; %length along the positive x-directon of the flow domain 
    breadth = 1;%length along the positive y-direction of the flow domain
    mesh_x = 200; %number of cells discretized along the x-direction
    mesh_y = 100; %number of cells discretized along the y-direction
end
dx = length/mesh_x; %cell size along the x-direction
dy = breadth/mesh_y; %cell size along the y-direction
dt=0.1; %time step size
% fluid properties
mu = 0.001; %viscosity of the fluid 
rho = 1; %density of the fluid
%iteration parameters and relaxation factors
omega_u = 0.7; %relaxation parameter for u momentum 
omega_v = 0.7; %relaxation parameter for v momentum
omega_p = 0.3; %relaxation parameter for pressure correction
beta=0.95; %blending factor for controlling interpolation in convection terms
outer_iterations = 50; %number of times we are going to iterate thru continuity, X and Y momentum
number_time_steps=1200;%number of time steps
iter_v = 10; %number of iterations for u and v momentum solvers
iter_p = 100; %number of iterations for pressure correction solver
%% Flow variables
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
dcu= zeros(mesh_y+2,mesh_x+2);%deferred correction upwind (first order interpolation)
dcc= zeros(mesh_y+2,mesh_x+2);%deferred corerction central (second order interpolation)
source = zeros(mesh_y+2,mesh_x+2); %source term at each stage to check mass conservation
%% BOUNDARY CONDITION SUBROUTINE %%
if (prob_type==1)
    %%CASE 1:fully developed channel flow with inlet velocity specified 
    u(mesh_y+2,:) = 0; %north boundary is the pipe wall
    u(1,:)=0; %south boundary is also the pipe wall
    u(:,1) = 1; %entry velocity is 1 m/s
    % %optional coding to add a jet in cross-flow 
%     v(1,mesh_y/2:mesh_y/2+5)=0.5;
    u(:,mesh_x+1)=1; %set the same velocity at outlet to conserve mass
elseif (prob_type==2)
    %%CASE 2: Pressure-driven Channel flow 
    u(mesh_y+2,:) = 0; %north boundary is the pipe wall
    u(1,:)=0; %south boundary is also the pipe wall
    u(:,1)=1; %inlet velocity specified
    pressure(:,mesh_x+1)=0; %fix outlet pressure
else
    %%CASE 3: Lid-driven cavity
    %initialize variables for the driven cavity problem
    u(mesh_y+2,:) = 1; %u velocity is fixed along the north boundary
    %other four sides are the walls hence the velocity is zero, this is our
    %default declaration
end
%% Solver using SIMPLE-algorithm
u_old = u;
v_old = v;tot=[];counter=[];residual_u=[];residual_v=[];
%time loop
for time=1:number_time_steps
    if (time==1)
        outer_iterations=6;
    else
        outer_iterations=6;
    end
    %Main (outer) loop
    for k = 1:outer_iterations
        %% X-Momentum subroutine
        %initialize coefficients for the u-momentum the centers 
        for j = 2:mesh_x
            for i=2:mesh_y+1
                %first initialize everything to first order accurate interpolation
                me=0.5*rho*dy*(u_old(i,j)+u_old(i,j+1));
                mw=0.5*rho*dy*(u_old(i,j)+u_old(i,j-1));
                mn=0.5*rho*dx*(v_old(i,j)+v_old(i,j+1));
                ms=0.5*rho*dx*(v_old(i-1,j)+v_old(i-1,j+1));
                ae(i,j) = max(-me,0); 
                aw(i,j) = max(mw,0); 
                an(i,j) = max(-mn,0);
                as(i,j) = max(ms,0); 
                apu(i,j)= ae(i,j)+aw(i,j)+as(i,j)+an(i,j)+me-mw+mn-ms;
                %compute the residual using first order
                dcu(i,j)=u(i,j)*apu(i,j)-(ae(i,j)*u(i,j+1)+aw(i,j)*u(i,j-1)+an(i,j)*u(i+1,j)+as(i,j)*u(i-1,j));
                %compute residual for second order
                dcc(i,j)=0.5*(me*(u(i,j)+u(i,j+1))-mw*(u(i,j)+u(i,j-1))+mn*(u(i,j)+u(i+1,j))-ms*(u(i,j)+u(i-1,j)));
                %add the viscous terms
                ae(i,j)=ae(i,j)+mu*dy/dx;
                aw(i,j)=aw(i,j)+mu*dy/dx;
                an(i,j)=an(i,j)+mu*dx/dy;
                as(i,j)=as(i,j)+mu*dx/dy;
            end
        end
        %correct the boundary values on the north and south boundary for the
        %X-momentum
        for j=2:mesh_x
            an(mesh_y+1,j) = max(-0.5*rho*dx*(v_old(mesh_y+1,j)+v_old(mesh_y+1,j+1)),0) + mu*dx/(dy/2);
            as(2,j) = max(0.5*rho*dx*(v_old(1,j)+v_old(1,j+1)),0) + mu*dx/(dy/2);
        end
        apu=ae+aw+as+an+mn-ms+me-mw+rho*dx*dy/dt;
        apu = apu/omega_u;
        if (iblock==1)
            apu(47:52,100:111)=1e30;
        end
        %iterate on the x-momentum equations
        for var = 1:iter_v
            for j = 2:mesh_x
                for i = 2:mesh_y+1
                    u(i,j) = (1-omega_u)*u_old(i,j) + (1/apu(i,j))*...
                        (ae(i,j)*u(i,j+1) + aw(i,j)*u(i,j-1) + an(i,j)*u(i+1,j)+ as(i,j)*u(i-1,j) -beta*(dcc(i,j)-dcu(i,j))+... 
                    dy*(pressure(i,j)-pressure(i,j+1))+rho*dx*dy*u_old(i,j)/dt);
                end
            end
            if (prob_type==1)
                u(:,mesh_x+1) = u(:,mesh_x); %du/dx=0 boundary condition at the exit
            end
        end
        if (prob_type==1)
            %ensure mass conservation each time for the exit
            in_mass_flow = 1; %as declared
            out_mass_flow =0;
            for i = 2:mesh_y
                out_mass_flow = out_mass_flow+ rho*dy*u(i,mesh_x+1);
            end
            u(:,mesh_x+1) = u(:,mesh_x+1)*in_mass_flow/out_mass_flow;
%             u(:,1) = u(:,mesh_x+1);%generate fully developed flow
        end
        %% RESIDUALS MONITOR DISPLAY
        if rem(k,1)==0
            res_u=0;
            for j = 2:mesh_x
                for i = 2:mesh_y+1
                    res_u = res_u+(u(i,j)-u_old(i,j)).^2;
                end
            end
            res_u=sqrt(res_u);
            counter=[counter,k];
            residual_u = [residual_u,res_u];
            figure(1);
            hold on;
            plot(counter,residual_u,'-b','LineWidth',2);grid on;
            title('Residuals of Mass Conservation with Iterations');
            xlabel('Iterations');
            ylabel('Residuals');        
        end
        %% Y-momentum subroutine
        %write the values of coefficients
        for j = 2:mesh_x+1
            for i = 2:mesh_y
                me=0.5*rho*dy*(u_old(i+1,j)+u_old(i,j));
                mw=0.5*rho*dy*(u_old(i+1,j-1)+u_old(i,j-1));
                mn=0.5*rho*dx*(v_old(i,j)+v_old(i+1,j));
                ms=0.5*rho*dx*(v_old(i,j)+v_old(i-1,j));
                ae(i,j) = max(-me,0)+ mu*dy/dx;
                aw(i,j) = max(mw,0) + mu*dy/dx;
                an(i,j) = max(-mn,0) + mu*dx/dy;
                as(i,j) = max(ms,0) + mu*dx/dy;
                apv(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+me-mw+mn-ms;
                %compute the residual using first order
                dcu(i,j)=v(i,j)*apv(i,j)-(ae(i,j)*v(i,j+1)+aw(i,j)*v(i,j-1)+an(i,j)*v(i+1,j)+as(i,j)*v(i-1,j));
                %compute residual for second order
                dcc(i,j)=0.5*(me*(v(i,j)+v(i,j+1))-mw*(v(i,j)+v(i,j-1))+mn*(v(i,j)+v(i+1,j))-ms*(v(i,j)+v(i-1,j)));
                %add the viscous terms
                ae(i,j)=ae(i,j)+mu*dy/dx;
                aw(i,j)=aw(i,j)+mu*dy/dx;
                an(i,j)=an(i,j)+mu*dx/dy;
                as(i,j)=as(i,j)+mu*dx/dy;            
            end
        end
        %overwrite east west boundary values
        for i = 2:mesh_y
            ae(i,mesh_y+1) = max(-0.5*rho*dy*(u_old(i+1,mesh_y+1)+u_old(i,mesh_y+1)),0) + mu*dy/(dx/2);
            aw(i,2) = max(0.5*rho*dy*(u_old(i+1,1)+u_old(i,1)),0) + mu*dy/(dx/2);
        end
        apv=ae+aw+as+an+me-mw+mn-ms+rho*dx*dy/dt;
        %iterative solver
        apv = apv/omega_v;
        if (iblock==1)
            apv(47:53,100:110)=1e30;
        end
        for var = 1:iter_v
            for j = 2:mesh_x+1
                for i = 2:mesh_y
                    v(i,j) = v_old(i,j)*(1-omega_v) + (1/apv(i,j))*...
                        (ae(i,j)*v(i,j+1) + aw(i,j)*v(i,j-1) + as(i,j)*v(i-1,j)+ an(i,j)*v(i+1,j) -beta*(dcc(i,j)-dcu(i,j))+...
                        dx*(pressure(i,j)-pressure(i+1,j))+rho*dx*dy*v_old(i,j)/dt);
                end
            end
            if (prob_type==1)
                v(:,mesh_x+1) = v(:,mesh_x); %dv/dx=0 boundary condition at the exit
            end
        end
        %% RESIDUALS MONITOR DISPLAY
        if rem(k,1)==0
            res_u=0;
            for j = 2:mesh_x
                for i = 2:mesh_y+1
                    res_u = res_u+(v(i,j)-v_old(i,j)).^2;
                end
            end
            res_u=sqrt(res_u);
            residual_v = [residual_v,res_u];
            plot(counter,residual_v,'-r','LineWidth',2);
            fprintf('\n');
        end
        %% Pressure subroutine
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
%         app(2,2)=1e30;%set reference cell value
        if (prob_type==2)
            app(:,mesh_x+1)=1e30;% FOR CASE 2: (Pressure driven flow):fix outlet pressure
        end
        p_prime = zeros(mesh_y+2,mesh_x+2);%initialize all pressure corrections to zero 
        source=0;%compute the source term
        for i = 2:mesh_y+1
            for j = 2:mesh_x+1
                source(i,j) = rho*dy*(u(i,j)-u(i,j-1)) + rho*dx*(v(i,j)-v(i-1,j));
            end
        end
        total = sqrt(sum(source.^2,'all'));
        if rem(k,1)==0
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
        for i =  2:mesh_y+1
            for j = 2:mesh_x
                u(i,j) = u(i,j) + (dy/apu(i,j))*(p_prime(i,j)-p_prime(i,j+1));
            end
        end
        %v-momentum
        for i = 2:mesh_y
            for j = 2:mesh_x+1
                v(i,j) = v(i,j)+(dx/apv(i,j))*(p_prime(i,j)-p_prime(i+1,j));
            end
        end
       % CASE 2 (fully developed channel flow)
        if (prob_type==2)
            for i=2:mesh_y+1
                u(i,mesh_x+1)=u(i,mesh_x)+(dx/dy)*(v(i-1,mesh_x+1)-v(i,mesh_x));
            end
        end
        %update velocities
        u_old = u;
        v_old = v;
        %% RESIDUAL MONITOR DISPLAY (for mass conservation)
        %recompute source term to check mass conservation
         for i = 2:mesh_y+1
            for j = 2:mesh_x+1
                source(i,j) = rho*dy*(u(i,j)-u(i,j-1)) + rho*dx*(v(i,j)-v(i-1,j));
            end
        end
        total = sqrt(sum(source.^2,'all'));
        fprintf('total residual = %0.4f\n',total);
%         monitor mass conservation
        if rem(k,1)==0
            tot = [tot,total];
            plot(counter,tot,'-g','LineWidth',2.5);legend('x-velocity','y-velocity','Mass-conservation');
        end
        hold off;
    fprintf('Iteration = %d of Navier Stokes\n',k)
    end
    %% END of the NAVIER STOKES SOLVER
    %% POST PROCESSING AND DISPLAY PLOTS
    %quiver-plot: a type of plot where the velocities are represented with
    %arrows give (x,y) components at a given coordinate or (u,v) at (x,y)
    %hence interpolate staggered mesh-velocities to the corners of the main
    %control volume after computing X, Y vectors as cell corner coordinates

    %X-coordinates computation
    X = [];
    for i = 1:mesh_x+1
            x = (i-1)*dx;
            X = [X,x];
    end

    %Y-coordinates computation
    Y = [];
    for j = 1:mesh_y+1
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
    U(mesh_y+1,1) = 0;
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
    %blocked cells
    %quiverplot for velocities
    if (iblock==1)
        U(47:52,100:111)=0;
        V(47:63,100:110)=0;
    end
%     figure(2);
%     quiver(X,Y,U,V,'Color','k','LineWidth',2);
    %% pressure terms computation
    P = zeros(mesh_y+1,mesh_x+1);
    %interior points
    for i = 2:mesh_y
        for j=2:mesh_x
            P(i,j) = 0.25*(pressure(i,j) + pressure(i+1,j) + pressure(i,j+1) + pressure(i+1,j+1));
        end
    end
    %boundaries except corners
%     %west boundary
%     P(2:mesh_y,1) = 0.5*(pressure(2:mesh_y,1)+ pressure(3:mesh_y+1,1));
%     %east 
%     P(2:mesh_y,mesh_x+1) = 0.5*(pressure(2:mesh_y,mesh_x+1) + pressure(3:mesh_y+1,mesh_x+1));
%     %north
%     P(mesh_y+1,2:mesh_x) = 0.5*(pressure(mesh_y+1,2:mesh_x) + pressure(mesh_y+1,3:mesh_x+1));
%     %south 
%     P(1,2:mesh_x) = 0.5*(pressure(2,2:mesh_x) + pressure(2,3:mesh_x+1));
%     %corner points
%     %southwest
%     P(1,1) = pressure(2,2);
%     %southeast
%     P(1,mesh_x+1) = pressure(2,mesh_x+1);
%     %northeast
%     P(mesh_y+1,mesh_x+1) = pressure(mesh_y+1,mesh_x+1);
%     %northwest
%     P(mesh_y+1,1) = pressure(mesh_y+1,2);
%     figure(3);
%     contourf(X,Y,U,'ShowText','off','LineColor','none','LevelStep',0.01);
%     colormap(jet);
%     title('X-velocity Contours for Channel flow');
%     figure(4);
%     contourf(X,Y,V,'ShowText','off','LineColor','none','LevelStep',0.01);
%     colormap(jet);
%     title('Y-velocity Contours for Channel flow');
    % figure(5);
    % contourf(X,Y,P,'ShowText','off','LineColor','none','LevelStep',1);
    % colormap(jet);
    % title('Pressure Contours for Channel flow');
%     if (rem(time,10)==0)
%         figure(2);
%         hold on;
%         contourf(X,Y,sqrt(U.^2+V.^2),'ShowText','off','LineColor','none','LevelStep',0.01);
%         colormap(jet);
%         title('Velocity magnitude Contours for Channel flow')
%     end
    if (rem(time,1)==0)
        fprintf('time = %d done\n', time);
        figure(2);
        hold on;
        contourf(X,Y,sqrt(U.^2+V.^2),'LevelStep',0.05,'LineColor','none');
        colormap(jet);set(gcf,'position',[10,10,5000,500]);
%         quiver(X,Y,U,V,5,'Color','k','LineWidth',0.5);
        title('Velocity-magnitude contours-Vortex shedding around a square');
    end
end
fprintf('done with computations!');

    

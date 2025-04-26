close all
clear all
clc

%%%%%%%%%%%%%%% My convention: Head -> 0, Toe -> 1;
% All material, mission and boundary properties needed for project three
% Material properties provided in header, along with plots for ease of
% visualization

% Assume inner core of robot is perfectly insulated, meaning heat cannot
% escape inner most layer

% Mission times: 1 hour, 3 hours and 7 hours

% Thermal protection material properties
% Thermal conductivity: 0.2 W/m*K
% Specific heat capacity: 1200 J/Kg*K
% Glass transition temperature: 1200K
% Density: 160 kg/m^3

% Outer surface: carbon-fiber
% Thermal conductivity: 500 W/m*K
% Specific heat capacity: 700 J/Kg*K
% Glass transition temperature: 350K

% Middle: glue
% Thermal conductivity: 200 W/m*K
% Specific heat capacity: 900 J/Kg*K
% Glass transition temperature: 400K
% Constant thickness: 0.1 cm

% Inner most: steel
% Thermal conductivity: 100 W/m*K
% Specific heat capacity: 500 J/Kg*K
% Max temperature: 800K

T = 3600;   % total time = 1hr = 3600s
C_thickness = CarbonThickness(false);
G_thickness = GlueThickness(false);
Fe_thickness = SteelThickness(false);
total_thickness = C_thickness + G_thickness + Fe_thickness;
u0 = InitialTemp(false);
dx = 0.0001;
dt = 1;
z = 1;
x = 0:dx:0.18;  % Spatial grid
Nx = length(x);
u = 300*ones(Nx,1);
u(1) = u0(end);
figure
for i = 1:dt:T
    [u_new, finish] = TempAfterCarbonLayer(Fe_thickness, z, dx, dt, u, i, 100);
    u = u_new;
    pause(0.01);
    if(finish)
        break
    end
end
%% Carbon layer thickness profile:
function y = CarbonThickness(Plot)
    f = 1;                  % Frequency in Hz
    cycles = 1;             % Number of cycles
    Fs = 1000;              % Sampling frequency in Hz
    T = cycles / f;         % Duration to complete 1.5 cycles
    A = 1.5;                % Amplitude in cm

    % Time vector
    t = 0:1/Fs:T;

    % Thickness 
    y = abs(A * sin(2*pi*f*t))+ 0.1;

    if(Plot)
        % Plotting
        figure;
        plot(t, y, 'b-', 'LineWidth', 1.5);
        xlabel('Non-dimensional distance l/L measured from robot head to toe');
        ylabel('Layer Thickness (cm)');
        title('Carbon fiber thickness starting from head to toe, Robot Length L = 2.5 m');
        grid on;
    end
end

%% Glue thickness profile:
function y = GlueThickness(Plot)
    % Domain
    x = linspace(0, 1, 1001);

    % Let’s determine A and B empirically
    A = 0.1;     % Controls the curve steepness
    B = 20;        % Controls how fast it levels off
    C = 0.01;     % Starting value

    % Thickness function
    y = A * log(B * x + 1) + C;

    if(Plot)
        % Plot
        figure;
        plot(x, y, 'r-', 'LineWidth', 2);
        xlabel('Non-dimensional distance l/L measured from robot head to toe');
        ylabel('Layer Thickness (cm)');
        title('Glue thickness starting from head to toe, Robot Length L = 2.5 m');
        grid on;
        grid on;
    end
end
%% Steel profile:
function y_shifted = SteelThickness(Plot)
    % Parameters
    f = 5;                   % Frequency in Hz
    A = 5;                   % Amplitude in cm
    T = 1/f;                 % Period = 5 seconds
    Fs = 1000;               % Sampling frequency in Hz
    duration = 1;            % Duration of signal in seconds (e.g., 3 full cycles)

    % Time vector
    t = 0:1/Fs:duration;

    % Generate sawtooth wave (range: -1 to 1), scale to amplitude
    y = A * sawtooth(2*pi*f*t);   % Full swing from -5 cm to +5 cm

    % Thickness function:
    y_shifted = (A/2) * (sawtooth(2*pi*f*t) + 1) + 0.1;  % Range: 0.1 to 5.1 cm

    if(Plot)
        % Plot
        figure;
        plot(t, y_shifted, 'm', 'LineWidth', 1.5);
        xlabel('Non-dimensional distance l/L measured from robot head to toe');
        ylabel('Layer Thickness (cm)');
        title('Steel thickness starting from head to toe, Robot Length L = 2.5 m');
        grid on;
    end
end

%% Temperature profile:
function T = InitialTemp(Plot)
    % Domain
    x = linspace(0, 1, 1001);

    % Empirically chosen parameters to fit T(0) = 900, T(1) ≈ 700
    A = -100;
    B = 8;
    C = 900;
    
    x_mean = mean(x);
    x = 2*x_mean-x;
    T = A * log(B*x + 1) + C;
    if(Plot)
        % Plot
        figure;
        plot(x, T, 'r', 'LineWidth', 2);
        xlabel('Non-dimensional distance l/L measured from robot toe to head');
        ylabel('Layer Thickness (cm)');
        title('Exhaust Gas Temperature profile starting from toe to head, Robot Length L = 2.5 m');
        grid on;
    end
end

% %% Temp dist. after carbon layer
% function [temp, finish] = TempAfterCarbonLayer(C_thickness, z, dx, dt, u, idx, stopTime)
%     finish = false;
%     
%     k = 0.2;           % Thermal conductivity [W/m-K]
%     rho = 160;        % Density [kg/m^3]
%     cp = 1200;          % Specific heat [J/kg-K]
%     alpha = k / (rho * cp); 
% 
%     % Spatial domain
%     z = 1000*z + 1; 
%     x = 0:dx:0.2;  % Spatial grid
%     Nx = length(x);
%     % Lambda
%     lambda = alpha * dt / dx^2;
% 
%     % Matrix coefficients
%     a_diag = (1 + 2 * lambda) * ones(Nx, 1);
%     b_off = -lambda * ones(Nx-1, 1);
%    
%     % Storage for results
%     temp = [];
%     temp_before = u;
%     u_new = thomas_algorithm(b_off, a_diag, b_off, u);
% %     A = diag(a_diag) + diag(b_off, -1) + diag(b_off, 1);
% %     u_new = A\u;
%     u_new(1) = 900; % Keep the hot boundary fixed
%     temp = u_new;
%     % Plot final temperature profile
% %     if(idx == stopTime)
%     plot(x, temp);
%     xlabel('Depth [m]');
%     ylabel('Temperature [K]');
%     title(sprintf("1D Heat Conduction after %d sample", int32(idx/dt+1)));
% %     end
% 
%     if(max(abs(temp - temp_before)) < 0.05)
%         finish = true;
%     end
% 
% end
% 

%% new implementation
function [u_new, finish] = TempAfterCarbonLayer(C_thickness, z, dx, dt, u_old, idx, stopTime)
    finish = false;
    
    % material properties (carbon fiber)
    k   = 0.2;    % W/(m·K)
    rho = 160;    % kg/m^3
    cp  = 1200;   % J/(kg·K)
    alpha = k/(rho*cp);

    x = 0:dx:0.18;
    Nx = length(x);
    
    lambda = alpha * dt / dx^2;

    % We solve only for u(2) through u(Nx) — u(1) is fixed at 900 K.
    n = Nx - 1;                   % number of unknowns per time step
    a = -lambda * ones(n-1,1);    % sub-diagonal
    b =  (1+2*lambda) * ones(n,1);% main diagonal
    c = -lambda * ones(n-1,1);    % super-diagonal
    b(end) = 1 + lambda;          % modify last row for insulated BC
    d = u_old(2:end);                % RHS = u^n for i=2:N
    
    % Adjust first equation for Dirichlet at left:
    %  -λ·u(1) + (1+2λ)u(2) - λu(3) = u_old(2)
    % ⇒ RHS(1) += λ·900
    d(1) = d(1) + lambda * 900;

    % Modify last equation for insulated right BC:
    % Instead of -λ u(N-1) + (1+2λ) u(N) - λ u(N+1) = u_old(N)
    % we do   -λ u(N-1) + (1+λ) u(N) = u_old(N)
    b(end) = 1 + lambda;

    % Solve tridiagonal system for u(2:Nx)
    u_unknown = thomas_algorithm(a, b, c, d);

    % Reconstruct full solution:
    u_new = zeros(Nx,1);
    u_new(1)    = 900;        % left Dirichlet
    u_new(2:end)= u_unknown;  % interior & right

    % Plot
    plot(x, u_new, 'LineWidth',1.5);
    xlabel('Depth [m]');
    ylabel('Temperature [K]');
    title(sprintf('1D BTCS after t = %.1f s', idx));
    drawnow;

    % Check convergence
    if max(abs(u_new - u_old)) < 0.05
        finish = true;
    end
end

%% Temp dist. after glue layer
function temp = TempAfterGlueLayer()
    
end

%% Temp dist. after steel layer
function temp = TempAfterSteelLayer()
    
end

%% Thomas Algo
function x = thomas_algorithm(a, b, c, d)
    % Inputs:
    % a -> sub-diagonal 
    % b -> main diagonal 
    % c -> super-diagonal 
    % d -> right-hand side vector 

    n = length(b);
    c_prime = zeros(n-1, 1);
    d_prime = zeros(n, 1);
    
    % Forward sweep
    c_prime(1) = c(1) / b(1);
    d_prime(1) = d(1) / b(1);

    for i = 2:n-1
        denom = b(i) - a(i-1) * c_prime(i-1);
        c_prime(i) = c(i) / denom;
        d_prime(i) = (d(i) - a(i-1) * d_prime(i-1)) / denom;
    end

    d_prime(n) = (d(n) - a(n-1) * d_prime(n-1)) / (b(n) - a(n-1) * c_prime(n-1));

    % Back substitution
    x = zeros(n, 1);
    x(n) = d_prime(n);
    
    for i = n-1:-1:1
        x(i) = d_prime(i) - c_prime(i) * x(i+1);
    end
end

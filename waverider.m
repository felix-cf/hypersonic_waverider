%%
clear

vRider=0.0073;
m=vRider*2500;

L=@(M, alt) force_m_alt_lift(M, alt); 
D=@(M, alt) force_m_alt_drag(M, alt);
setGlobalLD(1)
setGlobaldeltalast(1)
last = 1;
deltalast = 1;
gam=1.4;
R=287;
g=9.8;


[t,y] = ode45(@myodefun2,[0 1000] ,[0, 10*sqrt(atmosisa(40000)*R*gam), 40000, 0],[], m, D, L, R, gam, g);

%%


figure
plot(t, y(:, 1), 'b')
hold on
plot(t, y(:, 3), 'r')

%%
finalval = 491;
figure
plot(y(1:finalval,1), y(1:finalval,3), 'b')
title("Altitude vs Distance")
xlabel("Distance (m)")
ylabel("Altitude (m)")
Mach = sqrt(y(1:finalval, 2).^2  + y(1:finalval, 4).^2) ./ sqrt(atmosisa(y(1:finalval,3))*R*gam);

figure
plot(y(1:finalval,1), Mach, 'r')
title("Mach number vs Distance")
xlabel("Distance (m)")
ylabel("Mach")
%%
tiledlayout(1, 2)
nexttile
res = getGlobalLoverD;
plot(y(1:finalval, 1), res(1, 1:finalval), 'o', 'DisplayName','Values')
xlabel("Distance (m)")
ylabel("L/D")
title("L/D vs. Distance")
nexttile
plot(Mach(1:5:finalval, 1), res(1, 1:5:finalval), 'o', 'DisplayName','Values')
xlabel("Mach")
ylabel("L/D")
title("L/D vs. Mach")

%%
tiledlayout(1,2)
nexttile
res = getGlobaldelta;
plot(y(1:finalval, 1), (180/pi)*res(1:finalval), 'o')
xlabel("Distance (m)")
ylabel("Delta (degrees)")
title("Deflection angle vs. Distance")
nexttile
plot(y(1:finalval, 1), 90 - (180/pi)*(res(1:finalval)), 'o')
xlabel("Distance (m)")
ylabel("Resultant Force Angle (degrees)")
title("Resultant Force Angle vs. Distance")
%%

% Generating the pressure distributions of different states during the
% trajectory
force_m_alt(10, 40000)
force_m_alt(7, 10000)
%%

% Functions for taking data for L/D

% Set and get the LD index 
function setGlobalLD(val)
    global last
    last = val;
end
function r = getGlobalLD
    global last;
    r = last;
end

% Set and get the LD matrix
function setGlobalLoverD(val)
    global LoverD;
    mylast = getGlobalLD;
    LoverD(mylast) = val;
    setGlobalLD(mylast + 1);
end

function r = getGlobalLoverD
    global LoverD;
    r = LoverD;
end

% set and get delta index
function setGlobaldeltalast(val)
    global deltalast
    deltalast = val;
end

function r = getGlobaldeltalast
    global deltalast;
    r = deltalast;
end

% set and get delta matrix
function setGlobaldelta(d)
    global deltas;
    mylast = getGlobaldeltalast;
    deltas(mylast) = d;
    setGlobaldeltalast(mylast + 1);
end

function r = getGlobaldelta
    global deltas;
    r = deltas;
end

% y = [x dx y dy]
function ydot = myodefun2(t,y,m, D, L, R, gam, g)
    myL = L(sqrt(y(2)^2+y(4)^2)/sqrt(atmosisa(y(3))*R*gam), y(3));
    myD = D(sqrt(y(2)^2+y(4)^2)/sqrt(atmosisa(y(3))*R*gam), y(3));
    setGlobalLoverD(myL/myD);
    [f, delta] = force_m_alt(sqrt(y(2)^2+y(4)^2)/sqrt(atmosisa(y(3))*R*gam), y(3));
    setGlobaldelta(delta);
    ydot(1,1) = y(2);
    ydot(2,1) = -myD/m;
    ydot(3,1) = y(4);
    ydot(4,1) = myL/m-g;
end
%%
function f2 = force_m_alt_lift(M_inf,alt)
    [fint, delta] = force_m_alt(M_inf, alt);
    f2 = fint * cos(delta);
end
function f3 = force_m_alt_drag(M_inf,alt)
    [fint, delta] = force_m_alt(M_inf, alt);
    f3 = fint * sin(delta);
end 


function [f, delta] = force_m_alt(M_inf, alt)
    
    % Step 1
    
    % Solving for the deflection angle and Mach angle immediately after the 
    % shock, given an arbitrary shock wave angle
    
    % 20 degree shock wave angle
    theta_s = 20 * (pi/180);
    [T_1, a_1, P_1, rho_1] = atmosisa(alt);
    
    delta = atan( ...
        2*cot(theta_s)*( ...
            (M_inf^2*sin(theta_s)^2 - 1) / ...
            (M_inf^2 * (1.4 + cos(2*theta_s)) + 2) ...
            ) ...
        );
    
    % Angle in degrees of deflection
    delta * (180/pi);
    
    % M2 as a function of M of stream, shock angle, and deflection angle
    M2 = @(M, beta, theta) ( (2 + .4*(M*sin(beta))^2 ) / (sin(beta-theta)^2*(2*1.4*(M*sin(beta))^2 - .4) ) )^(.5);
    
    
    % Confirmed post-shock M 
    M_2 = M2(M_inf, theta_s, delta);
    
    %
    % Step 2
    % Find the radial and normal components of flow velocity
    M_1n = M_inf * sin(theta_s);
    T_2 = T_1 * (1 + ((2*1.4)/(1.4+1))*(M_1n^2-1) ) * ((2+(1.4-1)*M_1n^2)/((1.4+1)*M_1n^2));
    P_2 = P_1 * (1 + ((2*1.4)/(1.4+1))*(M_1n^2-1));
    a_2 = (1.4*287*T_2)^(.5);
    V = M_2 * a_2;
    V_max = sqrt((2/(1.4-1))*a_2^2 + V^2);
    V_prime = (2/((1.4-1)*M_2^2) + 1)^(-1/2);
    
    % Finding the radial and normal components of velocity 
    V_r_prime = V_prime * cos(theta_s - delta);
    V_theta_prime = -V_prime * sin(theta_s - delta);
    
    %
    % Step 3
    % Solve numerically for V_r_prime with the above value of V_r_prime as a
    % boundary condition. We will get a flow field of V's across theta
    
    % v = [v_r_prime v_theta_prime]
    [theta, v] = ode45(@myodefun,[20*(pi/180) 10*(pi/180)], [V_r_prime; V_theta_prime]);
    %figure
    %plot(theta*(180/pi), v(:,1), 'DisplayName', "V_r'")
    %hold on
    %plot(theta*(180/pi), v(:,2), 'DisplayName', "V_{theta}'")
    %title("Numerical Solution of V_r' and V_{theta}'")
    %legend()
    %xlabel("Theta")
    %ylabel("Velocity")
    
    %
    % Step 4
    % Find the zero of v_theta_prime
    % By inspection of previous graph, it's about 17.45 degrees. 
    % Verified here:
    theta_c = theta(find(abs(v(:,2)) < 0.005), :);
    [theta, v] = ode45(@myodefun,[20*(pi/180) theta_c(1)], [V_r_prime; V_theta_prime]);
    
    %figure
    %plot(theta*(180/pi), v(:,1), 'DisplayName', "V_r'")
    %hold on
    %plot(theta*(180/pi), v(:,2), 'DisplayName', "V_{theta}'")
    %legend()
    %xlabel("Theta")
    %ylabel("Velocity")
    %title("Numerical Solution of V_r' and V_{theta}'")
    theta_c = theta(41,:);
    % Approximate theta_c
    
    % V_r_prime at 22.36 is V at surface of cone. Since V_theta_prime is 0 at
    % this point, V' here is all the V_r_prime component
    V_surface = v(41, 1); % Last row, first column
    M_surface =  sqrt((2/(1.4-1)) * ((V_surface^2)/(1-V_surface^2)));
    
    %
    % Step 5
    % Velocity flow field
    V_field = sqrt(v(:,1).^2 + v(:,2).^2);
    M_field = sqrt((2/(1.4-1)) * ((V_field.^2)/(1-V_field.^2)));
    M_field = M_field(:, 41);
    
    %semilogy(theta*(180/pi), M_field)
    %set(gca,'XDir', 'reverse')
    %title("Mach Field from Shock to Cone angle")
    %xlabel("Theta")
    %ylabel("Mach")

    %
    cone = @(x, y) sqrt((tan(theta_s)*x)^2 - y^2);
    cone2 = @(x, y) -sqrt((tan(theta_s)*x)^2 - y^2);
    
    span = [0 10 -5 5];
    %figure
    %fsurf(cone, span)
    %hold on
    %fsurf(cone2, span)
    
    %figure
    dist = 5; % Distance downstream of cone peak
    length = 1; % Length of waverider
    offset = -tan(theta_s)*dist;
    u_ang = delta + pi/16; % angle of upper surf.
    bottomsurface = @(x, y) tan(theta_s)*(x-dist)*(sqrt((tan(delta)*(x-dist))^2 - (y)^2) + tan(theta_s)*offset) + offset;
    uppersurface = @(x, y) tan(theta_s)*(x-dist)*(sqrt((tan(u_ang)*(x-dist))^2 - (y/.55)^2) + tan(theta_s)*offset) + offset;
    offset_wave = -tan(delta)*length; 
    span_wave = [dist dist+length offset_wave -offset_wave];
    %fsurf(bottomsurface, span_wave, 'b', 'EdgeColor', 'none', 'LineStyle', 'none', 'Marker','o')
    %hold on
    %fsurf(uppersurface, span_wave, 'r', 'EdgeColor', 'none', 'LineStyle', 'none', 'Marker','o')
    
    
    % Surface area and pressures are idealized to speed up computation
    % Assume that the surface area of the wave rider is approximately equal
    % half of a cone that would define the surface of the waverider.
    a = tan(theta_c);
    surface_area = (1/2) * pi * a * (a + sqrt(length^2 + a^2));
    % We approximate the force due to pressure as a fraction of the
    % pressure generated by the bottom surface, since the upper surface has
    % a similar geometry and is physically very close to the bottom surface.
    % This reduces the computational load for our algorithm immensely, and
    % is a relatively sane simplifying assumption to make, as the following
    % are reasonable assumptions:

    % Force = SA_lower * P_lower - SA_upper * P_upper
    % P_lower ~  P_upper
    % SA_upper < SA_lower
    % SA_upper ~ c*SA_lower where c is slightly less than 1
    % Force ~ P * SA_lower * (1 - c)
    % Force ~ P * SA_lower * C` where C` << 1
    pres1 = pressure(theta_s, theta, M_field, P_2, M_2);
    f = pres1*surface_area*(1/500);

    
    % For generating pressure distribution data at different stages of the
    % trajectory

    %inc = theta_c;
    %for ii=1:1:41
        %int = abs(theta_c - theta_s);
        %inc = inc + int/41;
        %pressuredist(ii) = pressure(inc, theta, M_field, P_2, M_2); 
    %end
    
    %figure
    %plot(theta.*(180/pi), pressuredist)
    %title("Pressure Distribution at Mach " + M_inf + " and Altitude " + alt)
    %xlabel("Theta Angle (degrees)")
    %ylabel("Pressure (Pa)")

end

%%
%syms x y z
%bottomsurface = tan(theta_s)*(x-dist)*(sqrt((tan(delta)*(x-dist))^2 - (y)^2) + tan(theta_s)*offset) + offset;
%res = vpaintegral(vpaintegral(pressure(x, y, z, theta, M_field, P_2, M_2)*bottomsurface, y, -tan(delta), tan(delta)),x, 5, 6)


%%


%syms v vdot vddot gamma theta
%solve(((gamma-1)/2)*(1-v^2-vdot^2)*(2*v + vdot*cot(theta) + vddot) - vdot*(v*vdot + vdot*vddot) == 0, vddot)
function dvdt = myodefun(theta, v)
    dvdt(1) = v(2);
    dvdt(2) = -(v(1)*v(2)^2 + (1.4/2 - 1/2)*(2*v(1) + v(2)*cot(theta))*(v(1)^2 + v(2)^2 - 1))/((1.4/2 - 1/2)*(v(1)^2 + v(2)^2 - 1) + v(2)^2);
    dvdt = dvdt';
end

% theta wrt x,y,z
function p = pressure(theta_fun, theta_vals, M_vals, P_2, M_2)
    findval = find(abs(theta_vals - theta_fun) < 0.001);
    M_theta = M_vals(findval(1), :);
    p =P_2 * ( (2 + .4*M_2.^2) / (2 + .4*M_theta.^2) ).^(1.4/.4);
end


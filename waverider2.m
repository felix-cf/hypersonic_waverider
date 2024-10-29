% 20 degree shock wave angle
theta_s = 20 * (pi/180);
M_inf = 10;
[T_1, a_1, P_1, rho_1] = atmosisa(40000);

delta = atan( ...
    2*cot(theta_s)*( ...
        (M_inf^2*sin(theta_s)^2 - 1) / ...
        (M_inf^2 * (1.4 + cos(2*theta_s)) + 2) ...
        ) ...
    );

theta_s = 20 * (pi/180);
dist = 5; % Distance downstream of cone peak
length = 1; % Length of waverider
offset = -tan(theta_s)*dist;
u_ang = delta + pi/16; % angle of upper surf. 15 degrees between surfaces.
bottomsurface = @(x, y) tan(theta_s)*(x-dist)*(sqrt((tan(delta)*(x-dist))^2 - (y)^2) + tan(theta_s)*offset) + offset;
uppersurface = @(x, y) tan(theta_s)*(x-dist)*(sqrt((tan(u_ang)*(x-dist))^2 - (y/.46)^2) + tan(theta_s)*offset) + offset;
span_wave = [dist dist+length offset -offset];


% Here we perform an offline computation of the bounds for the region in
% which we are calculating the volume to reduce the computational load of
% the online guessing portion of the Monte Carlo method
yb=[];
c=1;

for x=5.2:0.01:6
    syms ybound 
    ybound=vpasolve(bottomsurface(x, ybound)-uppersurface(x, ybound), ybound, 0.5);
    yb(c)=ybound;
    c=c+1;
end

%%

% Monte Carlo method part 2, online guessing of points within volume 

%number of guesses
total_guesses=100000;
%number of guesses that land inside the wave rider
inside=0;
counti=1;
counto=1;
count=1;

for i=1:1:total_guesses
    
    %generate random points
    x=rand*(6-5.2) +5.2; %random number from 5.2 to 6
    y=rand-0.5; %random value from -0.5 to 0.5
    z=rand*0.2-2;%random value from -1.8 to -2

    if (bottomsurface(x, y)<z) && (z<uppersurface(x, y)) && (ybound*-1<y) && (ybound>y)
        inside=inside+1;
        xi(counti)=x;
        yi(counti)=y;
        zi(counti)=z;
        counti=counti+1;
    else 
        xo(counto)=x;
        yo(counto)=y;
        zo(counto)=z;
        counto=counto+1;
    end
    count = count + 1;

end
%%
clf
vTotal = 1*0.2*0.8;
plot3(xi, yi, zi, 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
hold on
plot3(xo, yo, zo, 'o', 'MarkerSize', 0.1, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
vRider=inside/total_guesses * vTotal

xlim([5.2 6])
zlim([-2 -1.8])

%%
max_m = 9;
for ii=1:max_m
    total_guesses = 10^ii;
    inside = 0;
    for j=1:total_guesses
        x = rand-1;
        y = rand-1;
        z = rand-1;

        if sqrt(x^2 + y^2 + z^2) < 1
            inside = inside+1;
        end
    end
    monte(ii) = (inside/total_guesses) * 2*2*2;
end
%%
for ii=1:max_m
    error(ii) = (abs((4/3) * pi - monte(ii)) / (4/3 * pi)) * 100;
end

approx = @(x) error(1) * x^(-4.5);
figure
semilogy(1:max_m, error, 'DisplayName', "Experimental Error Values")
hold on
fplot(approx, [1 9], 'DisplayName', "Approximate Convergence, x^{-4.5}")
title("Convergence of Error of Monte Carlo Method")
xlabel("Number of Guesses (10^x)")
ylabel("Relative Error (%)")
legend()
% README:
% 
% %% Set Environment: Diesen Abschnitt einmal ausführen um alle Gleichungen
% berechnen zu lassen. Hier kann die Periode, also die Zeit, die gebraucht
% wird einmal über der kompletten Kontur durch zulaufen, über die Variable
% 'period' geändert werden. Die Drehzahl wird automatisch berechnet. Die
% Abtastrate, also die Feinheit kann über die Schrittzahl-Variable 'steps'
% geändert werden. Je höher die Variable gewählt wird, desto feiner werden
% die Berechnungen, desto langsamer jedoch werden die Animationen.
% Wenn etwas geändert wird muss dieser Abschnitt nochmal ausgeführt werden.
% 
% %% Plot Images: Diesen Abschnitt ausführen wenn die Plots erstellt werden
% sollen.
% 
% %% Plot Animation: Diesen Abschnitt ausführen um die Animation laufen zu
% lassen. Hier kann die Skalierung geändert werden über 'scale_v' und
% 'scale_a'.
% 
% Wenn der Wert einer bestimmten kinematischen Gleichung zu einem
% bestimmten Zeitpunkt bestimmt werden soll: Im Command Window
% Beispielsweise: a_circular(find(t==0.370)) eingeben und ausführen, wobei
% 0.370 durch den gewünschten Zeitpunkt in Sekunden ersetzt wird.
% 
% Wenn der erreichte Winkel zu einem bestimmten Zeitpunkt bestimmt werden
% soll: Im Command Window: t*w_0 eingeben und ausführen, wobei t durch den
% gewünschten Zeitpunkt in Sekunden ersetzt wird

%% Set Environment
% set 'Periode' period (s), 'Schrittanzahl' steps (1), 'Abtastrate' sampleRate (s), 'Zeitvariable' time variable t (s), 'Winkelgeschw.' angular velocity w_0 (1/s)
period = 0.4;
steps = 5000;
sampleRate = period/steps;         %adjust in this line!
t = 0.0:sampleRate:period;
w_0 = (2*pi)/period;

% Radius of the big circle and the small circle
r_small = 4.84070738;
r_big = 8.10317546;

% Centers of all circles
x_1 = 25.74;
y_1 = 0.0;
x_2 = 18.09718091;
y_2 = 10.44841227;
x_3 = 12.86876268;
y_3 = 22.28935079;
x_4 = 0.0;
y_4 = 20.89681945;
x_5 = -12.86876268;
y_5 = 22.28935079;
x_6 = -18.09717651;
y_6 = 10.44840973;
x_7 = -25.73752542;
y_7 = 0.0;
x_8 = -18.09717681;
y_8 = -10.4484099;
x_9 = -12.8687635;
y_9 = -22.28935221;
x_10 = 0.0;
y_10 = -20.89681986;
x_11 = 12.86876362;
y_11 = -22.28935242;
x_12 = 18.09717686;
y_12 = -10.44840993;
x_13 = 25.73752722;
y_13 = 0.0;

% Equations of the upper and lower side of a circle 
r_unten = @(x, y, R) x.*cos(t*w_0) - (x.^2.*cos(t*w_0).^2 - y.^2.*cos(t*w_0).^2 + R.^2 - x.^2 + x.*y.*sin(2.*t*w_0)).^(1/2) + y.*sin(t*w_0);
r_oben = @(x, y, R)  x.*cos(t*w_0) + (x.^2.*cos(t*w_0).^2 - y.^2.*cos(t*w_0).^2 + R.^2 - x.^2 + x.*y.*sin(2.*t*w_0)).^(1/2) + y.*sin(t*w_0);

% Equation for setting sequences
seq = @(t_1, t_2) heaviside(t*w_0-t_1).*heaviside(-t*w_0+t_2);

% Assembled contour of all circle-sequences
r = r_unten(x_1, y_1, r_small).*seq(0.0, 0.169147) + r_oben(x_2, y_2, r_big).*seq(0.169147, 0.87805055) + r_unten(x_3, y_3, r_small).*seq(0.87805055, 1.21634458) + r_oben(x_4, y_4, r_big).*seq(1.21634458, 1.92524808) + r_unten(x_5, y_5, r_small).*seq(1.92524808, 2.26354213) + r_oben(x_6, y_6, r_big).*seq(2.26354213, 2.97244563) + r_unten(x_7, y_7, r_small).*seq(2.97244563, 3.31073968) + r_oben(x_8, y_8, r_big).*seq(3.31073968, 4.01964315) + r_unten(x_9, y_9, r_small).*seq(4.019643, 4.35793726) + r_oben(x_10, y_10, r_big).*seq(4.35793726, 5.06684069) + r_unten(x_11, y_11, r_small).*seq(5.06684069, 5.40513482) + r_oben(x_12, y_12, r_big).*seq(5.40513482, 6.11403824) + r_unten(x_13, y_13, r_small).*seq(6.11403824, 2*pi);
r(1) = r(2);
r(length(r)) = r(length(r) - 1);

% Convert mm to m
r = r/1000;

% Calculate derivatives
r_prime = diff(r, 1)/sampleRate;
r_doublePrime = diff(r_prime, 1)/sampleRate;

% Get all matrices in same size
r_prime = [r_prime r_prime(length(r_prime))];
r_doublePrime = [r_doublePrime r_doublePrime(length(r_doublePrime)) r_doublePrime(length(r_doublePrime))];

% Calculation of actual kinematic equations (according to http://www.tu-berlin.de/uploads/media/Vorlesungsnotizen_MeII.pdf)
v_radial = r_prime;
v_circular = r*w_0;
v_absolute = sqrt(v_radial.^2 + v_circular.^2);

a_radial = r_doublePrime - (r.*(w_0^2));
a_circular = 2*r_prime*w_0;
a_absolute = sqrt(a_radial.^2 + a_circular.^2);

%% Plot Images in Polar Coordinates

% Plot Contour
figure
polarplot(t*w_0, r)
title ('Kontur r(t)')

% Plot velocity
figure
subplot(1,3,1)
polarplot(t*w_0, v_radial)
rlim([-1.6, .8]);
title('Radiale Geschwindigkeit')

subplot(1,3,2)
polarplot(t*w_0, v_circular)
title('Zirkulare Geschwindigkeit')

subplot(1,3,3)
polarplot(t*w_0, v_absolute)
title('Absolute Geschwindigkeit')

% Plot acceleration
figure
subplot(1,3,1)
polarplot(t*w_0, a_radial)
rlim([-350, 350]);
title('Radiale Beschleunigung')

subplot(1,3,2)
polarplot(t*w_0, a_circular)
rlim([-50, 25]);
title('Zirkulare Beschleunigung')

subplot(1,3,3)
polarplot(t*w_0, a_absolute)
rlim([-300, 400]);
title('Absolute Beschleunigung')

%% Plot Images in Cartesian Coordinates

% Plot Contour
figure
plot(t*w_0, r)
title ('Kontur r(t)')

% Plot velocity
figure
subplot(3,1,1)
plot(t*w_0, v_radial)
% rlim([-1.6, .8]);
title('Radiale Geschwindigkeit')

subplot(3,1,2)
plot(t*w_0, v_circular)
title('Zirkulare Geschwindigkeit')

subplot(3,1,3)
plot(t*w_0, v_absolute)
title('Absolute Geschwindigkeit')

% Plot acceleration
figure
subplot(3,1,1)
plot(t*w_0, a_radial)
% rlim([-350, 350]);
title('Radiale Beschleunigung')

subplot(3,1,2)
plot(t*w_0, a_circular)
% rlim([-50, 25]);
title('Zirkulare Beschleunigung')

subplot(3,1,3)
plot(t*w_0, a_absolute)
% ylim([-300, 400]);
title('Absolute Beschleunigung')

%% Plot Animation
% Convert polar to cartesian
theta = t*w_0;
[x_r, y_r] = pol2cart(theta, r);

[x_v_radial, y_v_radial] = pol2cart(theta, v_radial);
[x_v_circular, y_v_circular] = pol2cart(theta +(pi/2), v_circular);     % ..and rotate 90°

[x_a_radial, y_a_radial] = pol2cart(theta, a_radial);
[x_a_circular, y_a_circular] = pol2cart(theta + (pi/2), a_circular);    % ..and rotate 90°

% Scales to fit it in plot
scale_v = .03;
scale_a = .0005;
pause on;

figure('units','normalized','outerposition',[0 0 1 1])
for t_ = 1:(period/sampleRate)-1
    
    % Set radial velocity vector
    v_rad_vector_x = [x_r(t_), x_r(t_) + x_v_radial(t_)*scale_v];
    v_rad_vector_y = [y_r(t_), y_r(t_) + y_v_radial(t_)*scale_v];
        
    % Set circular velocity vector
    v_cir_vector_x = [x_r(t_), x_r(t_) + x_v_circular(t_)*scale_v];
    v_cir_vector_y = [y_r(t_), y_r(t_) + y_v_circular(t_)*scale_v];
    
    % Set total velocity vector
    v_vector_x = [x_r(t_), x_r(t_) + (x_v_circular(t_) + x_v_radial(t_))*scale_v];
    v_vector_y = [y_r(t_), y_r(t_) + (y_v_circular(t_) + y_v_radial(t_))*scale_v];
    
    % Set radial acceleration vector
    a_rad_vector_x = [x_r(t_), x_r(t_) + x_a_radial(t_)*scale_a];
    a_rad_vector_y = [y_r(t_), y_r(t_) + y_a_radial(t_)*scale_a];
    
    % Set circular acceleration vector
    a_cir_vector_x = [x_r(t_), x_r(t_) + x_a_circular(t_)*scale_a];
    a_cir_vector_y = [y_r(t_), y_r(t_) + y_a_circular(t_)*scale_a];
    
    % Set total acceleration vector
    a_vector_x = [x_r(t_), x_r(t_) + (x_a_circular(t_) + x_a_radial(t_))*scale_a];
    a_vector_y = [y_r(t_), y_r(t_) + (y_a_circular(t_) + y_a_radial(t_))*scale_a];
    
    % Plot velocity
    subplot(2,2,1);
    plot(x_r, y_r, v_rad_vector_x, v_rad_vector_y, v_cir_vector_x, v_cir_vector_y);
    title_ = ['Radialer und zirkularer Geschwindigkeitsvektor in mm/s (Faktor: ', num2str(scale_v), ')'];
    title(title_);
    pbaspect([1 1 1]);
    xlim([-.04,.04]);
    ylim([-.04,.04]);
    xlabel('mm');
    ylabel('mm');
    
    subplot(2,2,2);
    plot(x_r, y_r, v_vector_x, v_vector_y);
    xlabel('mm');
    ylabel('mm');
    title_ = ['Geschwindigkeitsvektor in mm/s (Faktor: ', num2str(scale_v) ,')'];
    title(title_);
    pbaspect([1 1 1])
    xlim([-.04,.04]);
    ylim([-.04,.04]);
    xlabel('mm');
    ylabel('mm');
    
    % Plot acceleration
    subplot(2,2,3);
    plot(x_r, y_r, a_rad_vector_x, a_rad_vector_y, a_cir_vector_x, a_cir_vector_y);
    title_ = ['Radialer und zirkularer Beschleunigungsvektor in mm/s^2 (Faktor: ', num2str(scale_a), ')'];
    title(title_);
    pbaspect([1 1 1])
    xlim([-.04,.04]);
    ylim([-.04,.04]);
    xlabel('mm');
    ylabel('mm');
    
    subplot(2,2,4);
    plot(x_r, y_r, a_vector_x, a_vector_y);
    title_ = ['Beschleunigungsvektor in mm/s^2 (Faktor: ', num2str(scale_a), ')'];
    title(title_);
    pbaspect([1 1 1])
    xlim([-.04,.04]);
    ylim([-.04,.04]);
    xlabel('mm');
    ylabel('mm');
    
    pause(sampleRate);
end

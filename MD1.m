%% Description

% This is a script to solve the exercise 1 of Machine Dynamics

%% Output

% Formated figures and required values

%% Version

%           author: Ricardo Fitas (rfitas99@gmail.com, University of Porto)
%        copyright: 2021 Ricardo Fitas University of Porto
%    creation date: 30/10/2021
%   Matlab version: R2020b
%          Version: 1.1.1

%% Revision

% V1.0   | 30/10/2021 | Ricardo Fitas | creation
% V1.1   | 01/11/2021 | Ricardo Fitas | exercises resolution
% V1.1.1 | 02/11/2021 | Ricardo Fitas | revision

%% Program
clear       % Delete Workspace
clc         % Clear Command Window
close all   % Close all figures

%% 1) Exercise 1 - Machine Dynamics

%% 1.1.) Constants definition
a = 0.100;         % position of the damper [m]
b = 0.500;         % position of the spring [m]
d = 0.550;         % position of the motor [m]
l = 0.715;         % position of the bar [m]
mb = 1.68;         % mass of the bar [kg]
mm = 0.772;        % mass of the motor [kg]
me = 0.100;        % mass of the eccentric mass [kg]
e = 0.010;         % distance of the eccentric mass to the centre of the motor [m]

% Unknown variables are calculated explicitally in further exercises
syms c             % damping coefficient [Ns/m]
syms k             % stiffness coefficient [N/m]
syms omega         % angular frequency [rad/s]

% system variables
t = 0:0.01:4;                         % time [s]
syms x;                               % displacement [rad]
x_dot = diff(x,1);                    % velocity [m/s]
x_2dot = diff(x,2);                   % acceleration [m/s^2]   

%% 1.a) Establish the (linearized) differential equation

m_eq = (mb*l^2)/3+mm*d^2;               % equivalent mass [kg*m^2]
c_eq = c*a^2;                           % equivalent damping [N*s*m/rad]
k_eq = k*b^2;                           % equivalent stiffness [N*m/rad]
Amplitude = d*me*(omega^2)*e;           % force amplitude [N*m]
f_eq = Amplitude*cos(omega*t);          % equivalent force [N*m]

%Governing equation
m_eq * x_2dot + c_eq * x_dot + k_eq * x == f_eq;

%% 1.b) Establish the damping ratio (xi)
vect_amp = [11.4; 10.3; 9.0; 7.8; 6.8; ...
    5.8;4.9; 3.9; 3.2; 2.4; 1.7];                    % vector storing all first 11 peak values
vect_decr_log = log(vect_amp(1)./vect_amp);          % logarithmic decrement
[m_and_b,more_data] = polyfit(0:10,vect_decr_log,1); % linear regression
delta = m_and_b(1);                                  % the slope represents a more 
                                                     % representative value for the logarithmic decrement 
xi = delta / sqrt(4*pi^2+delta^2);                   % damping ratio

%% 1.c)

T = 0.1909;                         % period [s]
omega_n = 2*pi/(sqrt(1-xi^2)*T);    % undamped natural frequency [rad/s]

%% 1.d)

c_eq = xi*2*m_eq*omega_n;
c = c_eq/(a^2);

%% 1.e)

k_eq = (omega_n^2)*m_eq;
k = k_eq/(b^2);

%% 1.f)

x_0 = 0.050;                    % initial extremety displacement [m]
theta_0 = x_0/l;                % initial angular displacement [rad]
theta_dot_0 = 0;                % initial velocity [rad/s]

% For the equation 
% x(t) = A*exp(-xi*omega_n*t)*cos(omega_d*t-phi)

A = sqrt(((theta_dot_0+xi*omega_n*theta_0)/(omega_n*sqrt(1-xi^2)))^2+theta_0^2);
omega_d = omega_n*sqrt(1-xi^2);
phi = atan2(theta_dot_0+xi*omega_n*theta_0,omega_n*sqrt(1-xi^2)*theta_0);

x = real(A*exp(-xi*omega_n*t).*cos(omega_d*t-phi));
%% 1.f.1.) Plot

% simplified version:
plot_simple(t,x);
% animated version:
% plot_special(t,x,[0, diff(x)]);

%% 1.g)

a = 0.100:0.100:0.7;                % a is now variable
c_eq = c * a.^2;
xi = c_eq/(2*m_eq*omega_n);


%% 1.g.1) Plot

%% 1.g.1.1.1.) Definitions
clr = [236/255 237/255 237/255];
unts = 'normalized';
lnwidth = 2;
fntsz = 18;

%% 1.g.1.1.2.) Positions, titles and labels
pos_fig = [0.01 0.25 .75 .65];
title_graph = 'Damping ratio vs damper position';
xlabel_graph = 'Damper position /m';
ylabel_graph = 'Damping ratio';

%% 1.g.1.2) Plot
figure('color',clr,'units',unts,'position',pos_fig);
graph_plot = scatter(1,1,'*');
axes_graph(1) = gca;
set(axes_graph(1),'FontSize',fntsz);
xlabel(axes_graph(1),xlabel_graph,'fontsize',fntsz);
ylabel(axes_graph(1),ylabel_graph,'fontsize',fntsz);
xlim([0, 1])
title(title_graph,'fontsize',fntsz);
set(graph_plot(1),'linewidth',lnwidth);
set(graph_plot(1),'Parent',axes_graph(1),'XData',a,'YData',xi);
hold on
graph_plot = plot(1,1);
axes_graph(1) = gca;
set(axes_graph(1),'FontSize',fntsz);
xlabel(axes_graph(1),xlabel_graph,'fontsize',fntsz);
ylabel(axes_graph(1),ylabel_graph,'fontsize',fntsz);
xlim([0, 1])
title(title_graph,'fontsize',fntsz);
set(graph_plot(1),'linewidth',lnwidth);
set(graph_plot(1),'Parent',axes_graph(1),'XData',0:1,'YData',ones(size(0:1)));
hold off
legend("Damping ratios for" + newline + "the specified damper" + newline +...
    "positions","Critically-damped" + newline + "situation");

%% 1.h)


t = 0:0.001:1;
x = NaN(size(t,2),size(a,2));
for k = 1:size(a,2)
    c_eq = c * a(k)^2;
    xi = c_eq/(2*m_eq*omega_n);
    if xi > 1
        A1 = theta_0;
        A2 = (xi*omega_n*theta_0+theta_dot_0)./(omega_n*sqrt(xi.^2-1));
        x(:,k) = real(exp(-xi*omega_n*t').*(A1*cosh(omega_n*sqrt(xi^2-1).*t')+A2*sinh(omega_n*sqrt(xi^2-1)*t')));
    else
        omega_d = omega_n*sqrt(1-xi^2);
        A = sqrt(((theta_0*omega_n*xi + theta_dot_0)./(omega_n*sqrt(1-xi^2))).^2+theta_0^2);
        phi = atan2(real(theta_dot_0+xi*omega_n*theta_0),real(omega_n*sqrt(1-xi^2)*theta_0));
        x(:,k) = real(A*exp(-xi*omega_n*t').*cos(omega_d*t'-phi));
    end
end
plot_simple_holded(t,x)
for k = 1:size(a,2)
    c_eq = c * a(k)^2;
    xi = round(c_eq/(2*m_eq*omega_n),2);
    legenda(k) = 'a = ' + string(a(k)) + '; \xi = ' + string(xi);
end
legend(legenda)

%% 1.i) Stationary response of the harmonic-forced regime

%% 1.i.1) Definitions

syms a       % damper position [m]
syms omega   % angular frequency [rad/s]

%% 1.i.2) Other definitions

c_eq = c * a^2;                 % equivalent damping [N*s*m/rad]
xi = c_eq/(2*m_eq*omega_n);     % damping ratio
beta = omega/omega_n;           % frequency ratio

%% 1.i.3) Expression components

dy_amp_fct = ((1-beta^2)^2+(2*xi*beta)^2)^(-.5);     % dynamic amplification factor 
phi_stat = atan2(2*xi*beta,1-beta^2);                % phase angle
Amplitude = d*me*(omega^2)*e;                        % static force amplitude [N*m]
Res_amplitude = Amplitude/k_eq;                      % static response amplitude [rad]
Res_amp = Res_amplitude*dy_amp_fct;                  % Response amplitude [rad]

%% 1.i.4) Response expression
x_stat = Res_amp*cos(omega*t-phi_stat);

%% 1.j)

%% 1.j.1.) For a general damper position a
syms a
%% 1.j.2.) For the damper position of a = 0.100 m
a = 0.100;
c_eq = c * a^2;                           % equivalent damping [N*s*m/rad]
xi = c_eq/(2*m_eq*omega_n);               % damping ratio

beta_omega_r = 1/sqrt(1-2*xi^2);          % frequency ratio for a maximum response amplitude
% Alternative
% beta_omega_r = ...
%    peak_ex(omega_n,xi);                 % frequency ratio for a maximum spring force amplitude
omega_r = omega_n * beta_omega_r;         % Resonance frequency
dy_amp_fct_max = ...
    ((1-beta_omega_r^2)^2 ...
+(2*xi*beta_omega_r)^2)^(-.5);            % dynamic amplification factor 
Amplitude = d*me*(omega_r^2)*e;           % static force amplitude [N*m]
Res_amplitude = Amplitude/k_eq;           % static response amplitude [rad]
Amp_max = Res_amplitude * dy_amp_fct_max; % Maximum response amplitude
omega_amp_max = beta_omega_r*omega_n;     % Frequency for the maximum amplitude

%% 1.k)

%% 1.k.1) Expression of f_k(t)
syms omega
k = k_eq/(b^2);
f_k = k*b*Res_amp*cos(omega*t - phi_stat);     % spring force

%% 1.k.2) Expression of the maximum spring amplitude
F_k_max = k*b*Amp_max;                         % maximum spring amplitude

%% 1.k.3) Plot amplitudes as a function of the frequency

%% 1.k.3.1) Variables

k = k_eq/(b^2);

a = 0.100;                      % damper position
beta = 0:0.01:4;                % frequency ratio

c_eq = c * a^2;                 % equivalent damping [N*s*m/rad]
xi = c_eq/(2*m_eq*omega_n);     % damping ratio
omega = omega_n*beta;           % frequency [rad/s]

dy_amp_fct = ((1-beta.^2).^2+(2*xi*beta).^2).^(-0.5); % dynamic amplification factor
Amplitude = d*me*(omega.^2)*e;                        % static force amplitude [N*m]
Res_amplitude = Amplitude/k_eq;                       % Response amplitude [rad]
Amp_k = k*b*Res_amplitude.*dy_amp_fct;                % spring force amplitude

%% 1.k.3.2) Plot

%% 1.k.3.2.1.) Definitions
clr = [236/255 237/255 237/255];
unts = 'normalized';
lnwidth = 2;
fntsz = 18;

%% 1.k.3.2.2.) Positions, titles and labels
pos_fig = [0.01 0.25 .75 .65];
title_graph = 'Spring amplitude vs Frequency ratio';
xlabel_graph = 'Frequency ratio (\beta)';
ylabel_graph = 'Spring amplitude F_k(\beta)';

%% 1.k.3.2.3.) Plot
figure('color',clr,'units',unts,'position',pos_fig);
graph_plot = plot(1,1,'color','k');
axes_graph(1) = gca;
set(axes_graph(1),'FontSize',fntsz);
xlabel(axes_graph(1),xlabel_graph,'fontsize',fntsz);
ylabel(axes_graph(1),ylabel_graph,'fontsize',fntsz);
title(title_graph,'fontsize',fntsz);
set(graph_plot(1),'linewidth',lnwidth);
set(graph_plot(1),'Parent',axes_graph(1),'XData',beta,'YData',Amp_k);

%% 1.l)

%% 1.l.1) Expression of f_c(t)
syms omega
f_c = omega*c*a*Res_amp*cos(omega*t - phi_stat + pi/2);     % damper force

%% 1.l.2) Expression of the maximum damper amplitude
%% 1.l.2.1.) Calculation of peak value

a = 0.100;

c_eq = c * a^2;                           % equivalent damping [N*s*m/rad]
xi = c_eq/(2*m_eq*omega_n);               % damping ratio

[beta_omega_rc,expressao] = ...
    peak_ex(omega_n,xi,'damper');         % frequency ratio for a maximum damper force amplitude
omega_rc = omega_n * beta_omega_rc;       % Resonance frequency
dy_amp_fct_max_c = ...
    ((1-beta_omega_rc^2)^2+...
    (2*xi*beta_omega_rc)^2)^(-.5);        % Maximum dynamic amplification factor
Amplitude = d*me*(omega_rc^2)*e;          % static force amplitude [N*m]
Res_amplitude = Amplitude/k_eq;           % static response amplitude [rad]
Amp_max_c = Res_amplitude *...
    dy_amp_fct_max_c;                     % Maximum response amplitude
omega_amp_max_c = beta_omega_rc*omega_n;  % Frequency for the maximum amplitude

%% 1.l.2.2.) Maximum damper amplitude
F_c_max = omega_amp_max_c*c*a*Amp_max_c;

%% 1.l.3) Plot amplitudes as a function of the frequency

%% 1.l.3.1) Variables

a = 0.100;                      % damper position
beta = 0:0.01:4;                % frequency ratio

c_eq = c * a^2;                 % equivalent damping [N*s*m/rad]
xi = c_eq/(2*m_eq*omega_n);     % damping ratio
omega = omega_n*beta;           % frequency [rad/s]

dy_amp_fct = ((1-beta.^2).^2+(2*xi*beta).^2).^(-0.5); % dynamic amplification factor
Amplitude = d*me*(omega.^2)*e;                        % static force amplitude [N*m]
Res_amplitude = Amplitude/k_eq;                       % Response amplitude [rad]
Amp_c = omega.*c*a.*Res_amplitude.*dy_amp_fct;        % damper force amplitude

%% 1.l.3.2) Plot

%% 1.l.3.2.1.) Definitions
clr = [236/255 237/255 237/255];
unts = 'normalized';
lnwidth = 2;
fntsz = 18;

%% 1.l.3.2.2.) Positions, titles and labels
pos_fig = [0.01 0.25 .75 .65];
title_graph = 'Damper amplitude vs Frequency ratio';
xlabel_graph = 'Frequency ratio (\beta)';
ylabel_graph = 'Damper amplitude F_c(\omega)';

%% 1.l.3.2.3.) Plot
figure('color',clr,'units',unts,'position',pos_fig);
graph_plot = plot(1,1,'color','k');
axes_graph(1) = gca;
set(axes_graph(1),'FontSize',fntsz);
xlabel(axes_graph(1),xlabel_graph,'fontsize',fntsz);
ylabel(axes_graph(1),ylabel_graph,'fontsize',fntsz);
title(title_graph,'fontsize',fntsz);
set(graph_plot(1),'linewidth',lnwidth);
set(graph_plot(1),'Parent',axes_graph(1),'XData',beta,'YData',Amp_c);

%% 1.m)

syms omega
phi_stat_deg = 135;                              % phase angle [º]
phi_stat_rad = phi_stat_deg*pi/180;              % phase angle [rad]

tan_phi = tan(phi_stat_rad);
beta = omega/omega_n;

assume(omega>=0)
omega = solve(tan_phi == 2*xi*beta/(1-beta^2),omega);     % solve phase equation for beta
omega = double(omega);                                    % final frequency [rad/s]

%% 1.n)

a = 0.100:0.100:0.700;                      % damper position
beta = 0:0.001:4;                           % frequency ratio

omega = omega_n*beta;                       % frequency [rad/s]
amp_ger = NaN(size(beta,2),size(a,2));

for k = 1:size(a,2)
    c_eq = c * a(k)^2;
    xi = c_eq/(2*m_eq*omega_n);
    dy_amp_fct = ((1-beta.^2).^2+(2*xi*beta).^2).^(-0.5); % dynamic amplification factor
    Amplitude = d*me*(omega.^2)*e;                        % static force amplitude [N*m]
    Res_amplitude = Amplitude/k_eq;                       % static response amplitude [rad]
    Res_amp = Res_amplitude.*dy_amp_fct;                  % Response amplitude [rad]
    amp_ger(:,k) = (Res_amp')*l;                          % Displacement at C [m]
end
plot_simple_holded(beta,amp_ger)
for k = 1:size(a,2)
    legenda(k) = 'a = ' + string(a(k));
end
ylim([0,max(max(amp_ger))*1.05]);
title('Response amplitude at C vs Frequency ratio');
xlabel('Frequency ratio (\beta)');
ylabel('Response amplitude at C, X(\omega)/m');
legend(legenda)

%% 1.n.1) Limit calculation

syms beta_s xi_s

x_s = limit((beta_s^2)/...
    sqrt((1-beta_s^2)^2+(2*xi_s*beta_s)^2),beta_s,inf);     % = 1
X_inf = double(d*me*e*x_s*omega_n^2*l/k_eq);                % = .75652 mm             

%% 1.o)

a = 0.100:0.100:0.700;                      % damper position
beta = 0:0.001:4;                           % frequency ratio

omega = omega_n*beta;                       % frequency [rad/s]
amp_ger = NaN(size(beta,2),size(a,2));
k2 = k_eq/(b^2);                            % another variable for the stiffness constant k
for k = 1:size(a,2)
    c_eq = c * a(k)^2;
    xi = c_eq/(2*m_eq*omega_n);
    dy_amp_fct = ((1-beta.^2).^2+(2*xi*beta).^2).^(-0.5); % dynamic amplification factor
    Amplitude = d*me*(omega.^2)*e;                        % static force amplitude [N*m]
    Res_amplitude = Amplitude/k_eq;                       % static response amplitude [rad]
    Res_amp = Res_amplitude.*dy_amp_fct;                  % Response amplitude [rad]
    F_k = Res_amp * k2*b;                                 % Spring force amplitude [N*m/rad]
    amp_ger(:,k) = F_k;                                   
end
plot_simple_holded(beta,amp_ger)
for k = 1:size(a,2)
    legenda(k) = 'a = ' + string(a(k));
end
ylim([0,max(max(amp_ger))*1.05]);
title('Spring force amplitude vs Frequency ratio');
xlabel('Frequency ratio (\beta)');
ylabel('Spring force amplitude, F_k(\omega)/Nm*rad^{-1}');
legend(legenda)

%% 1.o) Limit calculation

F_k_inf = double(k2*b*d*me*e*x_s*omega_n^2/k_eq);             % = 1.1926 Nm

%% 1.p)

a = 0.100:0.100:0.700;                      % damper position
beta = 0:0.001:4;                           % frequency ratio
omega = omega_n*beta;                       % frequency [rad/s]
amp_ger = NaN(size(beta,2),size(a,2));

for k = 1:size(a,2)
    c_eq = c * a(k)^2;
    xi = c_eq/(2*m_eq*omega_n);
    dy_amp_fct = ((1-beta.^2).^2+(2*xi*beta).^2).^(-0.5); % dynamic amplification factor
    Amplitude = d*me*(omega.^2)*e;                        % static force amplitude [N*m]
    Res_amplitude = Amplitude/k_eq;                       % static response amplitude [rad]
    Res_amp = Res_amplitude.*dy_amp_fct;                  % Response amplitude [rad]
    F_c = Res_amp .* omega * c*a(k);                      % Damper force amplitude [N*m*s/rad]
    amp_ger(:,k) = F_c;                                   
end
plot_simple_holded(beta,amp_ger)
for k = 1:size(a,2)
    legenda(k) = 'a = ' + string(a(k));
end
ylim([0,max(max(amp_ger))*1.05]);
title('Damper force amplitude vs Frequency ratio');
xlabel('Frequency ratio (\beta)');
ylabel('Damper force amplitude, F_c(\omega)/Nms*rad^{-1}');
legend(legenda)
pause
%% 2) Save figures

saveas(figure(1),'Fig1.png');
saveas(figure(1),'Fig1.fig');
saveas(figure(2),'Fig2.png');
saveas(figure(2),'Fig2.fig');
saveas(figure(3),'Fig3.png');
saveas(figure(3),'Fig3.fig');
saveas(figure(4),'Fig4.png');
saveas(figure(4),'Fig4.fig');
saveas(figure(5),'Fig5.png');
saveas(figure(5),'Fig5.fig');
saveas(figure(6),'Fig6.png');
saveas(figure(6),'Fig6.fig');
saveas(figure(7),'Fig7.png');
saveas(figure(7),'Fig7.fig');
saveas(figure(8),'Fig8.png');
saveas(figure(8),'Fig8.fig');

%% Functions

function plot_special(time,x_t,v_t)

%% 3.) Plot
%% 3.1.) Initialise figures
run('Exercise_2_3_initialize_figures.m');

%% 3.2.) Draw ground
hold on
run('Exercise_2_4_draw_ground.m');

%% 3.3.) Draw mass
run('Exercise_2_4_draw_mass.m');

%% 3.4.) Draw spring
run('Exercise_2_5_draw_spring.m');

%% 3.5.) Draw damper
run('Exercise_2_6_draw_damper.m');

%% 3.6.) Draw coordinate system
run('Exercise_2_7_draw_cos.m');

%% 3.7.) Plot Animation
% Initialise vectors

x_t_length =length(x_t);
t_plot = NaN(1,x_t_length);
x_t_plot = NaN(1,x_t_length);
v_t_plot = NaN(1,x_t_length);

u = 1;
for k = 1:x_t_length
    cla
    
    % Plot Graph
    
    t_plot(k) = time(u);
    x_t_plot(k) = x_t(u);
    v_t_plot(k) = v_t(u);
    
    set(graph_plot(1),'Parent',axes_graph(1),'XData',t_plot,'YData',x_t_plot);
    set(graph_plot(2),'Parent',axes_graph(2),'XData',t_plot,'YData',v_t_plot);
    
    % Plot ground
    plotcube(axes_ani,dimension_g,position_g,clr_g)
    
    % Plot
    position_m = [x_t(u) 0 0];
    plotcube(axes_ani,dimension_m,position_m,clr_m);
    
    % Plot spring
    spring_head = x_t(u) + dimension_m(1)/2;
    x_pos_spring = phi_s/phi_max*(spring_head - spring_foot) + spring_foot;
    plot3(axes_ani,x_pos_spring,y_pos_spring,z_pos_spring,'b','linewidth',lnwidth);
    
    % Plot Damper
    damper_head = x_t(u) + dimension_m(1)/2;
    plotdamper(stroke_length_max,damper_foot,damper_head,y_offset_d,clr_d,lnwidth);
    
    % Plot cos
    plotcos(x_ar,variable_cos,clr_cos,lnwidth,fntsz);
    
    % Rotate Views
    view(90,-90);
    
    title(title_ani,'fontsize',fntsz);
    xlabel(xlabel_ani,'fontsize',fntsz);
    
    drawnow
    
    u = u + 1;
end

end

function plot_simple(t,x,v)

%% 0.) Error handle
if nargin ~= 2 && nargin ~= 3
    error('Incorrect number of input arguments')
end
%% 1.) Definitions
clr = [236/255 237/255 237/255];
unts = 'normalized';
lnwidth = 2;
fntsz = 18;

%% 1.1.) Positions, titles and labels
pos_fig = [0.01 0.25 .75 .65];
if nargin == 3
    title_graph = 'Displacement and velocity vs time';
else
    title_graph = 'Displacement vs time';
end

xlabel_graph = 'Time t/s';
ylabel_graph{1} = 'Displacement x/rad';
ylabel_graph{2} = 'Velocity v/rad*s^{-1}';

%% 2.) Plot
%% 2.1.) Initialise Figures
fig = figure('color',clr,'units',unts,'position',pos_fig);
graph_plot = plot(1,1,1,1);
set(graph_plot(1),'color','k','linewidth',lnwidth);
if nargin == 3
    set(graph_plot(2),'color','r','linewidth',lnwidth);
end
axes_graph(1) = gca;
set(axes_graph(1),'FontSize',fntsz);
if nargin == 3
    axes_graph(2) = axes('Position',axes_graph(1).Position,'YAxisLocation','right','YColor','r','Color','none','XTickLabel',[],'fontsize',fntsz);
end
xlabel(axes_graph(1),xlabel_graph,'fontsize',fntsz);
ylabel(axes_graph(1),ylabel_graph(1),'fontsize',fntsz);
if nargin == 3
    ylabel(axes_graph(2),ylabel_graph(2),'fontsize',fntsz);
end
title(title_graph,'fontsize',fntsz);

x_t_max_limit = max(abs(x)) + 0.05*max(abs(x));
ylim(axes_graph(1),[-x_t_max_limit,x_t_max_limit]);

if nargin == 3
    v_t_max_limit = max(abs(v)) + 0.05*max(abs(v));
    ylim(axes_graph(2),[-v_t_max_limit,v_t_max_limit]);
end

xlim(axes_graph(1),[t(1) t(end)]);

set(graph_plot(1),'Parent',axes_graph(1),'XData',t,'YData',x);

if nargin == 3
    set(graph_plot(2),'Parent',axes_graph(2),'XData',t,'YData',v);
end

end

function plot_simple_holded(t,x)

%% 0.) Error handle
if nargin ~= 2
    error('Incorrect number of input arguments')
end
%% 1.) Definitions
clr = [236/255 237/255 237/255];
unts = 'normalized';
lnwidth = 2;
fntsz = 18;

%% 1.1.) Positions, titles and labels
pos_fig = [0.01 0.25 .75 .65];
title_graph = 'Displacement vs time';
xlabel_graph = 'Time t/s';
ylabel_graph{1} = 'Displacement x/rad';

%% 2.) Plot
%% 2.1.) Initialise Figures
fig = figure('color',clr,'units',unts,'position',pos_fig);
for k = 1:size(x,2)
    graph_plot = plot(1,1);
    hold on
    axes_graph(1) = gca;
    set(axes_graph(1),'FontSize',fntsz);
    xlabel(axes_graph(1),xlabel_graph,'fontsize',fntsz);
    ylabel(axes_graph(1),ylabel_graph(1),'fontsize',fntsz);
    title(title_graph,'fontsize',fntsz);
    x_t_max_limit = max(max(abs(x)) + 0.05*max(max(abs(x))));
    ylim(axes_graph(1),[-x_t_max_limit,x_t_max_limit]);
    xlim(axes_graph(1),[t(1) t(end)]);
    set(graph_plot(1),'linewidth',lnwidth);
    set(graph_plot(1),'Parent',axes_graph(1),'XData',t,'YData',x(:,k));
end
end


function [y,y2] = peak_ex(on,xi2,tex) 

%% 0.) Input evaluation
if nargin ~= 3
    error('Incorrect number of input arguments')
end
if ~(strcmp(tex, 'spring') || strcmp(tex,'damper'))
    warning('Response / spring max amplitude will be assumed as output')
end

%% 1.) Definitions
syms d me omega omega_n e xi k b c a

%% 2.) Components
beta = omega/omega_n;                                 % frequency ratio
dy_amp_fct = ((1-beta^2)^2+(2*xi*beta)^2)^(-0.5);     % dynamic amplification factor
Force_amp = d*me*(omega^2)*e;                         % static force amplitude [N*m]
Res_amplitude = Force_amp/(k*b^2);                    % static response amplitude [rad]
Res_amp = Res_amplitude*dy_amp_fct;                   % Response amplitude [rad]

%% 3.) Expression selection
if strcmp(tex, 'spring')
    %% 3.1.) For the spring force

    F_k_amp = k*b*Res_amp;                            % Spring force amplitude
    assume(omega>0)                                   % frequencies are always positive
    assumeAlso(omega_n>0)
    omega_der = solve(0==diff(F_k_amp,omega), omega);
    y2  = omega_der/omega_n;                          % expression for the frequency ratio
end

if strcmp(tex, 'damper')
    %% 3.2.) For the damper force

    F_c_amp = c*a*omega*Res_amp;                      % Damper force amplitude
    assume(omega>0)                                   % frequencies are always positive and real
    assumeAlso(omega_n>0)
    assumeAlso(omega,'real')
    assumeAlso(omega_n,'real')
    assumeAlso(xi>0)                                  % damping ratio is always positive and real
    assumeAlso(xi,'real')
    omega_der = solve(0==diff(F_c_amp,omega), omega);
    y2  = omega_der/omega_n;                          % expression for the frequency ratio
end
omega_n = on;
xi = xi2;
% If two solutions are obtained, the true solution is the minimum of both
% frequency ratio values
y = min(subs(y2));      

end
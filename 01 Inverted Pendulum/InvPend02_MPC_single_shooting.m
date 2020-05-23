%% Single shooting MPC for the classical swing-up inverted pendulum problem
% fork from https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi
% all MPC with CASADI code, credits to Mohamed

% inv. pendulum by Helon Ayala, 05/2020

% animation credits to Renato Coral, UnB
% link to thesis: https://repositorio.unb.br/bitstream/10482/34629/1/2018_RenatoCoralSampaio.pdf

clear
close all
clc

import casadi.* % import casadi libs (make sure its on path)

% define the MPC - optimal control problem

% pendulum parameters
n_states = 4;
n_controls = 1;
M = 10;
m = 80;
c = 0.1;
J = 100;
l = 1;
gamma = 0.01;
g = 9.8;

% simulation parameters
x0 = [0 ; pi; 0; 0];    % initial condition  = pendulum facing down
xs = zeros(n_states,1); % setpoint = desired final position (pendulum up)
T = 0.01; % sampling time [s]
sim_tim = 20; % Maximum simulation time
u_max = 800; u_min = -u_max;        % force constraint (N)

% MPC parameterization
N = 200; % prediction horizon (2 seconds)
Q = zeros(n_states); % weighing matrices (states)
Q(2,2) = 100; % greater penalty for theta
Q([1 3 4],[1 3 4]) = 1; 
R = zeros(n_controls); % weighing matrices (controls)
Rf = 0; % weighing matrices delta U

% modeling with Casadi:
p      = SX.sym('p'); 
th     = SX.sym('th');
p_dot  = SX.sym('p_dot'); 
th_dot = SX.sym('th_dot');
u = SX.sym('u');
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + n_states); % parameters (which include the initial and the reference state of the robot)
X = SX.sym('X',n_states,(N+1)); % A Matrix that represents the states over the optimization problem.

states = [p;
          th;
          p_dot;
          th_dot];

controls = u; 

% system r.h.s (state equations)
rhs = [p_dot;
       th_dot;
       -(c*p_dot - u + l*m*th_dot^2*sin(th) + (l*m*cos(th)*(gamma*th_dot - g*l*m*sin(th)))/(m*l^2 + J))/(M + m - (l^2*m^2*cos(th)^2)/(m*l^2 + J));
       -(gamma*th_dot - g*l*m*sin(th) + (l*m*cos(th)*(l*m*sin(th)*th_dot^2 - u + c*p_dot))/(M + m))/(J + l^2*m - (l^2*m^2*cos(th)^2)/(M + m))];

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)

% compute solution symbolically
X(:,1) = P(1:n_states); % initial state
for k = 1:N
    st = X(:,k);
    con = U(:,k);
    f_value  = f(st,con);
    st_next  = st + (T*f_value); % Euler derivative approximation
    X(:,k+1) = st_next;
end
% this function to get the optimal trajectory knowing the optimal solution
ff=Function('ff',{U,P},{X});

obj = 0; % Objective function
g = [];  % constraints vector

% compute objective
for k=1:N
    st = X(:,k);  
    con = U(:,k);
    obj = obj+(st-P(n_states+1:end))'*Q*(st-P(n_states+1:end)) + con'*R*con; % calculate obj
end
obj = obj + diff(U)*Rf*diff(U)'; % penalize delta u

% compute constraints % no constraints for the states
% for k = 1:N+1   % box constraints due to the map margins
%     g = [g ; X(1,k)];   %state x
%     g = [g ; X(3,k)];   %state y
% end

% make the decision variables one column vector
OPT_variables = reshape(U,N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter    = 100;
opts.ipopt.print_level = 0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-4;
opts.ipopt.acceptable_obj_change_tol = 1e-4;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;
% % inequality constraints (state constraints)
% args.ubg = [0.5 pi];   % upper bound of the states x and y 
% args.lbg = -args.ubg;  % lower bound of the states x and y
args.lbg = [];  % lower bound of the states x and y
args.ubg = [];   % upper bound of the states x and y 

% input constraints
args.lbx(1:N,1) = u_min;
args.ubx(1:N,1) = u_max;

% run MPC
t0 = 0;

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,1);  % two control inputs 

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-2 and the number of mpc steps is less than its maximum
% value.
main_loop = tic;
while(norm((x0-xs),2) > 1e-2 && mpciter < sim_tim / T)
% while (mpciter < sim_tim / T)
    args.p   = [x0;xs]; % set the values of the parameters vector
    args.x0 = reshape(u0',N,1); % initial value of the optimization variables
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
                 'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    %toc
    u = reshape(full(sol.x)',1,N)';
    ff_value = ff(u',args.p); % compute OPTIMAL solution TRAJECTORY
    xx1(:,1:n_states,mpciter+1)= full(ff_value)';
    
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    [t0, x0, u0] = f_shift(T, t0, x0, u,f); % get the initialization of the next optimization step
    
    xx(:,mpciter+2) = x0; 
    J(:,mpciter+1) = full(sol.f);   % keeps cost function value
    
    
    fprintf('Iter %d out of %d\n',mpciter,sim_tim / T)
    mpciter = mpciter + 1;
end
main_loop_time   = toc(main_loop);
ss_error         = norm((x0-xs),2);
average_mpc_time = main_loop_time/(mpciter+1);

xx = xx(:,1:end-1);

%% plot
close all
ttl = {'p','th','p dot','th dot'};

fname = ['MPC_sing_shooting_N_' num2str(N) '_T_' num2str(T)];

figure(1)
for i =1:n_states
    subplot(n_states+1,1,i)
    plot(t,xx(i,:),'k','linewidth',2), hold on
    yline(xs(i),'r--','linewidth',1.5),
    grid on
    hold off
    title(ttl{i})
end
subplot(n_states+1,1,i+1)
plot(t,u_cl,'k','linewidth',2), hold on
yline(u_max,'r--','linewidth',1.5),
yline(0,'r--','linewidth',1.5),
yline(u_min,'r--','linewidth',1.5),
grid on
title('u')
% saveas(gcf,[fname '.png'])

% f_animate(t,xx,u_cl,u_min,u_max,J,[fname '.avi']) % generate animation
% 
% implay(fname,60)

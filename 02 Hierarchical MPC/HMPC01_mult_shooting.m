%% Multiple shooting MPC for hierarchical MPC
% fork from https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi
% all MPC with CASADI code, credits to Mohamed

% hierarchical MPC by Helon Ayala, 05/2020

clear
close all
clc

% addpath(genpath('/home/helon/Documents/MATLAB/casadi/'));
import casadi.* % import casadi libs (make sure its on path)

% define the MPC - optimal control problem

% model parameters -- single-link robot (Astrom p. 156) - sist lento ts = 20
J1 = 10/9;
J2 = 10;
d = 0.1;
ki = 1;
w0 = 1;
alpha = J1/(J1+J2);
beta1 = d/(J1*w0);
beta2 = d/(J2*w0);
gamma = ki/(J1*w0);
% delta = 1/(J1*w0);

A = w0 * [0 1 -1;
          alpha-1   -beta1    beta1;
          alpha      beta2   -beta2]; 
B = [0;gamma;0];
C = [0 0 w0];
D = 0;

Mo = [A B;
    -C 0]

det(Mo)

mu = [-5 -20 -21 -22];

% projeto com integrador no ramo MA
Ah = [A zeros(3,1);-C 0]; 
Bh = [B;0];
Kh = place(Ah, Bh,mu);
K  = Kh(1:3);
ki = -Kh(4);

% sist MF
A2 = [A-B*K   B*ki;
      -C 0];
B2 = [0;0;0;1];
C2 = [eye(4);   % x, epsilon
      -K ki];   % u
D2 = zeros(5,1);


% %%
% mu = [-5 -5 -20];
% k = 50; m = 20;
% A = [0 1;-k/m 0];
% B = [0;1/m];
% C = [1 0];
% D = 0;
% 
% % projeto com integrador no ramo MA
% Ah = [A zeros(2,1);-C 0]; 
% Bh = [B;0];
% Kh = acker(Ah, Bh,mu);
% K  = Kh(1:2);
% ki = -Kh(3);
% 
% % sist MF
% A2 = [A-B*K   B*ki;
%       -C 0];
% B2 = [0;0;1];
% C2 = [eye(3);   % x, epsilon
%       -K ki];   % u
% D2 = zeros(4,1);

sysMA = ss(A,B,C,D);
sysMF = ss(A2,B2,C2,D2);

step(sysMA,sysMF,5);
legend('Open loop','CL')

%%
n_states = 5;
n_controls = 1;

% simulation parameters
x0 = zeros(n_states,1);   % initial condition  = pendulum facing down
w2s = 1; % 1 rad/s
xs = w2s;   % setpoint = desired velocity in the robot arm

T = 0.01; % sampling time [s]
sim_tim = 20; % Maximum simulation time
u_max = 10; u_min = -u_max;        % force constraint (N)

% MPC parameterization
N = 200; % prediction horizon (2 seconds)
Q = eye(n_states); % weighing matrices (states)
Q(3,3) = 100; % greater penalty for w3
R = zeros(n_controls); % weighing matrices (controls)
Rf = 0; % weighing matrices delta U

% modeling with Casadi:
p      = SX.sym('p'); 
th     = SX.sym('th');
p_dot  = SX.sym('p_dot'); 
th_dot = SX.sym('th_dot');

x = SX.sym('x',n_states);
u = SX.sym('u',n_controls);
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + n_states); % parameters (which include the initial and the reference state of the robot)
X = SX.sym('X',n_states,(N+1)); % A Matrix that represents the states over the optimization problem.

states = x;
controls = u; 

% system r.h.s (state equations)
rhs = A2*x + B2*u;

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
% 
% % compute solution symbolically
% X(:,1) = P(1:n_states); % initial state
% for k = 1:N
%     st = X(:,k);
%     con = U(:,k);
%     f_value  = f(st,con);
%     st_next  = st + (T*f_value); % Euler derivative approximation
%     X(:,k+1) = st_next;
% end
% % this function to get the optimal trajectory knowing the optimal solution
% ff=Function('ff',{U,P},{X});

obj = 0; % Objective function
g = [];  % constraints vector

% compute objective
st  = X(:,1); % initial state
g = [g;st-P(1:n_states)]; % initial condition constraints
for k=1:N
    st = X(:,k);  
    con = U(:,k);
    obj = obj+(st-P(n_states+1:end))'*Q*(st-P(n_states+1:end)) + con'*R*con; % calculate obj
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end
% obj = obj + diff(U)*Rf*diff(U)'; % penalize delta u

% make the decision variable one column vector
OPT_variables = [reshape(X,n_states*(N+1),1);reshape(U,n_controls*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter    = 2000;
opts.ipopt.print_level = 0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;

args.lbg = zeros(1,n_states*(N+1)); % Equality constraints
args.ubg = zeros(1,n_states*(N+1)); 

args.lbx = zeros(1,n_states*(N+1) + n_controls*N);
args.ubx = zeros(1,n_states*(N+1) + n_controls*N);

args.lbx(1:n_states*(N+1)) = -inf;  % state constraints
args.ubx(1:n_states*(N+1)) = inf;     
args.lbx(n_states*(N+1)+1:end) = u_min; % input constraints
args.ubx(n_states*(N+1)+1:end) = u_max;

% run MPC
t0 = 0;

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

% initialization of the decision variables
u0 = zeros(N,1);  % two control inputs 
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

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
    args.x0  = [reshape(X0',n_states*(N+1),1);reshape(u0',n_controls*N,1)];
    
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
                 'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    

    X0 = reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)';     % get solution TRAJECTORY
    u  = reshape(full(sol.x(n_states*(N+1)+1:end))',n_controls,N)'; % get only controls from the solution
    
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    [t0, x0, u0] = f_shift(T, t0, x0, u,f); % get the initialization of the next optimization step
    
    xx(:,mpciter+2) = x0; 
    J(:,mpciter+1) = full(sol.f);   % keeps cost function value
        
    X0 = [X0(2:end,:);X0(end,:)]; % Shift trajectory to initialize the next step
    
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

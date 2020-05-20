%% get the state equations for the inverted pendulum

clear
close all
clc

syms p_dot p_ddot th_dot th_ddot p th u
syms M m c J l gamma g b

% dynamic eqns.
% see https://www.cds.caltech.edu/~murray/books/AM05/pdf/am08-complete_22Feb09.pdf
% p. 36, eq. (2.9)
eq1 = (M+m)*p_ddot - m*l*cos(th)*th_ddot + c*p_dot + m*l*sin(th)*th_dot^2 - u
eq2 = -m*l*cos(th)*p_ddot + (J+m*l^2)*th_ddot + gamma*th_dot - m*g*l*sin(th) 

% isolate accelerations (higher order derivatives)
p_ddot_is  = solve(eq1,p_ddot)   
th_ddot_is = solve(eq2,th_ddot)  

% inspect
pretty(p_ddot_is)
pretty(th_ddot_is)

% manipulate eqns (get only 1 eq. / acceleration)
eq1_a = subs(eq1,th_ddot,th_ddot_is)    
eq1_final = solve(eq1_a,p_ddot)  
eq2_a = subs(eq2,p_ddot,p_ddot_is)    
eq2_final = solve(eq2_a,th_ddot) 

x = [p;th;p_dot;th_dot];


% finally, state equations
xp =subs([ p_dot;
    th_dot;
    eq1_final;    
    eq2_final]);

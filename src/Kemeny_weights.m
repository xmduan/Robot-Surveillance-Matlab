clc;
clear all;
close all;



n = 4;
% A = [1 1 0 1;1 1 1 0;0 1 1 1;1 0 1 1];
A = [1 1 1 0;1 1 0 1;1 0 1 1;0 1 1 1];
pi_v = [1/4 1/4 1/4 1/4]';

% pi_v = rand(4,1);
% pi_v = pi_v/sum(pi_v);
PI = diag(pi_v);
q = pi_v.^(1/2);

%%
% cvx_begin
% %     cvx_solver sedumi
%     variable Q(n,n) symmetric
%     variable X(n,n)
%     minimize(trace(X))
%     subject to
%     trace(X) >= trace_inv(eye(n,n)-Q+q*q');    
%     PI^(-1/2)*Q*PI^(1/2)*ones(n,1) == ones(n,1);
%     PI^(-1/2)*Q*PI^(1/2).*(ones(4,4)-A) == 0;
%     Q >= zeros(n,n);
% cvx_end
% P1 = PI^(-1/2)*Q*PI^(1/2);
% Q1 = Q;






%%


cvx_begin
%     cvx_solver sedumi
    variable Q(n,n) symmetric
    minimize(trace_inv(eye(n,n)-Q+q*q'))
    subject to
    PI^(-1/2)*Q*PI^(1/2)*ones(n,1) == ones(n,1);
    PI^(-1/2)*Q*PI^(1/2).*(ones(4,4)-A) == 0;
    Q >= zeros(n,n);
cvx_end
P1 = PI^(-1/2)*Q*PI^(1/2);
Q1 = Q;


%% Warning: the variable has to be PI^(1/2)*P*PI^(-1/2).
epsilon = 10^(-10);
cvx_begin sdp
%     cvx_solver SeDuMi
    cvx_precision best
    cvx_solver SDPT3
    variable Q(n,n) symmetric
    variable X(n,n) semidefinite
    minimize(trace(X))
    subject to        
        [eye(n,n)-Q+q*q', eye(n,n);
            eye(n,n), X] >= -epsilon * eye(2*n,2*n);
        PI^(-1/2)*Q*PI^(1/2)*ones(n,1) == ones(n,1);
        PI^(-1/2)*Q*PI^(1/2).*(ones(4,4)-A) == 0;        
        Q >= zeros(n,n);
cvx_end
P2 = PI^(-1/2)*Q*PI^(1/2);
Q2 = Q;
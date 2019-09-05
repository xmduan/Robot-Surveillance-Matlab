function P = RandomMC(n,A,b,C,d,number)


% sampling from the polytope defined by
% Ax <=b; Cx=d
% A: m*n; C:p*n
% b: m*1, d:p*1
% A,b,C,d should be the input


% n = 3; % number of components in x

% C = [1 1 1]; % coeffieicnt matrix of the equality constraints Cx = d
% d = 1;       % right-hand side of equality Cx = d
% A = -eye(n); % coefficient matrix of the inequality constriants. Ax <= b, nonegative constraints
% A = [A;-1 1 0;0 -1 1]; % other constraints
% b = [0;0;0;0;0];        % right-hand side of inequality Ax <= b

P = [];

C = sym(C);

C_pinverse = pinv(C);   % pseudo-inverse of C

F = eye(n)-C_pinverse*C; % F = I - C^* C

r = rank(F);  %Compute the rank of I-C^* C

if r ==0 % which means that A_eq is invertible and there is only one solution
    P = double(C_pinverse * d);
else
    
    
    % [Q R] = qr(F);%QR decomposition of I-C^* C
    % D = Q(:,1:r); % Pick the first r columns of Q matrix of I-C^* C
    
    D = colspace(F);
    
    D = double(D);
    C_pinverse = double(C_pinverse);
    
    
    % x = C^* d + D*y, note that y \in R^r
    % sampling in r-dimensional space
    % implementing HAR in r dimensional space with the following constraints
    % Ax <=b   ->  A(C^* d + D*y) <= b
    % rearraging: AD*y <= b - AC^* d
    % Or A_r y <= b_r
    % A_r has
    % b_r has the same dimension as b
    A_r = A*D;
    b_r = b - A*C_pinverse*d;
    
    %% sampling in A_r y <= b_r uniformly using HAR
    
    
    % generating an interior point:y0, such that A_r y0 <= b_r
    % Finding an initial point using a slack formultion
    % max \delta
    % s.t. A_r x + e = b_r
    %      \delta <= e_i
    % Also can implementing a randomized version by setting e_i to be c_i e_i,
    % where c_i \in (0,1] is random
    [tmp_m,tmp_n] = size(A_r);
    f = zeros(1,tmp_n + tmp_m + 1);
    f(tmp_n + tmp_m + 1) = -1;
    A_ineq = [zeros(tmp_m,tmp_n) -eye(tmp_m,tmp_m) ones(tmp_m,1)];
    b_ineq = zeros(tmp_m,1);
    A_eq = [A_r eye(tmp_m,tmp_m) zeros(tmp_m,1)];
    b_eq = b_r;
    
    
    
    A_eq = double(A_eq);
    b_eq = double(b_eq);
    
    x = linprog(f,A_ineq,b_ineq,A_eq,b_eq);
    y0 = x(1:tmp_n,1);
    
    
    
    
    k=0;
    
    while k <number
        
        
        % for simplicity, our set of sampling direction using the standard unit
        % vector.
        %Actually, the sampling direction can also be generated randomly:
        % generate a vector whose components are standard normal r.v.
        % Then normalize this vectore. The resulting vector is uniformly
        % distributed on the hypersphere
        %
        % d_set = eye(r);
        % d_index = ceil(r*rand);
        % d_sampling = d_set(:,d_index);
        
        % ramdon direction
%         d_sampling = zeros(r,1);
%         for i = 1:r
%             d_sampling(i) = normrnd(1,1);
%         end
%         d_sampling = d_sampling/sqrt(sum(d_sampling.^2));
%         
        d_sampling = normrnd(1,1,r,1);
        d_sampling = d_sampling/sqrt(sum(d_sampling.^2));
        
        
        
        % determining the bounds for l such that y0+ld_sampling is still in the polytope
        % l\in [L0,L1];
        v = b_r - A_r*y0;
        u = A_r * d_sampling;
        count1 = 1;
        count2 = 1;
        ratio_pos = zeros(tmp_m,1);
        ratio_neg = zeros(tmp_m,1);
        for i = 1:tmp_m
            if u(i)>0
                ratio_pos(count1) = v(i)/u(i);
                count1 = count1 +1;
            elseif u(i)<0
                ratio_neg(count2) = v(i)/u(i);
                count2 = count2 + 1;
            end
        end
        % We should have been more careful here since if the polytope is not
        % bounded, L0 may be negative infinity, and L1 may be positive infinity,
        % but we should be fine here because all the polytopes we considered are
        % bounded.
        
        L0 = max(ratio_neg(1:count2-1));
        L1 = min(ratio_pos(1:count1-1));
        l = rand*(L1-L0)+L0;
        
        y = y0+ l * d_sampling;
        
        x = C_pinverse * d + D*y;
        P = [P x];
        
        k = k+1;
        y0=y;
        
    end
    
end
function P_op=Rushabh_op_sdp(A,PI,W)
% if size(phi,2)~=1
%     error('the stationary distribution you have input is not a column vector');
% end
% if single(sum(phi))~=1
%     error('the stationary distribution you have input is illegal');
% end
% n=size(phi,1);
% switch nargin
%     case 1
%         W=ones(n,n);
%     case 0
%         error('please input your stationary distribution of markov chain');        
% end
n=size(A,2);
A_inv=ones(n,n);
A_inv(A>0)=0;

cvx_begin sdp
%     variable Y(n,n) symmetric
    variable Q(n,n) symmetric
    variable X(n,n) semidefinite
    variable t
    minimize(trace(X))
%     tmp=[t*(eye(n)+sqrt(PI)*sqrt(PI)')-sqrt(diag(PI))*Y*diag(PI)^(-1/2) eye(n);eye(n) X];
      tmp=[t*(eye(n)+sqrt(PI)*sqrt(PI)')-Q eye(n);eye(n) X];
    subject to
    tmp+tmp'>=0;
%     for i=1:n
%         Y(:,i)>=0;
%     end
%     for i=1:n
%         Y(:,i)-t*ones(n,1)<=0;
%     end
Q>=zeros(n,n);
 diag(PI)^(-1/2)*Q*diag(PI)^(1/2)*ones(n,1) == t*ones(n,1);
diag(PI)^(-1/2)*Q*diag(PI)^(1/2).*A_inv == 0;
%     Y*ones(n,1)==t*ones(n,1);
    t>=0;
%     Y.*repmat(PI,[1,n])==Y'.*(repmat(PI',[n,1]));
%     X>=0;
%     Y.*A_inv==zeros(n,n);
    PI'*(Q.*W)*ones(n,1)==1;

cvx_end
P_op=diag(PI)^(1/2)*Q*diag(PI)^(-1/2)/t;
end
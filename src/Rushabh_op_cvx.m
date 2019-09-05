function P_op=Rushabh_op_cvx(A,PI,W)
n=size(A,2);
A_inv=ones(n,n);
A_inv(A>0)=0;

cvx_begin 
    variable P(n,n)
    
%     part1=PI'*(P.*W)*ones(n,1);
%     part2=trace_inv(eye(n)+sqrt(PI)*sqrt(PI)'-sqrt(diag(PI))*P*diag(PI)^(-1/2))
    minimize(PI'*P.*W*ones(n,1)*trace_inv(eye(n)+sqrt(PI)*sqrt(PI)'-sqrt(diag(PI))*P*diag(PI)^(-1/2)))
    subject to
    for i=1:n
        P(:,i)>=0;
    end
    for i=1:n
        P(:,i)<=ones(n,1);
    end
    P*ones(n,1)==ones(n,1);
    P.*repmat(PI,[1,n])==P'.*(repmat(PI',[n,1]));
%     for i=1:n
%         X(i,:)>=0;
%     end
    P.*A_inv==zeros(n,n);
cvx_end
P_op=P;
end
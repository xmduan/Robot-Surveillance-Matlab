function P_op=Rushabh_op_norm_cvx(A,PI)

n=size(A,2);
A_inv=ones(n,n);
A_inv(A>0)=0;

cvx_begin 
    variable P(n,n) 
%     variable X(n,n) 
%     variable t
    minimize( trace_inv((eye(n)+sqrt(PI)*sqrt(PI)')-sqrt(diag(PI))*P*diag(PI)^(-1/2)))
%     tmp=[(eye(n)+sqrt(PI)*sqrt(PI)')-sqrt(diag(PI))*P*diag(PI)^(-1/2) eye(n);eye(n) X];
    Ones=ones(n,n);
    Zeros=zeros(n,n);
    subject to
    P(:)<=Ones(:);
    P(:)>=Zeros(:);
    P*ones(n,1)==ones(n,1);
%     t>=0;
    P.*repmat(PI,[1,n])==P'.*(repmat(PI',[n,1]));

    P.*A_inv==zeros(n,n);
cvx_end
P_op=P;
end
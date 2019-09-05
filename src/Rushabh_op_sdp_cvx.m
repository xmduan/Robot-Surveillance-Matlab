function P_op=Rushabh_op_sdp_cvx(A,PI)

n=size(A,2);
A_inv=ones(n,n);
A_inv(A>0)=0;

cvx_begin sdp
    variable P(n,n) 
    variable X(n,n) 
%     variable t
    minimize( trace(X))
    tmp=[(eye(n)+sqrt(PI)*sqrt(PI)')-sqrt(diag(PI))*P*diag(PI)^(-1/2) eye(n);eye(n) X];
    Ones=ones(n,n);
    Zeros=zeros(n,n);
    subject to
    tmp+tmp'>=0;
    P(:)<=Ones(:);
    P(:)>=Zeros(:);
    P*ones(n,1)==ones(n,1);
%     t>=0
    P.*repmat(PI,[1,n])==P'.*(repmat(PI',[n,1]));
    P.*A_inv==zeros(n,n);
cvx_end
P_op=P;
end
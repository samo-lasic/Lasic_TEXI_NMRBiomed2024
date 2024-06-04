function Dw = fDwSphere(w2,R,D0,alpha,Nterms)

Ks = BesselKernelsSphere(Nterms);

% B=2*(R./Ks).^2./(Ks.^2-2);
% a=(Ks/R).^2;

Ks4 = Ks.^4;
C =2./(Ks.^2-2);

wmat2 = repmat(w2(:),1,length(R));
Rmat = repmat(R(:),1,length(w2))';
wR2 = D0.^2./Rmat.^4;

Dw = zeros(size(wmat2));
for n = 1:Nterms
    Dw = Dw + C(n)*wmat2./(Ks4(n)*wR2 + wmat2);
    %Dw(i)=D0*alpha+D0*(1-alpha)*sum(C*w(i)^2./(a.^2*D0^2+w(i)^2));
end

Dw = D0*alpha+D0*(1-alpha)*Dw;


function f=BesselKernelsSphere(n)
for l=1:n
f(l)=fzero(inline('x*besselj(1/2,x)-2*besselj(3/2,x)'),2+(l-1)*pi);
end


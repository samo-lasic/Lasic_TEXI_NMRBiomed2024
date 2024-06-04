function Dw = fDwCylinder(w2,R,D0,alpha,Nterms)

Kc = BesselKernelsCylinder(Nterms);

% B=2*(R./Kc).^2./(Kc.^2-1);
% a=(Kc/R).^2;

Kc4 = Kc.^4;
C =2./(Kc.^2-1);

wmat2 = repmat(w2(:),1,length(R));
Rmat = repmat(R(:),1,length(w2))';
wR2 = D0.^2./Rmat.^4;

Dw = zeros(size(wmat2));
for n = 1:Nterms
    Dw = Dw + C(n)*wmat2./(Kc4(n)*wR2 + wmat2);
    %Dw(i)=D0*alpha+D0*(1-alpha)*sum(a.*B*w(i)^2./(a.^2*D0^2+w(i)^2));
end

Dw = D0*alpha+D0*(1-alpha)*Dw;

function g=BesselKernelsCylinder(n)
for l=1:n
g(l)=fzero(inline('besselj(0,x)-besselj(1,x)/x'),2+(l-1)*pi);
end
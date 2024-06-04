function Dw = fDwPlane(w2,d,D0,alpha,Nterms)
% d = distance between planes
k = 1:Nterms;
% B = 8*d^2./(2*k-1).^4/pi^4;
% a = (2*k-1).^2*pi^2/d^2;

wmat2 = repmat(w2(:),1,length(d));
dmat = repmat(d(:),1,length(w2))';
wd2 = D0.^2./dmat.^4;

K2 = (2*k-1).^2;
C1 = 8./K2/pi^2;
C2 = K2.^2*pi^4;

Dw = zeros(size(wmat2));
for n = 1:Nterms
    %f(i)=D0*alpha+D0*(1-alpha)*sum(a.*B*w(i)^2./(a.^2*D0^2+w(i)^2));
    Dw = Dw + C1(n)*wmat2./(C2(n)*wd2 + wmat2);
end;

Dw = D0*alpha+D0*(1-alpha)*Dw;
function s = power_spectrum(qIN, dt, NFT, ind_PS)
q0pad = zeros(NFT-length(qIN),1);

%         figure(1),clf
%         subplot(2,1,1), plot(gp)
%         subplot(2,1,2), plot(qp)

%         bp(m,n) = sum(qp.^2)*dt;
%         bn(m,n) = sum(qn.^2)*dt;

q = [qIN; q0pad];

s = abs(fftshift(fft(q))).^2;

s = s(ind_PS);
s = s/sum(s);

%         figure(1),clf
%         subplot(2,1,1), plot(fPS,s)
%         subplot(2,1,2), plot(fPS,cumsum(s))

%ADC = Dw*s;
%Vw(cnt) = (2*pi)^2*sum(s.*fPS.^2);
end
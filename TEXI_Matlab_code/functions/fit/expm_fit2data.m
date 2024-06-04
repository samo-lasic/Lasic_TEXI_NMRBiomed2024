function s = expm_fit2data(m, xps)

k  = m(1);
f1 = m(2);
D1 = m(3);
D2 = m(4);

f = [f1 1-f1]';
D = [D1 0 ; 0 D2];

k12 = k * (1-f1);
k21 = k * f1;
K = [k12 -k21 ; -k12 k21];

NT = size(xps.q2,1);


dt = xps.dt;

if (0)
for nwfm = 1:xps.n
    S = f;
    for nt = 1:NT
        S = expm(-(xps.q2(nt,nwfm)*D + K) * dt) * S;
    end
    s(nwfm) = sum(S);
end

else

I = ones(size(xps.q2(1,:)));



S = f;

S1 = f(1) * I;
S2 = f(2) * I;
for nt = 1:NT
    d11 = (-xps.q2(nt,:)*D1 - k12)*dt;
    d12 = k21*dt*I;    
    d21 = k12*dt*I;
    d22 = (-xps.q2(nt,:)*D2 - k21)*dt;

    dS1 = d11.*S1 + d12.*S2;
    dS2 = d21.*S1 + d22.*S2;

    S1 = S1 + dS1;
    S2 = S2 + dS2;
end

s = S1+S2;
end




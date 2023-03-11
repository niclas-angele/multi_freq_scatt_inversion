
function u=solv(h,xbis,kk,grilx)
    y=xbis;
    r=max(xbis)-min(xbis);
    x=[0];
    M=x'*ones(1,grilx)-ones(1,1)*y;
    N=G(M,kk);
    A=N*(h(xbis).');
    u=r*A/grilx;
end

function y=G(x,kk)
    y=1i/2/kk.*exp(1i*kk*abs(x));
end
function s=ff(S,k,y,pas)
%%%% function which calculates the Fourier data of a source S given in 
% the form of a vector with its values on y, at the frequency k and it 
% is also necessary to give the discretization step of y
    M=(ones(length(y),1)*k).*(y'*ones(1,length(k)));
    V=exp(1i*M);
    s=(S*V).*(1i/2*pas);
end
function s=solve1Dcarreb(d,k,a,b,gril2,lambda)
%%% inverse fourier of d at frequency k with a penalization lambda. give
%%% the initial function on the discretization a:gril2:b
    x=linspace(a,b,gril2);
    pas=(b-a)/(gril2-1);
    S=0*x;
    D=Fadj(d,k,x,pas);
    for i=1:1000
        FS=F(S,k,x,pas);
        grad=Fadj(FS,k,x,pas)-D+lambda*Gadj(G(S));
        Fgrad = F(grad,k,x,pas);
        Ggrad =G(grad);
        rho=norm(grad)^2/(norm(Fgrad)^2+norm(Ggrad)^2);
        S=S-rho*grad;
        S=real(S);
        c=norm(F(S,k,x,pas)-d)^2+lambda*norm(G(S))^2;
        if c<10^(-10)
            break
        end
    end
    s=S;
end

function s=G(u)
    s=u-circshift(u,-1);
end 

function s=Gadj(u)
    s=u-circshift(u,1);
end

function s=Fadj(u,k,y,pas)
    M=(ones(length(k),1)*y).*(k'*ones(1,length(y)));
    V=exp(-1i*M);
    S=-1i*pas/2.*u;
    s=S*V;
end
 
 
function s=F(S,k,y,pas)
    M=(ones(length(y),1)*k).*(y'*ones(1,length(k)));
    V=exp(1i*M);
    s=(S*V).*(1i/2*pas);
end
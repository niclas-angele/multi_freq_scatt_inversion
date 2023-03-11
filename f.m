%%% fonction qui calcule s_n en donnant la fonction s, le mode n et la
%%% grille x sur laquelle on veut faire le calcul. 

function temp=f(x,n,s)
    temp=0*x;
    for i=0:99
        y=0.01*i.*ones(1,length(x));
        temp=temp+0.01*s(x,y).*cos(n*pi.*y);
    end
    if n>0
        temp=sqrt(2).*temp;
    end
end
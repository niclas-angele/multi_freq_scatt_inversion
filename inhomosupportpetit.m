clear; close all;
%%% choice of the inhomogeneity 
% h=@(x,y) (1.*(((x-4)/0.2).^2+((y-0.7)/0.15).^2)<1).*(1-(((x-4)/0.2).^2+((y-0.7)/0.15).^2))+(1.*(((x-4.5)/0.1).^2+((y-0.25)/0.1).^2)<1).*(1-(((x-4.5)/0.1).^2+((y-0.25)/0.1).^2))+(1.*(((x-5.2)/0.15).^2+((y-0.5)/0.1).^2)<1).*(1-(((x-5.2)/0.15).^2+((y-0.5)/0.1).^2));
%h=@(x,y) (1.*(((x-5)/0.2).^2+((y-0.5)/0.1).^2)<1).*(1-(((x-5)/0.2).^2+((y-0.5)/0.1).^2));
h=@(x,y) 0.05*(1.*(((x-4)/0.05).^2+((y-0.6)/0.15).^2)<1).*(1-(((x-4)/0.05).^2+((y-0.6)/0.15).^4));


%considered frequencies 
miin=0.01;
maax=150;
grilk=200;
k=linspace(miin,maax,grilk);
%maximum number of modes considered 
N=min(floor(sqrt(maax^2/pi^2)),20);


%domain 
nodes2 = [-1 0; 3 0; 6 0; 10 0;10 1;6 1; 3 1; -1 1];
edges2 = {{1,2}, {2,3}, {3,4}, {4,5}, {5,6}, {6,7}, {7,8}, {8,1},{2,7,'lr'},{3,6,'lr'}};
dx = 0.2; %pas
domain2 = Domain(nodes2,edges2);
mesh2 = Mesh(domain2,dx,'save');
mesh2=mesh2.submesh;
mesh2=mesh2.submesh;
%refine mesh to improve results 
% mesh2=mesh2.submesh;
% mesh2=mesh2.submesh;
% mesh2=mesh2.submesh;
% mesh2=mesh2.submesh;
%measurement point 
x_0=3;


%%% spactial parameters
a=3.8;
b=4.2;
grilx=max(floor(1000*(b-a)),100);
grily=100;
x=linspace(a,b,grilx);
y=linspace(0,1,grily);
[X,Y]=meshgrid(x,y);


%matrix containing the defect
xx = mesh2.tricenters(:,1);
yy = mesh2.tricenters(:,2);

%storage data
DONNEE=[];
%PML parameters 
alph=@(k) k.*(k>5)+1*(k<=5);
for i=1:length(k)
    kk=k(i);
    c_dampb = -1i*alph(kk)*kk*mesh2.P0( '1*(x>6.5).*(x-6.5) + 1*(x<2.5).*(-x+2.5)'); 
    bc1 = {[4,8],1,-1i*kk,0};
    bc2 = {[5,6,7],1,0,0};
    bc3 = {[1,2,3],1,0,0};
    v = mesh2.solve(1,-kk^2*(1+h(xx,yy))+c_dampb,kk^2*h(xx,yy).*exp(1i*kk*xx),{bc1,bc2,bc3});
    %measurements 
    I = mesh2.boundary(9);
    xgauche = mesh2.nodes(I,:);
    xgauche = xgauche(:,2);
    vgauche=v(I);
    G=sortrows([xgauche vgauche]);
    DONNEE=[DONNEE,G];
end

%computation of the Fourier transform 
Z=X*0+Y*0;
for n=0:10
    kn=real(sqrt(k.^2-n^2*pi^2));
    freq=[];
    U=[];
    for i=1:length(k)
        xsi=kn(i); 
        if xsi~=0 %only propagative modes 
        G=DONNEE(:,2*i-1:2*i);
        temp1=G(:,1);
        temp2=G(:,2);
        %computation of u_n 
        q=0;
        for j=1:length(temp1)-1
            q=q+(temp1(j+1)-temp1(j))*temp2(j)*cos(n*pi*temp1(j));
        end
        if n>0
            q=sqrt(2)*q;
        end       
        U=[U,xsi/k(i)/k(i)*exp(1i*xsi*x_0)*q]; %value of u_n(0)
        freq=[freq,k(i)+xsi]; %value of k_n+k
        end
    end
    %computation of h_n
    z=solve1Dcarreb(U,freq,a,b,grilx,0.002);
    %modal sum 
    for i=1:grilx
        for j=1:grily
            temp=z(i)*cos(n*pi*y(j)); 
            if n>0
                temp=sqrt(2).*temp;
            end
            Z(j,i)=Z(j,i)+temp;
        end
    end
end


ZZ=0*X;
for n=0:N %reste des modes
    V=[];
    freq2=[];
    kn=real(sqrt(k.^2-n^2*pi^2));
    hh=@(x) f(x,n,h); %mode n de la source
    freq2=[]; %init des freq
    for i=1:length(k)
        xsi=kn(i); 
        if xsi~=0 %on ne regarde que les k_n positifs (modes propagatifs)
            s1= @(x) k(i)^2*hh(x).*exp(1i*k(i)*x); 
            V=[V,xsi/k(i)/k(i)*(solv(s1,x,xsi,grilx))]; %valeur de u_n(0)
            freq2=[freq2,k(i)+xsi]; %valeur de k_n+k
        end
    end
    zz=solve1Dcarreb(V,freq2,a,b,grilx,0.08); %moindres carrés pour trouver la forme de h_n
    for i=1:grilx %reconstruction de la somme des modes
        for j=1:grily
            temp=zz(i)*cos(n*pi*y(j)); 
            if n>0
                temp=sqrt(2).*temp;
            end
            ZZ(j,i)=ZZ(j,i)+temp;
        end
    end
end


subplot 211
temp=1.*h(X,Y);
surf(X,Y,temp,'edgecolor','none');
view(0,90)
xlim([3.8,4.2])
colormap;
colorbar
caxis([-0.01,0.07]);
title('h')

subplot 212
surf(X,Y,real(Z),'edgecolor','none');
view(0,90)
colormap;
xlim([3.8,4.2])
colorbar
title('reconstruction of h kmax=350')
caxis([-0.01,0.07]);






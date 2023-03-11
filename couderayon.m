clear; close all;
%parameters of the bend 
theta=pi/50;
r=40;
xc=3;
%needed function to identify r 
amplitude = @(y) (y-1).*(1./(r+y)+1./(r+1)); 

%%%% domain 
nodes2 = [-19 0; 1 0; 7 0;27  0;27 1;7 1; 1 1; -19 1];
edges2 = {{1,2}, {2,3}, {3,4}, {4,5}, {5,6}, {6,7}, {7,8}, {8,1},{2,7,'lr'},{3,6,'lr'}};
dx = 0.3; %pas
domain2 = Domain(nodes2,edges2);
mesh2 = Mesh(domain2,dx,'save');
mesh2=mesh2.submesh;
mesh2=mesh2.submesh;
%%%improve the reconstruction by increasing the mesh 
% mesh2=mesh2.submesh;
% mesh2=mesh2.submesh;
% mesh2=mesh2.submesh;

%used frequencies 
miin=0.01;
maax=40;
grilk=100;
k=linspace(miin,maax,grilk);

%%%spatial parameters 
a=1;
b=7;
x_0=1;
grilx=floor(100);
grily=100;
x=linspace(a,b,grilx);
y=linspace(0,1,grily);
[X,Y]=meshgrid(x,y);

%storage matrices 
Z=X*0+Y*0;
freq=[];
U=[];


%matrices to include the defect 
xx = mesh2.tricenters(:,1);
yy = mesh2.tricenters(:,2);
alph=@(k) (k.*(k>5)+(k<=5));
S = [(xx<xc+theta*(r+1)).*(xx>xc).*(r+1)./(r+yy)+(xx>xc+theta*(r+1))+(xx<xc) 0*xx 0*xx (xx<xc+theta*(r+1)).*(xx>xc).*(r+yy)./(r+1)+(xx>xc+theta*(r+1))+(xx<xc)];
d=(xx>xc).*(xx<xc+theta*(r+1)).*(r+yy)/(r+1)+(xx>xc+theta*(r+1))+(xx<xc);

for i=1:length(k)
    kk=k(i);
    c_dampb = -1i*alph(kk)*kk*mesh2.P0( '1*(x>8).*(x-8) + 1*(x<0).*(-x)');
    h = @(x,y) (x<xc+theta*(r+1)).*(x>xc).*kk^2.*exp(1i*kk*x).*(y-1).*(1./(r+y)+1./(r+1));
    bc3 = {[4,8],1,-1i*kk,0};
    w = mesh2.solve(S,-kk^2*d+c_dampb,mesh2.P0(h),{{[1,2,3,5,6,7],1,0,0},bc3}); %apres Born
    I = mesh2.boundary(9);
    xgauche = mesh2.nodes(I,:);
    xgauche = xgauche(:,2);
    vgauche=w(I);
    %measurements at the left of the waveguide 
    G=sortrows([xgauche vgauche]);
    temp1=real(G(:,1));
    temp2=G(:,2);
    %computation of u_0
    p=0;
    for j=1:length(temp1)-1
        p=p+(temp1(j+1)-temp1(j))*temp2(j);
    end
    U=[U,exp(1i*kk*x_0)/kk*p]; %value of u_0(0)
    freq=[freq,2*kk]; %value of k_n+k
end

%inverse fourier transform with penalized algorithm
zz=solve1Dcarre(U,freq,a,b,grilx); 
z=abs(U);
t=imag(U);

%identification of the three parameters
g=@(p) norm(zz-(1-1/2/(p(2)+1)-(p(2)+1)*log((p(2)+1)/p(2)))*(x>p(3)).*(x<p(3)+(p(2)+1)*p(1)))^2; 
p=fminsearch(g,[theta,r,xc]); 
plot(x,zz)
hold on
plot(x,(1-1/2/(p(2)+1)-(p(2)+1)*log((p(2)+1)/p(2)))*(x>p(3)).*(x<p(3)+(p(2)+1)*p(1)))

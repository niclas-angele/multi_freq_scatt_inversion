clear; close all ; 

%domain 
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

%shape of the defect 
amp=1;
h=@(x) amp*5/4/4*(x>3.2).*(x<4.2).*(x-3.2).^2.*(x-4.2).^2+1;
hprim=@(x) amp*5/2/4*(x>3.2).*(x<4.2).*(x-3.2).*(x-4.2).*(2*x-7.4);
g=@(x) -7*amp*5/4/4*(x>3.4).*(x<4).*(x-3.4).^2.*(x-4).^2;
gprim=@(x) -7*amp*5/2/4*(x>3.4).*(x<4).*(x-3.4).*(x-4).*(2*x-7.4);

%parameters
a=3; %left of the waveguide
b=4.5; %right of the waveguide 
miin=0.01; %min freq 
maax=70; %max freq 
grilk=300; %discretization
k=linspace(miin,maax,grilk); %frequencies 
grilx=floor(100*(b-a)); %spatial discretization
x=linspace(a,b,grilx); 

%storage variables 
V=[];
W1=[];
freq=[];
freq1=[];
%coordinate of the measurement section
x_0=1;

%matrix to include the change of variable 
xx = mesh2.tricenters(:,1);
yy = mesh2.tricenters(:,2);
%parameters of the PML
alph=@(k) (k.*(k>8)+(k<=8));
S = [h(xx)-g(xx), -(hprim(xx)-gprim(xx)).*yy+gprim(xx) , -(hprim(xx)-gprim(xx)).*yy-gprim(xx), ((hprim(xx)-gprim(xx)).*yy+gprim(xx)).^2./(h(xx)-g(xx))+1./(h(xx)-g(xx))];
d=h(xx)-g(xx);
for i=1:length(k)
    kk=k(i);
    %computation of the wavefield 
    c_dampb = -1i*alph(kk)*kk*mesh2.P0( '1*(x>8).*(x-8) + 1*(x<0).*(-x)'); 
    bord_u0=@(x,y) 1i*kk*exp(1i*kk*x).*hprim(x);
    bord_u2=@(x,y) 1i*kk*exp(1i*kk*x).*gprim(x);
    bc1 = {[4,8],1,-1i*kk,0};
    bc2 = {[5,6,7],1,0,bord_u0};
    bc3 = {[1,2,3],1,0,bord_u2};
    v = mesh2.solve(S,-kk^2*d+c_dampb,0,{bc1,bc2,bc3}); %resolution
    I = mesh2.boundary(9);
    xgauche = mesh2.nodes(I,:);
    xgauche = xgauche(:,2);
    vgauche=v(I);
    %measurements at the left of the waveguide 
    G=sortrows([xgauche vgauche]); 
    temp1=real(G(:,1));
    temp2=G(:,2);
    %computation of u_0
    p=0;
    for j=1:length(temp1)-1
        p=p+(temp1(j+1)-temp1(j))*temp2(j);
    end
    %computation of u_1
    q=0;    
    for j=1:length(temp1)-1
        q=q+(temp1(j+1)-temp1(j))*temp2(j)*cos(pi*temp1(j));
    end
    %matrix containing the Fourier transform
    k_1=sqrt(kk^2-pi^2);
    if (2.6>kk) ||(kk>3.5) && (kk<6) || (kk>6.5) && (kk<9) || (kk>11)
        V=[V,-1i*exp(1i*kk*x_0)*p];
        freq=[freq,2*kk]; 
    end
    if (kk>3.5) && (kk<6) || (kk>6.5) && (kk<9) || (kk>11)
        W1=[W1,-1i*k_1/kk*exp(1i*k_1*x_0)*q];
        freq1=[freq1,kk+k_1];
    end    
end

%inverse fourier transform with penalized algorithm
zz=solve1Dcarreb(V,freq,a,b,grilx,0.08);
tt=solve1Dcarreb(W1,freq1,a,b,grilx,0.08);

%reconstruction of h and g derivatives
hhprim=(zz-tt)/2;
ggprim=(zz+tt)/2;
hhprim=hhprim-mean(hhprim);
ggprim=ggprim-mean(ggprim);

%recontruction of h ang g
hhg=0*hhprim;
hhd=0*ggprim;

for i=2:length(zz)
    hhg(i)=hhg(i-1)+hhprim(i)*(x(i)-x(i-1));
    hhd(i)=hhd(i-1)+ggprim(i)*(x(i)-x(i-1));
end

U=ff(hprim(x)+gprim(x),freq,x,(b-a)/(grilx-1));
U2=ff(hprim(x)-gprim(x),freq1,x,(b-a)/(grilx-1));

%%% plots 
subplot 221
plot(x,hprim(x))
hold on 
plot(x,hhprim)
title("h'")
legend("h'","h'_{app}")
subplot 223
plot(x,gprim(x))
hold on 
plot(x,ggprim)
title("g'")
legend("g'","g'_{app}")
subplot 222
plot(x,h(x)-1)
hold on 
plot(x,hhg)
title('h')
legend('h','h_{app}')
subplot 224
plot(x,g(x))
hold on 
plot(x,hhd)
title('g')
legend('g','g_{app}')



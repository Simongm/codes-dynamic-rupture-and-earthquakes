function dynamic_fem


% %% initialize movie:
% close all;figure
% myVideo = VideoWriter('dynamic.avi'); % create object
% myVideo.FrameRate = 30;  % frame rate Default 30
% open(myVideo);


%% define constants:
G=3e9; %shear modulus
nu=0.35; %poisson ratio
rho=3500; %density
E=2*G*(1+nu);
vs=(G/rho)^0.5;
vp=((2*G*(1-nu))/(rho*(1-2*nu)))^0.5;
lambda=(2*G*nu)/(1-2*nu);
% or lambda2=vp^2*rho-2*G;


%% construct the grid
nz=600;nx=600; %nb of elements
x=linspace(0,0.3,nx);
z=linspace(0,0.3,nz);
dx=abs(x(2)-x(1));
dz=abs(z(2)-z(1));
[xx,zz]=meshgrid(x,z);
dt=1/(vp*sqrt(1/dx^2+1/dz^2)); %timestep according to courant stability criterion
% nt=200; %nb of iterations
% t=0:nt;t=t*dt;
% tmax=dt*nt;
tmax=0.002;
t=linspace(0,tmax,round(tmax/dt));nt=length(t);
%% initializes the vectors
vx_old=zeros(size(xx));
vz_old=zeros(size(xx));
txx_old=zeros(size(xx));
tzz_old=zeros(size(xx));
txz_old=zeros(size(xx));
vx_new=vx_old;
vz_new=vz_old;
txx_new=txx_old;
tzz_new=tzz_old;
txz_new=txz_old;


%% source terms
% amplitude of gaussian pulse for point source
f0=50/(tmax); 
t0=2/f0;
amp=1*exp(-f0^2*(t-t0).^2); amp=amp-amp(1);


%% prestress to start dynamic rupt
size_patch=0.008; % size half patch in m
nb_pt_patch=round(size_patch/dx);
nxs=round(nx/2);nzs=round(nz/2);Dc=0.1e-5;tau_0=1.5e6;tau_p=1.8e6;
rup=zeros(1,nx);posrup=zeros(nt,nx-1);
prestress(1:nx)=tau_0;prestress(nxs-nb_pt_patch:nxs+nb_pt_patch)=1.800*10^6;
u_interface_old=zeros(1,nx);
tau_f=0.9e6;
cc=hot(nt);
figure;
% hold on;
stress_init=zeros(size(xx))+prestress(1,1);
colors = jet(600);
for i=1:600;%nt
    loading_rate=50e6;
    color = colors(i,:);
%     txz_old=txz_old+loading_rate*dt;
    [vx_new,vz_new,txx_new,tzz_new,txz_new]=tnew_virieux(xx,zz,vx_old,vz_old,txx_old,tzz_old,txz_old,vx_new,vz_new,txx_new,tzz_new,txz_new,lambda,G,dt,rho);

    % src term for kinematic ruputre
%     vx_new(round(length(x)/2),round(length(z)/2))=amp(i);
%     tyz_new(100,:)=amp(i);


    %% or spontaneous rupture
    u_interface_new=u_interface_old+vx_new(nzs,:)*dt;
    rup(rup==0 & (abs(txz_new(nzs,:)+prestress)>tau_p))=1;
    posrup(i,:)=diff(rup,1,2);
    txz_new(nzs,rup==1)=tau_p-(tau_p-tau_f)/Dc*u_interface_new(rup==1)-prestress(rup==1);
    txz_new(nzs,txz_new(nzs,:)<(tau_f-tau_0))=tau_f-tau_0;
    u_interface_old=u_interface_new;
    %to increase the loading 
    tau_p=tau_p-loading_rate*dt;
    tau_f=tau_f-loading_rate*dt;
    
     %% plots
%for stress at interface-------------
    if rem(i*dt,30*dt)==0
        plot(xx(nzs,:),tau_0+txz_new(nzs,:)+(i-1)*loading_rate*dt,'color', color)%,'color',cc(round(i*nt/500),:))
        ylim([tau_f,tau_p])
%         stress2plt=txz_new+stress_init+(i-1)*loading_rate*dt;
%         surf(xx,zz,stress2plt,'edgecolor','none')
        title(['t=',num2str(t(i)),'s'])
%         view(2)
%         colorbar()
%     zlim([0.9e6+(i-1)*loading_rate*dt, 1.9e6])
%     caxis([1.4e6+(i-1)*loading_rate*dt, 1.7e6])
%     xlabel('x (m)')
%     ylabel('y (m)')
%     zlabel('stress (MPa)')
        drawnow()
%         writeVideo(myVideo, getframe); 
        hold on
    end
%----------------------------------

    
    
    %%update vectors
    vx_old=vx_new;
    vz_old=vz_new;
    txx_old=txx_new;
    tzz_old=tzz_new;
    txz_old=txz_new;
    
    
    
end
% 
figure(2)
% surf(posrup,'edgecolor','none');view(2)
[t1,x1]=find(posrup==-1);
[t2,x2]=find(posrup==1);
plot(t1,x1,t2,x2)
% figure(3)
% plot(t1,(x1-x2))
% close(myVideo)

function [vx_new,vz_new,txx_new,tzz_new,txz_new]=tnew_virieux(xx,zz,vx_old,vz_old,txx_old,tzz_old,txz_old,vx_new,vz_new,txx_new,tzz_new,txz_new,lambda,G,dt,rho)

%for velocities we use the old stress variables
dtxxdxx=diff(txx_old,1,2)./diff(xx,1,2);dtxxdxx=dtxxdxx(1:end-1,:);dtxzdzz=diff(txz_old,1,1)./diff(zz,1,1);dtxzdzz=dtxzdzz(:,2:end);
vx_new(1:end-1,2:end)=vx_old(1:end-1,2:end)+dt/rho*(dtxxdxx+dtxzdzz);

dtxzdxx=diff(txz_old,1,2)./diff(xx,1,2);dtxzdxx=dtxzdxx(2:end,:);dtzzdzz=diff(tzz_old,1,1)./diff(zz,1,1);dtzzdzz=dtzzdzz(:,1:end-1);
vz_new(2:end,1:end-1)=vz_old(2:end,1:end-1)*dt/rho*(dtxzdxx+dtzzdzz);


% %% sponge boundaries  
%  for j=1:nz;for i=1:nn;sx(j,i)=sx(j,i)*(ax1*i*dx+cx1);end;end;
%  for j=1:nz;for i=nx-nn:nx;sx(j,i)=sx(j,i)*(ax2*i*dx+cx2);end;end;
%  for j=nz-nn:nz;for i=1:nx;sz(j,i)=sz(j,i)*(az*j*dx+cz);end;end;
 
 
%for stresses we use the newly calculated velocities
dvxdxx= diff(vx_new,1,2)./diff(xx,1,2);dvxdxx=dvxdxx(1:end-1,:);dvzdzz=diff(vz_new,1,1)./diff(zz,1,1);dvzdzz=dvzdzz(:,1:end-1);
txx_new(1:end-1,1:end-1)=txx_old(1:end-1,1:end-1)+dt*((lambda+2*G)*dvxdxx+lambda*dvzdzz);
tzz_new(1:end-1,1:end-1)=tzz_old(1:end-1,1:end-1)+dt*(lambda*dvxdxx+(lambda+2*G)*dvzdzz);

dvxdzz=diff(vx_new,1,1)./diff(zz,1,1);dvxdzz=dvxdzz(:,2:end);dvzdxx=diff(vz_new,1,2)./diff(xx,1,2);dvzdxx=dvzdxx(2:end,:);
txz_new(2:end,2:end)=txz_old(2:end,2:end)+dt*G*(dvxdzz+dvzdxx);

%% impose BC

% vx_new(1,:)=0;vx_new(:,1)=0;vx_new(end,:)=0;vx_new(:,end)=0;
% vz_new(1,:)=0;vz_new(:,1)=0;vz_new(end,:)=0;vz_new(:,end)=0;

txx_new(1,:)=0;txx_new(:,1)=0;txx_new(end,:)=0;txx_new(:,end)=0;
tzz_new(1,:)=0;tzz_new(:,1)=0;tzz_new(end,:)=0;tzz_new(:,end)=0;
txz_new(1,:)=0;txz_new(:,1)=0;txz_new(end,:)=0;txz_new(:,end)=0;
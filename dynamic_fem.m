function dynamic_fem


% %% initialize movie:
% close all;figure
% myVideo = VideoWriter('dynamic_supshear.avi'); % create object
% myVideo.FrameRate = 30;  % frame rate Default 30
% open(myVideo);


%% define constants:
G=30e9; %shear modulus
nu=0.25; %poisson ratio
rho=3000; %density
E=2*G*(1+nu);
vs=(G/rho)^0.5;
vp=((2*G*(1-nu))/(rho*(1-2*nu)))^0.5;
lambda=(2*G*nu)/(1-2*nu);
% or lambda2=vp^2*rho-2*G;


%% construct the grid
% nz=3000;nx=3000; %nb of elements
nx=300; %nb of elements
x=linspace(0,20*10^3,nx);
nz=300;
z=linspace(0,20*10^3,nz);
dx=abs(x(2)-x(1));
dz=abs(z(2)-z(1));
[xx,zz]=meshgrid(x,z);
dt=0.8/(vp*sqrt(1/dx^2+1/dz^2)); %timestep according to courant stability criterion

%% % Attenuation at the boundaries
% Attnx=ones(1,length(x));
% width_Bound=max(x)/5;
% Attnx(x<width_Bound)=(sin(pi().*x(x<width_Bound)/width_Bound-pi()/2)+0.5)/4;
% Attnx(x>(max(x)-width_Bound))=(sin(pi()*(x(x>(max(x)-width_Bound))-(2*max(x)-width_Bound))/width_Bound-pi()/2)+1)/2;
% plot(x,Attnx)
% % 
% Attnz=ones(1,length(z));
% width_Bound=max(z)/5;
% Attnz(z<width_Bound)=(sin(pi().*z(z<width_Bound)/width_Bound-pi()/2)+1)/2;
% Attnz(z>(max(z)-width_Bound))=(sin(pi()*(z(z>(max(z)-width_Bound))-(2*max(z)-width_Bound))/width_Bound+pi()/2)+1)/2;
% % c
% Attn=zeros(nz,nx);
% for zit=1:nz
%     Attn(zit,:)=Attnx;
% end
% 
% for xit=1:nx
%     size(Attn(:,xit))
%     size(Attnz)
%     Attn(:,xit)=Attn(:,xit).*Attnz';
% end

% surf(xx,zz,Attn)
% % nt=200; %nb of iterations
% % t=0:nt;t=t*dt;
% % tmax=dt*nt;
%% end attenuation box %% 


tmax=4;
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
size_patch=1.2*10^3; % size half patch in m
nb_pt_patch=round(size_patch/dx);
nxs=round(nx/2);nzs=round(nz/2);Dc=0.05;tau_0=3.6e6;tau_p=5e6;
rup=zeros(1,nx);posrup=zeros(nt,nx-1);
prestress(1:nx)=tau_0;prestress(nxs-nb_pt_patch:nxs+nb_pt_patch)=tau_p*1.2;
u_interface_old=zeros(1,nx);
tau_f=1e6;
cc=hot(nt);
figure;
% hold on;
stress_init=zeros(size(xx))+prestress(1,1);
tic
for i=1:nt
    loading_rate=0.1e6;
    [vx_new,vz_new,txx_new,tzz_new,txz_new]=tnew_virieux(xx,zz,vx_old,vz_old,txx_old,tzz_old,txz_old,vx_new,vz_new,txx_new,tzz_new,txz_new,lambda,G,dt,rho);

    % src term for kinematic ruputre
%    vx_new(round(length(x)/2),round(length(z)/2))=vx_new(round(length(x)/2),round(length(z)/2))+amp(i);
%     vz_new(round(length(x)/2),round(length(z)/2))=vz_new(round(length(x)/2),round(length(z)/2))+amp(i);
 %    txx_new(round(length(x)/2),round(length(z)/2))=txx_new(round(length(x)/2),round(length(z)/2))+1e-10*amp(i);
 %    tzz_new(round(length(x)/2),round(length(z)/2))=tzz_new(round(length(x)/2),round(length(z)/2))+1e-10*amp(i);

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
%         
%      %% plots
%for stress at interface-------------
    if rem(i*dt,1*dt)==0 % to have a plot every x iterations
    toc;
    fprintf([num2str(toc*(nt-i)/100/3600) 'hr ' 'left'])
%         plot(xx(nzs,:),tau_0+txz_new(nzs,:)+(i-1)*loading_rate*dt)%,'color',cc(round(i*nt/500),:))
%         ylim([tau_f,tau_p])
%         stress2plt=sqrt(vx_new.^2+vz_new.^2);
        stress2plt=txz_new+stress_init+(i-1)*loading_rate*dt;
        surf(xx,zz,stress2plt,'edgecolor','none')
        title(['t=',num2str(t(i)),'s'])
        view(2)
%         colorbar()
%     zlim([0.9e6+(i-1)*loading_rate*dt, 1.9e6])
      caxis([tau_f,tau_p])
      axis equal
%     xlabel('x (m)')
%     ylabel('y (m)')
%     zlabel('stress (MPa)')
      drawnow()
%         writeVideo(myVideo, getframe); 
%         hold on
    tic
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
% figure(2)
% surf(posrup,'edgecolor','none');view(2)
[it1,ix1]=find(posrup==-1);
[it2,ix2]=find(posrup==1);
x1=x(ix1);x2=x(ix2);t1=t(it1);t2=t(it2);
% plot(t1,x1,'.',t2,x2,'.')

% interpolate the data on a common time axis
index=1:10:length(t1); %index to downsampled the data
t1=t1(index);t2=t2(index);x2=x2(index);x1=x1(index);
T=linspace(min(t1),max(t1),length(t1)); X1=interp1(t1,x1,t1);X2=interp1(t2,x2,t1);
figure(3)
vel_crack=(diff(X1-X2)./(T(2)-t(1)))/2;
plot(T(1:end-1),vel_crack)
m=zeros(2,length(X1)-1);m(1,:)=T(1:end-1);m(2,:)=vel_crack; %there is a factor 2 because we take the crack expansion velocity
% close(myVideo)
figure(4)
plot(T,X1,'.',T,X2,'.')
csvwrite([num2str(loading_rate/1e6) 'mpas' '.dat'],m)

function [vx_new,vz_new,txx_new,tzz_new,txz_new]=tnew_virieux(xx,zz,vx_old,vz_old,txx_old,tzz_old,txz_old,vx_new,vz_new,txx_new,tzz_new,txz_new,lambda,G,dt,rho)

%for velocities we use the old stress variables
dtxxdxx=diff(txx_old,1,2)./diff(xx,1,2);dtxxdxx=dtxxdxx(1:end-1,:);dtxzdzz=diff(txz_old,1,1)./diff(zz,1,1);dtxzdzz=dtxzdzz(:,2:end);
vx_new(1:end-1,2:end)=vx_old(1:end-1,2:end)+dt/rho.*(dtxxdxx+dtxzdzz);

dtxzdxx=diff(txz_old,1,2)./diff(xx,1,2);dtxzdxx=dtxzdxx(2:end,:);dtzzdzz=diff(tzz_old,1,1)./diff(zz,1,1);dtzzdzz=dtzzdzz(:,1:end-1);
vz_new(2:end,1:end-1)=vz_old(2:end,1:end-1)+dt/rho.*(dtxzdxx+dtzzdzz);

%for stresses we use the newly calculated velocities
dvxdxx= diff(vx_new,1,2)./diff(xx,1,2);dvxdxx=dvxdxx(1:end-1,:);dvzdzz=diff(vz_new,1,1)./diff(zz,1,1);dvzdzz=dvzdzz(:,1:end-1);
txx_new(1:end-1,1:end-1)=txx_old(1:end-1,1:end-1)+dt*((lambda+2*G)*dvxdxx+lambda*dvzdzz);
tzz_new(1:end-1,1:end-1)=tzz_old(1:end-1,1:end-1)+dt*(lambda*dvxdxx+(lambda+2*G)*dvzdzz);

dvxdzz=diff(vx_new,1,1)./diff(zz,1,1);dvxdzz=dvxdzz(:,2:end);dvzdxx=diff(vz_new,1,2)./diff(xx,1,2);dvzdxx=dvzdxx(2:end,:);
txz_new(2:end,2:end)=txz_old(2:end,2:end)+dt*G.*(dvxdzz+dvzdxx);

%% impose BC

vx_new(1,:)=0;vx_new(:,5)=0;vx_new(end,:)=0;vx_new(:,end-5)=0;
vz_new(1,:)=0;vz_new(:,5)=0;vz_new(end,:)=0;vz_new(:,end-5)=0;
txx_new(5,:)=0;txx_new(:,1)=0;txx_new(end-5,:)=0;txx_new(:,end)=0;
tzz_new(5,:)=0;tzz_new(:,1)=0;tzz_new(end-5,:)=0;tzz_new(:,end)=0;
txz_new(5,:)=0;txz_new(:,1)=0;txz_new(end-5,:)=0;txz_new(:,end)=0;


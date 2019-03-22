import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from copy import deepcopy
##-------------------------------PRE PROCESSING---------------------------------##
## meshing --------------------------------------------------------
class setup:
    def __init__(self,nelx,nely,width,heigth,E,nu,densty):
        self.el_width= width/nelx 
        self.el_heigth= heigth/nely
        self.heigth=heigth
        self.width=width
        self.E=E
        self.nu=nu
        self.densty=densty
        self.nelx=nelx
        self.nely=nely
        self.nnodes=(nelx+1)*(nely+1)
        self.nels=nelx*nely
        self.p=10*np.max(np.array(E))/(self.el_width*(1-np.min(np.array(nu))**2))  
        
        self.coordx=np.linspace(0,width,nelx+1)
        self.coordy=np.insert(np.linspace(0,heigth,nely+1), int(nely/2)+1, heigth/2)
        coord=np.zeros(((nelx+1)*(nely+2),2))
        etpl=np.zeros(((nelx)*(nely),4))
        all_dofs= np.arange((nelx+1)*(nely+2)*2)
        
        # calculate coordinates of global nodes :
        idx=0
        for j in range(nely+2):
            for i in range(nelx+1):
                coord[idx,:]=np.array([self.coordx[i],self.coordy[j]])
                idx=idx+1        
        # calculate connectivity matrix etpl :
        idx=0
        etplall_bot=np.arange((1+int(nely/2))*(1+nelx)).reshape(1+int(nely/2),nelx+1)
        etplall_top=np.arange((1+int(nely/2))*(1+nelx)).reshape(1+int(nely/2),nelx+1)+np.max(etplall_bot)+1
        for i in range(int(nely/2)):
            for j in range(nelx):        
                etpl[idx]=np.array([etplall_bot[i,j],etplall_bot[i,j+1],etplall_bot[i+1,j+1],etplall_bot[i+1,j]])
                idx=idx+1

        for i in range(int(nely/2)):
            for j in range(nelx):        
                etpl[idx]=np.array([etplall_top[i,j],etplall_top[i,j+1],etplall_top[i+1,j+1],etplall_top[i+1,j]])
                idx=idx+1   
        etpl=etpl.astype(int)

        ## material properties ##################################
        matnodes=np.zeros((len(coord),1)) #--> E, nu
        matelems=np.zeros((len(etpl),1))#--> E, nu, rho
        coordxc=(coord[(etpl[:,0]),0]+coord[(etpl[:,3]),0])/2
        coordyc=(coord[(etpl[:,0]),1]+coord[(etpl[:,3]),1])/2
        ## nodes corresponding to steel holders: ## set to 0 
        matnodes[:,:]=1
        matelems[:,:]=1
        ## nodes corresponding to granite: ## rmove - 1 and 1
        x1=0.015-1
        x2=0.055+1
        matnodes[(coord[:,0]>x1) & (coord[:,0]<x2) & (coord[:,1]>0.006) & (coord[:,1]<0.014),:]=1
        matelems[(coordxc>x1) & (coordxc<x2) & (coordyc>0.006) & (coordyc<0.014),:]=1
        ## nodes corresponding to silicon spacers: ## set to 2
        matnodes[((coord[:,0]>0.005) & (coord[:,0]<x1) & (coord[:,1]>0.01015)) | \
                 ((coord[:,0]>=x2) & (coord[:,0]<0.065) & (coord[:,1]<0.00985))] \
                = 1
        matelems[((coordxc>0.005) & (coordxc<x1) & (coordyc>0.01015)) | \
                 ((coordxc>=x2) & (coordxc<0.065) & (coordyc<0.00985))] \
                = 1
        ##_END_Material properties_##############################

        self.etpl=etpl
        self.coord=coord
        self.all_dofs=all_dofs
        int_nodes=np.where(np.logical_and(coord[:,1]==(heigth/2),np.transpose(matnodes[:,0])==1))[0]
        dofs_int=np.append(int_nodes*2,int_nodes*2+1).reshape(4,int(len(int_nodes)/2))
        self.dofs_int=dofs_int        
        self.int_nodes=int_nodes
        self.coordint=(self.coordx)[np.where(np.logical_and(self.coordx>x1,self.coordx<x2))[0]]
        self.allgp=(1/3)**0.5*np.array([[-1,1,1,-1],[-1,-1,1,1]])
        self.J=np.array([[(0.5*width/nelx)/1,0],[0,(0.5*heigth/nely)/1]])
        self.Det=np.linalg.det(self.J)
        self.matnodes=matnodes
        self.matelems=matelems
        
    def D(self,E,nu):
        return E/((1+nu)*(1-2*nu))*np.array([[1-nu,nu,0],[nu,1-nu,0],[0,0,(1-2*nu)/2]])
    
##------------------------------- PROCESSING ------------------------------------##
## shape functions derivatives matrix -------------------------
    def N(self,x,y):
        n1=1/4 * (x-1)*(y-1)
        n2=-1/4  * (x+1)*(y-1)
        n3=1/4 * (x+1)*(y+1)    
        n4=-1/4 * (x-1)*(y+1)
        output=np.array([[n1,0,n2,0,n3,0,n4,0],[0,n1,0,n2,0,n3,0,n4]])
        return output
   ## derivatives matrix  :
    def B(self,x,y):
        output= 1/4 * \
                np.array([[y-1,0,1-y,0,y+1,0,-1-y,0], \
                          [0,x-1,0,-1-x,0,x+1,0,1-x], \
                     [x-1,y-1,-1-x,1-y,x+1,y+1,1-x,-1-y]])
        output[0,:]=output[0,:]/(self.el_width/2) 
        output[1,:]=output[1,:]/(self.el_heigth/2)
        output[2,::2]=output[2,::2]/(self.el_heigth/2)
        output[2,1::2]=output[2,1::2]/(self.el_width/2)
        return output

## element stiffness matrix ----------------------------------------------------
    def ke(self):
        ke=np.zeros((8,8,3))
        for gp in range(4):
            x=self.allgp[0,gp]
            y=self.allgp[1,gp]
            ke[:,:,0]=ke[:,:,0]+1*self.Det*(np.transpose(self.B(x,y))@(self.D(E[0],nu[0])@self.B(x,y)))
            ke[:,:,1]=ke[:,:,1]+1*self.Det*(np.transpose(self.B(x,y))@(self.D(E[1],nu[1])@self.B(x,y)))
            ke[:,:,2]=ke[:,:,2]+1*self.Det*(np.transpose(self.B(x,y))@(self.D(E[2],nu[2])@self.B(x,y)))
        return ke

    def me(self):
        mei=np.zeros((8,8,3))
        for gp in range(4):
            x=self.allgp[0,gp]
            y=self.allgp[1,gp]
            mei[:,:,0]=mei[:,:,0]+1*densty[0]*self.Det*(np.transpose(self.N(x,y))@self.N(x,y))
            mei[:,:,1]=mei[:,:,1]+1*densty[1]*self.Det*(np.transpose(self.N(x,y))@self.N(x,y))
            mei[:,:,2]=mei[:,:,2]+1*densty[2]*self.Det*(np.transpose(self.N(x,y))@self.N(x,y))
        return mei

    def kbase(self):
        Kbase=np.zeros(((self.nelx+1)*(self.nely+2)*2,(self.nelx+1)*(self.nely+2)*2))
        ke=self.ke() #put at the begining to call it only once
        etpl=self.etpl
        matelems=self.matelems
        for i in range(self.nels): #nels
            for row in range(4):
                for col in range(4):
                    matelem=int(matelems[i,0])
                    indr=etpl[i][row]*2; indc=etpl[i][col]*2
                    Kbase[indr:(indr+2),indc:(indc+2)]=Kbase[indr:(indr+2),indc:(indc+2)] + \
                                                   (ke[:,:,matelem])[row*2:(row*2+2),col*2:(col*2+2)]
        return Kbase

    def mbase(self):
        Mbase=np.zeros(((self.nelx+1)*(self.nely+2)*2,(self.nelx+1)*(self.nely+2)*2))
        me=self.me() #put at the begining to call it only once
        matelems=self.matelems
        etpl=self.etpl
        for i in range(self.nels): #nels
            matelem=int(matelems[i,0])
            for row in range(4):
                for col in range(4):
                    indr=etpl[i][row]*2; indc=etpl[i][col]*2
                    Mbase[indr:(indr+2),indc:(indc+2)]=Mbase[indr:(indr+2),indc:(indc+2)] + \
                                                   (me[:,:,matelem])[row*2:(row*2+2),col*2:(col*2+2)]
        return Mbase


### dynamic solver:
class output_dyn:
    def __init__(self,setup,K,Mbase_inv,BC_dofs,BC_val,P_horiz,P_vert,d1,d2,p,cond,F_int, \
                 Utold,Ut,dt,condrup,it,disp,vel,condition_simul,coord,D,nu,E,heigth, \
                width,nelx,nely,all_dofs):
        
        
 ### in case we need the theoretical deforation of a rectangular solid:
        exx=nu*np.average(P_vert)/E; # strain in horizontal direction
        eyy=np.average(P_vert)/E; # strain in vertical direction
        dxy=heigth*np.average(P_horiz)/(E/(2*(1+nu))) # shear displacement (horizontal direction)

        ## initialize force vector
        F=np.zeros(((nelx+1)*(nely+2)*2,1))

        ## assemble degrees of freedom
        bot_nodes=np.where(coord[:,1]==0)[0]
        dofs0=bot_nodes*2  # dofs x direction
        dofs1=bot_nodes*2+1  # dofs y direction 
        right_nodes=np.where(coord[:,0]==width)[0]
        dofs2=right_nodes*2  # dofs x direction
        dofs3=right_nodes*2+1  # dofs y direction         
        top_nodes=np.where(coord[:,1]==heigth)[0]
        dofs4=top_nodes*2  # dofs x direction
        dofs5=top_nodes*2+1  # dofs y direction 
        left_nodes=np.where(coord[:,0]==0)[0]
        dofs6=left_nodes*2  # dofs x direction
        dofs7=left_nodes*2+1  # dofs y direction 

        #extra_dofs:
        nodes_11=np.where(coord[:,0]==0)[0]
        nodes_12=np.where(coord[nodes_11,1]<=d1)[0]
        nodes_1=nodes_11[nodes_12]
        dofs8=nodes_1*2  # dofs x direction
        dofs9=nodes_1*2+1  # dofs y direction

        nodes_21=np.where(coord[:,0]==width)[0]
        nodes_22=np.where(coord[nodes_21,1]>=(heigth-d2))[0]
        nodes_2=nodes_21[nodes_22]
        dofs10=nodes_2*2  # dofs x direction
        dofs11=nodes_2*2+1  # dofs y direction
        
        dofs_BC=[dofs0,dofs1,dofs2,dofs3,dofs4,dofs5,dofs6,dofs7,dofs8,dofs9,dofs10,dofs11] ##all dofs on Boundaries as shown on drawing
        # +4 extra dofs groups
        
        nmn_dofs=np.concatenate(np.array(dofs_BC)[BC_dofs[1]])
        f=[];
        for expr in BC_val[1]:  
            f=np.concatenate((f,eval(expr)))
        F[nmn_dofs,0]=f
        
#         #==============add interface forces===============#
        dofs_int=setup.dofs_int
#         #=================or add springs==================#
        if it==0: #add the srpings only at the begining
            for i,j in zip(np.append(dofs_int[0,:],dofs_int[2,:]), \
                           np.append(dofs_int[1,:],dofs_int[3,:])):
                K[i,j]=K[i,j]-p
                K[j,i]=K[j,i]-p
                K[i,i]=K[i,i]+p
                K[j,j]=K[j,j]+p   
            
        dofs_condrup=dofs_int[:,condrup[0]] 
        for i,j in zip(dofs_condrup[0,:],dofs_condrup[1,:]):
            K[i,j]=K[i,j]+p
            K[j,i]=K[j,i]+p
            K[i,i]=K[i,i]-p
            K[j,j]=K[j,j]-p
        
#         Ksolve[indnonzero,indnonzero]=K[indnonzero,indnonzero]+Kspring[indnonzero,indnonzero]
    ## it won't work because we add only the diagonal terms this way...
        fint=np.transpose(np.array([F_int]))
        F[dofs_int[0,:]]=+fint #bottom dofs, horiz
        F[dofs_int[1,:]]=-fint #top dofs, horiz
           #=================Displacement BCs==================# 
        d=[];
        for expr in BC_val[0]:
            d=np.concatenate((d,eval(expr)))
   
        fixed_dofs=np.concatenate(np.array(dofs_BC)[BC_dofs[0]])
        free_dofs=all_dofs[np.in1d(all_dofs, fixed_dofs,invert=True)]
        print(fixed_dofs.shape)
        #===============solve static or dynamic system==================#
        if condition_simul=='static':
            for id in range(len(fixed_dofs)):
                i=fixed_dofs[id];
                rows=K[free_dofs,i];
                F[free_dofs,0]=F[free_dofs,0]-rows*d[id];
                K[:,i]=0; K[i,:]=0; K[i,i]=1;
                F[i]=d[id];
            
            ## solve the system with a sparse matrix ------------------------------------------
            Ksp = csr_matrix(K)
            U= spsolve(Ksp,F)
            self.U=U
#         ## solve the system ------------------------------------------
        else:
            Ut_acc=deepcopy(Ut)
            for id in range(len(fixed_dofs)):
                i=fixed_dofs[id];
                Ut_acc[i]=d[id]
            F[fixed_dofs,0]=(K@Ut_acc)[fixed_dofs] ##tranforms imposed displacement into imposed force

#####             F[fixed_dofs,0]=((K[:,fixed_dofs][fixed_dofs,:])@(Ut[fixed_dofs]))#[fixed_dofs]... slower version
            acc=Mbase_inv@(F[:,0]-K@Ut_acc) ## uses K@Utold(modified to ut) to get acceleration when comparing old displacements
#            vel[free_dofs]=vel[free_dofs]+dt*acc[free_dofs]
#            Ut[free_dofs]=Ut[free_dofs]+dt*vel[free_dofs]
            Ut[free_dofs]=dt**2*acc[free_dofs]+2*Ut[free_dofs]-Utold[free_dofs]
            self.Ut=Ut
            self.vel=(Ut-Utold)/(2*dt)
            
def stressplotter(nelx,nely,Utnew,D,E,nu,el_heigth,el_width,matnodes):
    sigma_b=np.zeros((int(nely/2)+1,nelx+1,3))
    sigma_t=np.zeros((int(nely/2)+1,nelx+1,3))
    ux_mat=Utnew[::2].reshape((nely+2, nelx+1))
    uy_mat=Utnew[1::2].reshape((nely+2, nelx+1))
    ux_mat_b=ux_mat[0:int(nely/2)+1,:]
    uy_mat_b=uy_mat[0:int(nely/2)+1,:]
    ux_mat_t=ux_mat[int(nely/2)+1:,:]
    uy_mat_t=uy_mat[int(nely/2)+1:,:]

    ##botom
    Dub=np.gradient(ux_mat_b,setup1.el_heigth,setup1.el_width); #to check width or heigth first
    Dvb=np.gradient(uy_mat_b,setup1.el_heigth,setup1.el_width);
    Dxub=Dub[1]
    Dyub=Dub[0]
    Dxvb=Dvb[1]
    Dyvb=Dvb[0]
    Bstress=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,1]]);
    n_nd=0
    for i in range(int(nely/2)+1):
        for j in range(nelx+1):
            strain=np.matmul(Bstress,np.array([[Dxub[i,j]],[Dyvb[i,j]],[Dxvb[i,j]],[Dyub[i,j]]]));
            matnode=int(matnodes[n_nd,0])        
            stress=np.matmul(D(E[matnode],nu[matnode]),strain); 
            sigma_b[i,j,0]=stress[0]; # 0 for sxx
            sigma_b[i,j,1]=stress[1]; # 1 for syy
            sigma_b[i,j,2]=stress[2]; # 2 for tauxy
            n_nd=n_nd+1
    #top
    Dut=np.gradient(ux_mat_t,setup1.el_heigth,setup1.el_width); #to check width or heigth first
    Dvt=np.gradient(uy_mat_t,setup1.el_heigth,setup1.el_width);
    Dxut=Dut[1]
    Dyut=Dut[0]
    Dxvt=Dvt[1]
    Dyvt=Dvt[0]
    Bstress=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,1]]);
    for i in range(int(nely/2)+1):
        for j in range(nelx+1):
            strain=np.matmul(Bstress,np.array([[Dxut[i,j]],[Dyvt[i,j]],[Dxvt[i,j]],[Dyut[i,j]]]));
            matnode=int(matnodes[n_nd,0])        
            stress=np.matmul(D(E[matnode],nu[matnode]),strain); 
            sigma_t[i,j,0]=stress[0]; # 0 for sxx
            sigma_t[i,j,1]=stress[1]; # 1 for syy
            sigma_t[i,j,2]=stress[2]; # 2 for tauxy
            n_nd=n_nd+1
    return [sigma_b,sigma_t]

## define friction law:
## slip weakening friction (simplified)
def friction_sl(d0,condfric,tan_disp_i,fint,p,el_width, p_init,p_final,Dc):
    ##include the sign of the velocity later
    f_init=p_init*el_width
    f_final=p_final*el_width
    if len(fint[condfric[0]]>0):
        fint[condfric[0]] = -(f_init-(f_init-f_final)*np.abs((tan_disp_i[0][condfric[0]]-d0[condfric[0]])/Dc))
        fint[np.where(np.logical_and(fint>=-f_final,fint!=0))]=-f_final ##careful with the sign and >= or <=
    return fint

def friction_sl_2(TDold,TDi,VD,p,Dc,mu_s,mu_d,cum_disp_old,condfric,signf):
    ## cumulative displacement
    cum_disp[0,condfric[0]]=cum_disp_old[0,condfric[0]] + (np.abs(TDi-TDold))[0,condfric[0]]
    ## friction coefficient
    mu=mu_s-((mu_s-mu_d)/Dc)*cum_disp_old; # check if have to use old or new: more stable?
    mu[np.where(mu<mu_d)]=mu_d
    ## final tractions
    fint[condfric[0]]=-signf[0,condfric[0]]*mu[0,condfric[0]]*VD[0,condfric[0]]*p ## REMOVE *0
    return [fint,cum_disp]


## Boundary conditions ----------------------------------------------------------      

#      --4/5--  10/11  even numbers -> horizontal dofs
#      |      |   odd numbers -> vert dofs
#      6/7   2/3
#      |      |
# 8/9  --0/1-- 

## =============================static solver====================================:
import scipy.sparse as sparse
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
import time

#E=[200e9,50.0e9,0.1e9] ##stainless steel, granite, silicon/rubber
#nu=[0.3,0.25,0.4]
#densty=[7700,2700,1500]
# ## homogeneous
E=[50.0e9,50.0e9,50.0e9] ##stainless steel, granite, silicon/rubber
nu=[0.25,0.25,0.25]
densty=[2700,2700,2700]
##========================geometry of (initial) static model======================##
setup1=setup(70,20,0.07,0.02,E,nu,densty)
##=============================Boundary conditions================================##
P_horiz=1.48e6
P_vert=-3.8e6
disp=0
BC_dofs=[[6,1,2,5],[4]] #we always give the displacements first, then the forces
## here all the BC values must be expression giving vectors
## to give pressure we must normalize in those expressions 
## by th length of elements (height or width)
BC_val=[['0.0*np.ones((len(dofs6)))', \
         '0.0*np.ones((len(dofs1)))', \
         '0.0*np.ones((len(dofs2)))- 25e-6+disp' ,\
        '0.0*np.ones((len(dofs5)))- 8e-6'], \
        ['0.0*np.ones((len(dofs4)))']] ## a dummy force is added for the code to keep running
F_int='0*1e6*np.ones((len(dofs_int[0,:])))'
# ##oading o homogeneous materials
BC_dofs=[[8,1,10,5],[4]] #we always give the displacements first, then the forces
BC_val=[['(0.0*np.ones((len(dofs8))))', \
         '0.0*np.ones((len(dofs1)))', \
         '(0.0*np.ones((len(dofs10)))- 40e-6+disp)' ,\
        '0.0*np.ones((len(dofs5)))- 8e-6'], \
        ['0.0*np.ones((len(dofs4)))']] ## a dummy force is added for the code to keep running
E=np.array(setup1.E)
nu=np.array(setup1.nu)
densty=np.array(setup1.densty)
nelx=setup1.nelx
nely=setup1.nely
heigth=setup1.heigth
coord=setup1.coord
p=setup1.p
el_width=setup1.el_width
densty=np.array(setup1.densty);
G=E/(2*(1+nu));
vs=(G/densty)**0.5;
vp=((2*G*(1-nu))/(densty*(1-2*nu)))**0.5;
Kspring=np.zeros(setup1.kbase().shape)
all_dofs=setup1.all_dofs
D=setup1.D
coord=setup1.coord
width=setup1.width
matnodes=setup1.matnodes
## --vector initialisation
nb_int_nodes=len(setup1.dofs_int[0])
fint=np.zeros(nb_int_nodes)
Utold=np.zeros(nb_int_nodes)
Ut=np.zeros(nb_int_nodes)
force=np.zeros(nb_int_nodes)
dt=0;it=0
condrup=np.where(fint==1)
disp=0
vel=0
int_nodes=np.where(np.logical_and(coord[:,1]==(heigth/2),np.transpose(matnodes[:,0])==1))[0]
dofs_int=np.append(int_nodes*2,int_nodes*2+1).reshape(4,(nb_int_nodes))



output1=output_dyn(setup1,setup1.kbase(),setup1.mbase(),BC_dofs,BC_val,P_horiz,P_vert,0.009,0.009,p, \
                       'penalty',fint,Utold,Ut,dt,condrup,0,disp,vel,'static',coord,D,nu,E,heigth, \
                width,nelx,nely,all_dofs)


U=output1.U
coord=setup1.coord
f42=plt.figure(42)
ax42=f42.add_subplot(111)
tan=np.array([U[dofs_int[1,:]]-U[dofs_int[0,:]]])[0,:]*p/el_width
ax42.plot(setup1.coordint,-tan,'b')
normal=np.array([U[dofs_int[3,:]]-U[dofs_int[2,:]]])[0,:]*p/el_width;
corr_norm=np.ones(normal.shape);corr_norm[0]=2;corr_norm[-1]=2
ax42.plot(setup1.coordint,-normal*corr_norm,'--b')
ax42.set_ylabel('Stresses [MPa]')
ax42.set_xlabel('x [m]')
ax43=ax42.twinx()
ax43.plot(setup1.coordint,tan/(normal*corr_norm),'g')
ax43.set_ylabel('$\sigma_{xy}$'+'/'+'$\sigma_{yy}$')
ax42.legend(['$\sigma_{xy}$','$\sigma_{yy}$'])
ax42.yaxis.label.set_color('blue')
ax42.tick_params(axis='y', colors='blue')
ax43.yaxis.label.set_color('green')
ax43.tick_params(axis='y', colors='green')
f42.tight_layout()


from mpl_toolkits.axes_grid1 import make_axes_locatable
f61=plt.figure(61)
ax61=f61.add_subplot(111)
Utnew=output1.U
[sigma_b,sigma_t]=stressplotter(nelx,nely,Utnew,D,E,nu,setup1.el_heigth,setup1.el_width,setup1.matnodes)
factor=1
st2plt=2
stressmap1=ax61.pcolormesh(((coord[:,0]+Utnew[::2]*factor).reshape((nely+2, nelx+1)))[0:int(nely/2)+1,:], \
                      ((coord[:,1]+Utnew[1::2]*factor).reshape((nely+2, nelx+1)))[0:int(nely/2)+1,:], \
                      -sigma_b[:,:,st2plt]/1e6, shading='Gouraud')
stressmap2=ax61.pcolormesh(((coord[:,0]+Utnew[::2]*factor).reshape((nely+2, nelx+1)))[int(nely/2)+1:,:], \
                      ((coord[:,1]+Utnew[1::2]*factor).reshape((nely+2, nelx+1)))[int(nely/2)+1:,:], \
                      -sigma_t[:,:,st2plt]/1e6, shading='Gouraud')

stressmap1.set_clim(np.min(sigma_t[:,:,st2plt])/1e6,np.max(sigma_t[:,:,st2plt])/1e6)
stressmap2.set_clim(np.min(sigma_t[:,:,st2plt])/1e6,np.max(sigma_t[:,:,st2plt])/1e6)
# stressmap1.set_clim(0,60)
# stressmap2.set_clim(0,60)
divider = make_axes_locatable(ax61)
cax1 = divider.append_axes("right", size="5%", pad=0.15)

ax61.set_xlabel('x [m]')
ax61.set_ylabel('y [m]')
ax61.set_aspect('equal');
# Visualizing colorbar part -start
plt.colorbar(stressmap2, cax=cax1,label='$\sigma_{xy}$ [MPa]')


from tqdm import tqdm
# =========================simulation parameters and initialisation ==============##
## (we call all the fixed variables before the loop to make the code faster)
p_init=60e6
p_final=40e6
Dc=10000e-6
mu_s_lim=np.max(tan/(normal*corr_norm))*1.000000001
mu_d=mu_s_lim/1.2
## --mechanical and geometrical parameters--:
condition_simul='dynamic'
Kdyn=setup1.kbase()
Mdyn=np.linalg.inv(setup1.mbase())

## -- time steps parapmeters --:

dt=0.0005*setup1.el_width/np.max(vp); #timestep according to courant stability criterion
#tmax=0.01e-3
tmax=80*dt
t=np.linspace(0,tmax,np.round(tmax/dt));nt=len(t);

# total_disp=-0.001e-6
disp_rate=0.03e-6; # 1 micro meter/s
total_force=P_horiz*2*0.55
loading_rate=total_force/tmax # for uniform displacement rate
# disp_all=np.linspace(0,total_disp,nt)
disp_all=t*disp_rate
##--use the static displacement as initial conditions--:
Utold=output1.U
Ut=output1.U
## explicit time-stepping scheme==================#

vel=deepcopy(Ut*0)

disp_stored=np.array([Ut])
tan_disp=np.array([Ut[dofs_int[1,:]]-Ut[dofs_int[0,:]]]) ## they should be 0
condrup=np.where((tan_disp[0]*p/el_width)<=-60e6)
condfric=condrup
fint=np.zeros(nb_int_nodes)
d0=np.zeros(nb_int_nodes)
cum_disp=np.zeros((1,nb_int_nodes))
mu_s=np.zeros((1,nb_int_nodes))
signf=np.zeros((1,nb_int_nodes))
signf[0,:]=np.sign(tan);
##=======================print the simulation parameters============================##
print('loading rate = ' + str(loading_rate) + '[m/s]')
print('total time = ' + str(tmax) + '[s]')
print('nb of iterations = ' + str(nt))
print('time step = ' + str(dt) + '[s]')
print('distance reached by S-waves = ' + str(tmax*vs*100) + '[cm]')
print('time to reach interface = ' + str(0.075/vs) + '[s]')
## ======================= run explicit time-stepping scheme =======================##
if 'pbar' in locals():
    del pbar
    print('deleted')
plt.close('all')
plt.ion()
pbar = tqdm(total=nt)

for it in range(nt):
    
    disp=disp_all[it]
    TDold=np.array([(Ut[dofs_int[1,:]]-Ut[dofs_int[0,:]])]) 
    
    output3=output_dyn(setup1,Kdyn,Mdyn,BC_dofs,BC_val,P_horiz,P_vert,0,0,p, \
                       'penalty',fint,Utold,Ut,dt,condrup,it,disp,vel,condition_simul,coord,D,nu,E,heigth, \
                width,nelx,nely,all_dofs)

    vel=output3.vel
    Utold=Ut;
    Ut=output3.Ut

    TDi=np.array([(Ut[dofs_int[1,:]]-Ut[dofs_int[0,:]])])  
    VD=np.array([(Ut[dofs_int[3,:]]-Ut[dofs_int[2,:]])])  
#     condrup=np.where(np.logical_and((tan_disp_i[0]*p/el_width)<=-p_init,fint==0))
#     condfric=np.where((tan_disp_i[0]*p/el_width)<=-p_init)
    condrup=np.where(np.logical_and(np.abs(TDi[0]/VD[0])>mu_s_lim,fint==0))
    condfric=np.where(np.abs(TDi[0]/VD[0])>mu_s_lim)
    mu_s[:,condrup[0]]=np.abs(TDi/VD)[:,condrup[0]]
#    
#    d0[condrup]= TDi[0][condrup]
#    fint= friction_sl(d0,condfric,TDi,fint,p,el_width, p_init,p_final,Dc)
    ##alternative friction law:
    [fint,cum_disp]=friction_sl_2(TDold,TDi,VD,p,Dc,mu_s,mu_d,cum_disp,condfric,signf)
    signf=np.sign(TDi-TDold); signf[np.where(VD>0)]=0 

##=========uncomment those lines to have rupture:
    pbar.update(1)
    ## update stiffness matrix with spring stiffnesses
    if it % 1 ==0:
#         print(mu_s)
#         print(t[it])
        disp_stored=np.append(disp_stored,np.array([Ut]),axis=0)
        tan_disp=np.append(tan_disp,TDi,axis=0)
        
pbar.close()


f6=plt.figure(6)
ax6=f6.add_subplot(111)

ax6.plot(setup1.coordint,np.transpose(tan_disp)*setup1.p/setup1.el_width);
ax6.set_ylim([-60e6,0])
#ax6.axhline
#f7=plt.figure(7)
#ax7=f7.add_subplot(111)
#f7.show()
#f7.canvas.draw()
#plot_item='displacements '
#
#a=0
#for it in range(nt):
#    if it % 2 ==0:
#        ax7.clear()
#        Utnew=disp_stored[a]
#        if plot_item=='displacements':
#            im=ax7.plot(coord[:,0]+Utnew[::2]*100,coord[:,1]+Utnew[1::2]*1,'xr') #U[::2,0] if not using sparse
#            f7.canvas.draw()   # draw
#            time.sleep(0.1); 
#        else:
#
#            
#            [sigma_b,sigma_t]=stressplotter(nelx,nely,Utnew,D,E,nu,setup1.el_heigth,setup1.el_width,setup1.matnodes)
#            factor=30
#            stressmap1=ax7.pcolormesh(((coord[:,0]+Utnew[::2]*factor).reshape((nely+2, nelx+1)))[0:int(nely/2)+1,:], \
#                                  ((coord[:,1]+Utnew[1::2]*factor).reshape((nely+2, nelx+1)))[0:int(nely/2)+1,:], \
#                                  -sigma_b[:,:,2]*1e-6, shading='Gouraud')
#
#            stressmap2=ax7.pcolormesh(((coord[:,0]+Utnew[::2]*factor).reshape((nely+2, nelx+1)))[int(nely/2)+1:,:], \
#                                  ((coord[:,1]+Utnew[1::2]*factor).reshape((nely+2, nelx+1)))[int(nely/2)+1:,:], \
#                                  -sigma_t[:,:,2]*1e-6, shading='Gouraud')
#
##             stressmap1.set_clim(np.min(sigma_t[:,:,2]),np.max(sigma_t[:,:,2]))
##             stressmap2.set_clim(np.min(sigma_t[:,:,2]),np.max(sigma_t[:,:,2]))
#            stressmap1.set_clim(0,60)
#            stressmap2.set_clim(0,60)  
##             stressmap1.set_clim(0,0.5)
##             stressmap2.set_clim(0,0.5) 
#            ax7.set_aspect('equal')
#            cbar=f7.colorbar(stressmap1,ax=ax7)
#            f7.canvas.draw()   # draw
#            time.sleep(0.1); 
#            cbar.remove()
#            
#        a=a+1
##         f7.savefig("DGdynforce\\" + 'a' + '{0:05}'.format(it) + '.png')
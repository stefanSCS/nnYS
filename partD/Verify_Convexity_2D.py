
import autograd.numpy as anp
from autograd import elementwise_grad as egrad 
import time 


def genPointsSphere(ndim=3,nRandom=10**5):
    anp.random.seed(99)     
    vuRandom=anp.random.normal(0,1,(nRandom,ndim)) ##nPoints on (n-1)D unit sphere of (ndim)D
    vNorm=anp.sqrt(anp.sum(vuRandom**2,axis=1)).reshape((nRandom,1))
    vuRandom[:,:]/=vNorm
    return vuRandom



def func2D_Mises(x1,x2,x3,param):
    return x1*x1-x1*x2+x2*x2+3.0*x3*x3        



def func2D_Poly6(x1,x2,x3,param):
    yy=(param[0]*x1**6+param[1]*x1**5*x2+param[2]*x1**4*x2**2+param[3]*x1**3*x2**3+
        param[4]*x1**2*x2**4+param[5]*x1*x2**5+param[6]*x2**6)
    yy=yy+(param[7]*x1**4+param[8]*x1**3*x2+param[9]*(x1*x2)**2+param[10]*x1*x2**3+param[11]*x2**4)*(x3**2)
    yy=yy+(param[12]*x1**2+param[13]*x1*x2+param[14]*x2**2)*(x3**4)+param[15]*x3**6
    return yy
    

def func2D_Poly8(x1,x2,x3,param):
    yy=(param[0]*x1**8+param[1]*x1**7*x2+param[2]*x1**6*x2**2+param[3]*x1**5*x2**3+param[4]*(x1*x2)**4+
        param[5]*x1**3*x2**5+param[6]*x1**2*x2**6+param[7]*x1*x2**7+param[8]*x2**8)
    yy=yy+(param[9]*x1**6+param[10]*x1**5*x2+param[11]*x1**4*x2**2+param[12]*(x1*x2)**3+
           param[13]*x1**2*x2**4+param[14]*x1*x2**5+param[15]*x2**6)*(x3**2)
    yy=yy+(param[16]*x1**4+param[17]*x1**3*x2+param[18]*(x1*x2)**2+param[19]*x1*x2**3+param[20]*x2**4)*(x3**4)
    yy=yy+(param[21]*x1**2+param[22]*x1*x2+param[23]*x2**2)*(x3**6)+param[24]*x3**8
    return yy

def nnYS_readData(fName):
    nDim=3 ##input dimension
    vw=[]
    ff=open(fName)
    for line in ff:
        line=line.strip()
        if(line):vw.append(float(line))
    ff.close()
    param=[]
    nLayers=int(vw[0])
    kStart,kLayer,nHom=2*nLayers-1,1,1.0
    for kk in range(nLayers-1):
        nCol=int(vw[kLayer])
        W=anp.zeros((nDim,nCol))
        kLayer+=1; nHom*=vw[kLayer]; kLayer+=1
        for jj in range(nCol):
            kEnd=kStart+nDim
            W[:,jj]=vw[kStart:kEnd]
            kStart=kEnd
        nDim=nCol    
        param.append(W)
    W=anp.zeros((nDim,1))
    W[:,0]=vw[kStart:]
    param.append(W)    
    return nHom,param        
        

def nnYS_nb10(x1,x2,x3,W):
    vvx=anp.column_stack((x1,x2,x3))
    va=(anp.dot(vvx,W[0]))**2
    va=(anp.dot(va,W[1]))**2
    va=(anp.dot(va,W[2]))**2
    va=(anp.dot(va,W[3]))**2
    return anp.dot(va,W[4])

paramMises=None


paramPoly6_AA2090T3_IJP2008=[1.0, -1.1059, 2.5255, -5.1914, 6.1458, -4.3254, 1.7753, 14.190, 
            -4.9759, -4.3926, 3.4652, 15.806, 0.0000, -9.4916, 86.661, 116.42]

paramPoly8_AA2090T3_IJP2008=[
    1.0, -1.3376, 2.1967, -5.7867, 12.312, -16.000, 13.260, -7.0415, 2.1508, 
    12.697, -4.2719, 74.294, -31.487, -18.387, -11.096, 28.553, 116.65, 
    -238.83, 158.71, -57.545, 172.28, -0.2559, -8.2036, 558.15, 543.50]

paramPoly8_AA6016T4_CM2023=[
    1.0,-2.7575360419397117,11.530135957809993,-25.895206461605195,34.38295306726676,
    -27.913685793880834,13.264500537932427,-3.4446490737696935,1.1470222893313808,
    19.694927566274288,-22.728141338803983,97.14527508865311,-118.15219232701236,
    79.56140212110728,-22.863865241214093,26.876869442462294,72.53986457067393,
    -76.33248924405487,263.73730517826544,-61.653927738263796,77.5711377175994,
    144.52991577535954,-151.25820450588395,144.06559886753368,77.06411729917284]     


if(0):
    nHom=2
    param=paramMises
    func2D = func2D_Mises
if(0):    
    nHom=6
    param=paramPoly6_AA2090T3_IJP2008
    func2D = func2D_Poly6    
if(0):    
    nHom=8
    param=paramPoly8_AA2090T3_IJP2008
    func2D = func2D_Poly8
if(1):    
    nHom=8
    param=paramPoly8_AA6016T4_CM2023
    func2D = func2D_Poly8    
if(0):
    nHom,param=nnYS_readData('nb10_NNdata.txt')
    func2D = nnYS_nb10    


def func(x1,x2,x3,param=param):
    return func2D(x1,x2,x3,param)    



df1=egrad(func,0)
df2=egrad(func,1)
df3=egrad(func,2)
h11=egrad(df1,0)
h12=egrad(df1,1)
h13=egrad(df1,2)
h22=egrad(df2,1)
h23=egrad(df2,2)
h33=egrad(df3,2)


tA=time.monotonic()


nRandom=10**5
print('Checking convexity at {x} points...'.format(x=nRandom))

vx=genPointsSphere(ndim=3,nRandom=nRandom)


RnHom=1.0/nHom
vtau=(1.0/func(vx[:,0],vx[:,1],vx[:,2],param)**RnHom).reshape((nRandom,1))
vx[:,0]*=vtau[:,0];vx[:,1]*=vtau[:,0];vx[:,2]*=vtau[:,0]

vdf1=df1(vx[:,0],vx[:,1],vx[:,2],param).reshape((nRandom,1))
vdf2=df2(vx[:,0],vx[:,1],vx[:,2],param).reshape((nRandom,1))
vdf3=df3(vx[:,0],vx[:,1],vx[:,2],param).reshape((nRandom,1))
vh11=h11(vx[:,0],vx[:,1],vx[:,2],param).reshape((nRandom,1))
vh12=h12(vx[:,0],vx[:,1],vx[:,2],param).reshape((nRandom,1))
vh13=h13(vx[:,0],vx[:,1],vx[:,2],param).reshape((nRandom,1))
vh22=h22(vx[:,0],vx[:,1],vx[:,2],param).reshape((nRandom,1))
vh23=h23(vx[:,0],vx[:,1],vx[:,2],param).reshape((nRandom,1))
vh33=h33(vx[:,0],vx[:,1],vx[:,2],param).reshape((nRandom,1))


normGF=(vdf1[:,0]*vdf1[:,0]+vdf2[:,0]*vdf2[:,0]+vdf3[:,0]*vdf3[:,0]).reshape((nRandom,1))
trueNorm=anp.sqrt(normGF).reshape((nRandom,1))
vdf1[:,0]/=trueNorm[:,0];vdf2[:,0]/=trueNorm[:,0];vdf3[:,0]/=trueNorm[:,0]


HS11 = (vh22[:,0]*vh33[:,0] - vh23[:,0]*vh23[:,0]).reshape((nRandom,1))
HS12 = (-(vh12[:,0]*vh33[:,0] - vh13[:,0]*vh23[:,0])).reshape((nRandom,1))
HS13 = (vh12[:,0]*vh23[:,0] - vh13[:,0]*vh22[:,0]).reshape((nRandom,1))
HS22 = (vh11[:,0]*vh33[:,0] - vh13[:,0]*vh13[:,0]).reshape((nRandom,1))
HS23 = (-(vh11[:,0]*vh23[:,0] - vh12[:,0]*vh13[:,0])).reshape((nRandom,1))
HS33 = (vh11[:,0]*vh22[:,0] - vh12[:,0]*vh12[:,0]).reshape((nRandom,1))
HS = anp.zeros((nRandom,1))
HS[:,0] = HS11[:,0]*vdf1[:,0]*vdf1[:,0]+HS22[:,0]*vdf2[:,0]*vdf2[:,0]+HS33[:,0]*vdf3[:,0]*vdf3[:,0]
HS[:,0]+= 2.0*(HS12[:,0]*vdf1[:,0]*vdf2[:,0]+HS13[:,0]*vdf1[:,0]*vdf3[:,0]+HS23[:,0]*vdf2[:,0]*vdf3[:,0])
KG = HS[:,0]/normGF[:,0]
kmin,kmax=anp.argmin(KG),anp.argmax(KG)
KGmin,KGmax=KG[kmin],KG[kmax]


vdet1 = vh22[:,0]*vh33[:,0] - vh23[:,0]*vh23[:,0]
vdet2 = vh11[:,0]*vh33[:,0] - vh13[:,0]*vh13[:,0]
vdet3 = vh11[:,0]*vh22[:,0] - vh12[:,0]*vh12[:,0]
vdet = (vh11[:,0]*vdet1[:] - vh12[:,0]*(vh12[:,0]*vh33[:,0]-vh13[:,0]*vh23[:,0]) +
                             vh13[:,0]*(vh12[:,0]*vh23[:,0]-vh13[:,0]*vh22[:,0]))

lpd1, lpd2, lpd3 = anp.min(vh11), anp.min(vh22), anp.min(vh33)
lpd4, lpd5, lpd6 = anp.min(vdet1), anp.min(vdet2), anp.min(vdet3) 
lpd7 = anp.min(vdet)   

tB=time.monotonic()

print('execution time = {x} (s)'.format(x=tB-tA))
print('Principal minors:')
print(lpd1, lpd2, lpd3, lpd4, lpd5, lpd6, lpd7)
print('KGmin = ', KGmin)
print(vx[kmin,:])
print('KGmax = ', KGmax)
print(vx[kmax,:])
    
    
    
    
    
    




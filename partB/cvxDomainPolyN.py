


import numpy as np
from matplotlib import pyplot as plt



                          

def vPolyMons(degree):
    dd={'nQ':degree,'nMons':int((int(degree/2)+1)**2)}
    vQ=[];vQidx=[];jj=0
    while(jj<=degree):
        lv=len(vQ)
        vQidx.append([jj,[k for k in range(lv,lv+degree+1-jj)]])
        vQ+=[(degree-jj-k,k,jj) for k in range(degree+1-jj)]
        jj+=2 
    dd['vQ']=vQ;dd['vQidx']=vQidx
    dd['vD1Q']=[(mon[0]-1,mon[1],mon[2]) if(mon[0]>0) else (0,0,0) for mon in vQ]
    dd['vC1Q']=np.array([mon[0] if(mon[0]>0) else 0 for mon in vQ])
    dd['vD2Q']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in vQ]
    dd['vC2Q']=np.array([mon[1] if(mon[1]>0) else 0 for mon in vQ])
    dd['vD3Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in vQ]
    dd['vC3Q']=np.array([mon[2] if(mon[2]>0) else 0 for mon in vQ])
    dd['vH11Q']=[(mon[0]-1,mon[1],mon[2]) if(mon[0]>0) else (0,0,0) for mon in dd['vD1Q']]
    dd['vCH11Q']=np.array([cf*mon[0] for (cf,mon) in zip(dd['vC1Q'],dd['vD1Q'])])
    dd['vH12Q']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in dd['vD1Q']]
    dd['vCH12Q']=np.array([cf*mon[1] for (cf,mon) in zip(dd['vC1Q'],dd['vD1Q'])])
    dd['vH13Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD1Q']]
    dd['vCH13Q']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC1Q'],dd['vD1Q'])])
    dd['vH22Q']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in dd['vD2Q']]
    dd['vCH22Q']=np.array([cf*mon[1] for (cf,mon) in zip(dd['vC2Q'],dd['vD2Q'])])
    dd['vH23Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD2Q']]
    dd['vCH23Q']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC2Q'],dd['vD2Q'])])
    dd['vH33Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD3Q']]
    dd['vCH33Q']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC3Q'],dd['vD3Q'])])     
    return dd
    


def PolyNparam(deg,vvvCoeff,matProp=[],matName='',fPrint=False):
    nMon=int(deg/2)
    NN=nMon+1
    nCoeff=NN*NN
    vcf=np.zeros(2+12*nCoeff)
    vcf[0]=deg
    mxV=max([abs(xx) for xx in vvvCoeff])
    if(mxV<zroTol):print('All parameters are = 0\n Computations aborted');exit()
    if(False):
        xvcf[1]=mxV
        xvCoeff=[xx/mxV for xx in vvvCoeff]
    vcf[1]=nCoeff
    vCoeff=[xx for xx in vvvCoeff]    
    d1,d2=np.zeros(nCoeff),np.zeros(nCoeff)
    d11,d22,d12=np.zeros(nCoeff),np.zeros(nCoeff),np.zeros(nCoeff)
    kcf=0
    for powerXY in range(nMon):
        kdeg=deg-2*powerXY
        for jj in range(kdeg+1):
            ##calculate derivatives 
            d1[kcf]=(kdeg-jj)*vCoeff[kcf]
            d2[kcf]=jj*vCoeff[kcf]
            d11[kcf]=(kdeg-jj-1)*d1[kcf]
            d22[kcf]=(jj-1)*d2[kcf]
            d12[kcf]=jj*d1[kcf]
            kcf+=1
        ##shift backwards d/dSigmaY    
        d2[kcf-kdeg-1:kcf-1]=d2[kcf-kdeg:kcf]
        d2[kcf-1]=zro         
        ##shift backwards d2/(dSigmaY,dSigmaY)    
        d22[kcf-kdeg-1:kcf-2]=d22[kcf-kdeg+1:kcf]
        d22[kcf-1],d22[kcf-2]=zro,zro
        ##shift backwards d2/(dSigmaX,dSigmaY)    
        d12[kcf-kdeg-1:kcf-2]=d12[kcf-kdeg:kcf-1]
        d12[kcf-1],d12[kcf-2]=zro,zro
    if(1): ##make pure zero some awkward -0.0        
        for kcf in range(nCoeff):
            if(abs(d1[kcf])<zroTol): d1[kcf]=zro
            if(abs(d2[kcf])<zroTol): d2[kcf]=zro
            if(abs(d11[kcf])<zroTol): d11[kcf]=zro
            if(abs(d22[kcf])<zroTol): d22[kcf]=zro
            if(abs(d12[kcf])<zroTol): d12[kcf]=zro
    ###setup the reversed sets of coefficients 
    ###use indexing starting from 1 
    ###(to avoid the problem of vec[-1] which Python interprets as last position in vec)
    nOne=nCoeff+1
    zvcf=np.zeros(nOne)
    zd1,zd2=np.zeros(nOne),np.zeros(nOne)
    zd11,zd22,zd12=np.zeros(nOne),np.zeros(nOne),np.zeros(nOne)
    zvcf[1:],zd1[1:],zd2[1:],zd11[1:],zd22[1:],zd12[1:]=vCoeff,d1,d2,d11,d22,d12
    pvcf=np.zeros(nCoeff)
    pd1,pd2=np.zeros(nCoeff),np.zeros(nCoeff)
    pd11,pd22,pd12=np.zeros(nCoeff),np.zeros(nCoeff),np.zeros(nCoeff)
    kcf=1 ##index starts from 1
    for ii in range(0,nMon):
        kdeg=deg-2*ii
        kcf2=kcf+kdeg+1
        pvcf[kcf-1:kcf2-1:1]=zvcf[kcf2-1:kcf-1:-1] ##if kcf were 0 =>>Problem
        pd1[kcf-1:kcf2-2:1]=zd1[kcf2-2:kcf-1:-1]
        pd2[kcf-1:kcf2-2:1]=zd2[kcf2-2:kcf-1:-1]
        pd11[kcf-1:kcf2-3:1]=zd11[kcf2-3:kcf-1:-1]
        pd22[kcf-1:kcf2-3:1]=zd22[kcf2-3:kcf-1:-1]
        pd12[kcf-1:kcf2-3:1]=zd12[kcf2-3:kcf-1:-1]
        kcf=kcf2
    pvcf[-1]=vCoeff[-1]        
    vcf[2:]=np.concatenate((vCoeff,d1,d2,d11,d22,d12,pvcf,pd1,pd2,pd11,pd22,pd12),axis=0)    
    #return d1,d2,d11,d22,d12,pvcf,pd1,pd2,pd11,pd22,pd12
    if(fPrint):
        fName=matName+'_deg'+str(deg)+'_FEdata.txt'
        ff=open(figDirData+fName,'w')
        for item in matProp:
            ff.write('{:.12f}\n'.format(item))
        for item in vcf:
            ff.write('{:.12f}\n'.format(item))
        ff.close()        
    return vcf


def fYF(sx,sy,sxy,vcf):
    deg=int(vcf[0])
    nMon=int(deg/2)
    nCoeff=int(vcf[1])
    KKo=nCoeff
    KKo-=1
    msx,msy=abs(sx),abs(sy)
    if(msx+msy<gTol):
        #KKo+=2
        #tau=abs(sxy)
        #yf=vcf[KKo]*tau       
        return (vcf[KKo+2]**(1.0/float(deg)))*abs(sxy)    
    if(msx>=msy):
        tt=sx;rho=sy/sx;gamm=sxy/sx;gamm=gamm*gamm
        #kcfA=2;kcfB=kcfA+nCoeff 
        cc=vcf[2:nCoeff+2]
    else: 
        tt=sy;rho=sx/sy;gamm=sxy/sy;gamm=gamm*gamm
        kcfA=2+6*nCoeff;kcfB=kcfA+nCoeff 
        cc=vcf[kcfA:kcfB]  
    yf=cc[KKo]
    jj,jj2=nMon,nMon
    while(jj>0):
        KKo-=1
        bb=cc[KKo]
        jj-=1
        mm=deg-2*jj ##degree of multiplying polynomial in sx and sy
        while(mm>0):
            mm-=1
            KKo-=1
            bb=cc[KKo]+bb*rho
        yf=bb+yf*gamm  
    #zyf=(vcf[1]*yf)**odeg  ##in case of coefficients normalization 
    #zyf=yf**odeg
    #zyf=att*zyf    
    return abs(tt)*(yf**(1.0/float(deg)))


'''
fGYF: Not vectorized yet
'''
def fGYF(sx,sy,sxy,vcf):
    deg=int(vcf[0])
    odeg=1.0/float(deg)
    nMon=int(deg/2)
    nCoeff=int(vcf[1])
    KKo=nCoeff
    KKo-=1
    msx,msy=abs(sx),abs(sy)
    if(msx+msy<gTol):
        KKo+=2
        tau=abs(sxy)
        yf=(vcf[KKo]**odeg)*tau 
        tau=sxy/tau ##( = +/- 1)
        ##gyf=[d/dsx,d/dsy,d/dsxy]     
        return yf,[0.0,0.0,tau]
    if(msx>=msy):
        tt=sx;rho=sy/sx;gamm=sxy/sx;gamm=gamm*gamm
        kcfA=2;kcfB=kcfA+nCoeff 
        cc=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d1=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d2=vcf[kcfA:kcfB]
    else: 
        tt=sy;rho=sx/sy;gamm=sxy/sy;gamm=gamm*gamm
        kcfA=2+6*nCoeff;kcfB=kcfA+nCoeff 
        cc=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d1=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d2=vcf[kcfA:kcfB]
    d3=np.zeros(nMon)   
    yf,g1yf,g2yf=cc[KKo],d1[KKo],d2[KKo]
    jj,jj2=nMon,nMon
    while(jj>0):
        KKo-=1
        bb=cc[KKo]
        d1bb,d2bb=d1[KKo],d2[KKo]
        jj-=1
        mm=deg-2*jj ##degree of multiplying polynomial in sx and sy
        while(mm>0):
            mm-=1
            KKo-=1
            bb=cc[KKo]+bb*rho
            d1bb=d1[KKo]+d1bb*rho
            d2bb=d2[KKo]+d2bb*rho           
        d3[jj-1]=jj*bb
        yf=bb+yf*gamm
        g1yf=d1bb+g1yf*gamm
        g2yf=d2bb+g2yf*gamm  
    #zyf=(vcf[1]*yf)**odeg  ##in case of coefficients normalization 
    zyf=yf**odeg
    yval,att,y2val=zyf/(deg*yf),abs(tt),float(deg-1)/zyf
    sgntt=(tt/att)*yval
    g1yf=sgntt*g1yf
    g2yf=sgntt*g2yf  
    jj2-=1    
    g3yf=float(nMon)*cc[-1]
    while(jj2>0):
        jj2-=1
        g3yf=d3[jj2]+g3yf*gamm
    tt=(2.0*sxy)/tt    
    g3yf=sgntt*tt*g3yf
    zyf=att*zyf    
    return zyf,[g1yf,g2yf,g3yf]
    
    



'''
input data: 
data = [{q,theta,sigma,r-value},[...],]
Note: input theta is in degrees
Note: for sigma or r-value use a '*' to indicate a missing value
wR = overall weight for r-values
'''
def procData(data,wR):
    if(wR>1.0 or wR<0.0):print('wR must be in [0,1]\nCalculations aborted');exit()
    vData,nsig,nrva=[],0,0
    rad=np.pi/180.0
    for item in data:
        if(item['s'] !='*'):
            dd=[0,item['q'],rad*item['theta'],item['s'],1.0]
            nsig+=1
            vData.append(dd)
        if(item['r'] !='*'):
            dd=[1,item['q'],rad*item['theta'],item['r'],1.0]
            nrva+=1
            vData.append(dd)
    wS=1.0-wR        
    for kk in range(len(vData)):
        if(vData[kk][0]):###r-value data 
            vData[kk][-1]=wR/nrva
        else: #### stress data 
            vData[kk][-1]=wS/nsig
    return vData

        
    
def fixedParam(degree,r0,s90,r90):
    vFidx=np.array([0,1,degree-1,degree],dtype=np.int)
    vFcoeff=np.ones(4,dtype=np.double)
    vFcoeff[1]=(-degree*r0)/(1+r0)
    vFcoeff[2]=(-degree*r90)/((1+r90)*s90)
    vFcoeff[3]=1.0/s90
    return vFidx,vFcoeff 


def projXzero(vFidx,fixedP,ndim):
    #vp,vnn=np.zeros((1,ndim)),np.zeros((1,ndim))
    #norm=np.sqrt(fixedP.shape[0]-1)
    #for jj in range(fixedP.shape[0]):
    #    vp[0,vFidx[jj]]=fixedP[jj]  #;print(vp)
    #    vnn[0,vFidx[jj]]=1.0/norm
    #vnn[0,0]=0.0    
    vM=np.array([1.0, -2, 3, -2, 1, 6, -6, 6, 9],dtype=np.double).reshape((1,ndim))
    for jj in range(fixedP.shape[0]):
        vM[0,vFidx[jj]]=fixedP[jj]    
    #lbd=np.sum((vp-vM)*vnn,axis=1)[0];print(lbd);print((vp-vM)*vnn)
    return vM   ####vM+lbd*vnn

'''
fCost=0.5*(MM*x)*x + VV*x
vData = [[0 or 1, q, theta, r-value or sigma, weight],[...],]
vFidx = subset of fixed parameters (calculated directly from data: s0,r0,s90,r90)  
'''
def fCostMatrx(vData,vMonoms,vFidx,vFcoeff):
    degree,nDIM,nFcoeff=vMonoms['nQ'],vMonoms['nMons'],vFidx.shape[0]
    ndim=nDIM-nFcoeff
    MM,VV=np.zeros((ndim,ndim)),np.zeros((1,ndim))
    vIdx=np.zeros(ndim,dtype=int)
    ii=0
    for k in range(nDIM):
        if(k in vFidx):continue
        vIdx[ii]=k;ii+=1
    #print('vIdx:',vIdx)
    di,vvi=0.0,np.zeros((1,nDIM))
    vQ=vMonoms['vQ']
    vD1,vD2,vD3=vMonoms['vD1Q'],vMonoms['vD2Q'],vMonoms['vD3Q']
    vC1,vC2,vC3=vMonoms['vC1Q'],vMonoms['vC2Q'],vMonoms['vC3Q']
    D1,D2,D3=np.zeros(nDIM),np.zeros(nDIM),np.zeros(nDIM)
    eps=1.0e-8
    #print(vData)
    for item in vData:
        #if(-eps<item[2]<eps or 0.5*np.pi-eps<item[2]<0.5*np.pi+eps):continue ##skip biaxial data (s0,r0,s90,r90)
        cphi,sphi=np.cos(item[2]),np.sin(item[2])
        cphi2,sphi2,csphi=cphi*cphi,sphi*sphi,cphi*sphi
        sx,sy,sxy=cphi2+item[1]*sphi2,sphi2+item[1]*cphi2,(1-item[1])*csphi
        if(item[0]):##r-value equation
            for k in range(0,nDIM):
                D1[k]=vC1[k]*(sx**vD1[k][0])*(sy**vD1[k][1])*(sxy**vD1[k][2])
                D2[k]=vC2[k]*(sx**vD2[k][0])*(sy**vD2[k][1])*(sxy**vD2[k][2])
                D3[k]=vC3[k]*(sx**vD3[k][0])*(sy**vD3[k][1])*(sxy**vD3[k][2])
            vvi[0,:]=(item[3]+sphi2)*D1[:]+(item[3]+cphi2)*D2[:]-csphi*D3[:]
            for k in range(nFcoeff):
                di-=vFcoeff[k]*vvi[0,vFidx[k]]
            MM[:,:]+=(item[4]*vvi[0,vIdx].reshape((ndim,1)))*vvi[0,vIdx].reshape((1,ndim))
            VV[:]+=item[4]*di*vvi[0,vIdx].reshape((1,ndim))
            di=0.0
        else:
            for k in range(0,nDIM):
                vvi[0,k]=(sx**vQ[k][0])*(sy**vQ[k][1])*(sxy**vQ[k][2])
            di=1.0/(item[3]**degree)
            for k in range(nFcoeff):
                di-=vFcoeff[k]*vvi[0,vFidx[k]]
            MM[:,:]+=(item[4]*vvi[0,vIdx].reshape((ndim,1)))*vvi[0,vIdx].reshape((1,ndim))
            VV[:]+=item[4]*di*vvi[0,vIdx].reshape((1,ndim))
            di=0.0
    #print(VV);print(MM)
    return MM,VV
    
    

'''
vParam=[MM,VV]
Note: it is assumed that VV and va have shape (1,ndim)
'''
def fCostFunc(va,vParam):
    return 0.5*np.sum(va*np.matmul(va,vParam[0]),axis=1)[0]-np.sum(va*vParam[1],axis=1)[0]    
    
def fCostGrad(va,vParam):
    return np.matmul(va,vParam[0])-vParam[1]
    
    


def fTestCost(data,basicData,wR,degree,vp):
    vData=procData(data,wR)
    vMonoms=vPolyMons(degree)
    r0,s90,r90=basicData[:]
    vFidx,vFcoeff=fixedParam(degree,r0,s90,r90)
    ##print('fixedParam = ',vFcoeff)
    X0=projXzero(vFidx,vFcoeff,vMonoms['nMons'])
    print('projected X0 = ',X0)
    MM,VV=fCostMatrx(vData,vMonoms,vFidx,vFcoeff)
    va=np.zeros((1,vMonoms['nMons']-vFidx.shape[0]))
    jj=0
    for k in range(vMonoms['nMons']):
        if(k in vFidx):continue
        va[0,jj]=vp[k];jj+=1
    print('fCost(vp) = ',fCostFunc(va,[MM,VV]))
    print('grad of fCost(vp) = ',fCostGrad(va,[MM,VV]))
    jj=0
    for k in range(vMonoms['nMons']):
        if(k in vFidx):continue
        va[0,jj]=X0[0,k];jj+=1
    print('fCost(proj_X0) = ',fCostFunc(va,[MM,VV]))
    print('grad of fCost(proj_X0) = ',fCostGrad(va,[MM,VV]))


        



def centerMises(degree):
    nDIM=int(degree/2)+1;nDIM*=nDIM
    if(nDIM==16):##Mises Poly6
        ##vm=np.array([1, -3, 6, -7, 6, -3, 1, 9, -18, 27, -18, 9, 27, -27, 27, 27]).reshape((nDIM,1))
        vm=np.array([1, -3, 6, -7, 6, -3, 1, 9, -18, 27, -18, 9, 27, -27, 27, 27],dtype=np.double)
    elif(nDIM==25):##Mises Poly8
        vm=np.array([1, -4, 10, -16, 19, -16, 10, -4, 1, 12, -36, 72, -84, 72, -36, 12, 54, -108, 162, -108, 54, 108, -108, 108, 81],dtype=np.double)
    elif(nDIM==36):##Mises Poly10
        vm=np.array([1, -5, 15, -30, 45, -51, 45, -30, 15, -5, 1, 15, -60, 150, -240, 285, -240, 150, -60, 15, 90, -270, 540, -630, 540, -270, 90, 270,
                     -540, 810, -540, 270, 405, -405, 405, 243],dtype=np.double)
    elif(nDIM==49): ###Mises Poly12
        vm=np.array([1, -6, 21, -50, 90, -126, 141, -126, 90, -50, 21, -6, 1, 18, -90, 270, -540, 810, -918, 810, -540, 270, -90, 18, 135, -540, 1350, 
                     -2160, 2565, -2160, 1350, -540, 135, 540, -1620, 3240, -3780, 3240, -1620, 540, 1215, -2430, 3645, -2430, 1215, 1458, -1458, 1458, 
                     729],dtype=np.double)
    elif(nDIM==9):##Mises Poly4
        vm=np.array([1, -2, 3, -2, 1, 6, -6, 6, 9],dtype=np.double)
    else:print('Unknown option: degree must be one of {4,6,8,10,12}');exit() 
    return vm  ####np.concatenate(([degree,nDIM],vm))




def genConstraintsPoints2DOptPoints(nPoints):
    N1=int(0.25*nPoints)
    vt1=np.linspace(0,np.pi/2,N1+1)
    vc1=np.cos(vt1[1:]);vs1=np.sin(vt1[1:])
    vN2=np.int_(nPoints*vs1)     
    vP=np.zeros((np.sum(vN2)+1,3))
    vN2[:]+=1
    vP[0,2]=1.0
    jj=1
    for kk in range(0,len(vt1)-1):
        ct1=vc1[kk];st1=vs1[kk]
        N2=vN2[kk];jN2=jj+N2-1
        vt2=np.linspace(0,np.pi,N2)
        vc2=np.cos(vt2[0:N2-1]);vs2=np.sin(vt2[0:N2-1])
        vP[jj:jN2,0]=st1*vc2
        vP[jj:jN2,1]=st1*vs2
        vP[jj:jN2,2]=ct1
        jj=jN2
    return vP 
    
        


def sampleSphere3D(nRandom=10**5):
    np.random.seed(99)     
    vuRandom=np.random.normal(0,1,(nRandom,3)) ##nPoints on 2D unit sphere of 3D
    vNorm=np.sqrt(np.sum(vuRandom**2,axis=1))
    vuRandom[:,0]/=vNorm;vuRandom[:,1]/=vNorm;vuRandom[:,2]/=vNorm
    vu=genConstraintsPoints2DOptPoints(nPoints=300)
    vu=np.concatenate((vu[:,0:3],vuRandom),axis=0)
    return vu
    
        


def PolyN_GaussCheck22(vCoeff,vMonoms,vu,allKG=False):
    nPoints=vu.shape[0]
    G1,G2,G3=np.zeros(nPoints),np.zeros(nPoints),np.zeros(nPoints)
    H11,H12,H13=np.zeros(nPoints),np.zeros(nPoints),np.zeros(nPoints)
    H22,H23,H33=np.zeros(nPoints),np.zeros(nPoints),np.zeros(nPoints)
    #NN=vMonoms['nMons']
    vD11=vMonoms['vH11Q'];vC11=vMonoms['vCH11Q']
    vD12=vMonoms['vH12Q'];vC12=vMonoms['vCH12Q']
    vD13=vMonoms['vH13Q'];vC13=vMonoms['vCH13Q']
    vD22=vMonoms['vH22Q'];vC22=vMonoms['vCH22Q']
    vD23=vMonoms['vH23Q'];vC23=vMonoms['vCH23Q']
    vD33=vMonoms['vH33Q'];vC33=vMonoms['vCH33Q']
    vD1=vMonoms['vD1Q'];vD2=vMonoms['vD2Q'];vD3=vMonoms['vD3Q']
    vC1=vMonoms['vC1Q'];vC2=vMonoms['vC2Q'];vC3=vMonoms['vC3Q']
    for k in range(vMonoms['nMons']):
        G1[:]+=vCoeff[k]*vC1[k]*(vu[:,0]**vD1[k][0])*(vu[:,1]**vD1[k][1])*(vu[:,2]**vD1[k][2])
        G2[:]+=vCoeff[k]*vC2[k]*(vu[:,0]**vD2[k][0])*(vu[:,1]**vD2[k][1])*(vu[:,2]**vD2[k][2])
        G3[:]+=vCoeff[k]*vC3[k]*(vu[:,0]**vD3[k][0])*(vu[:,1]**vD3[k][1])*(vu[:,2]**vD3[k][2])
        H11[:]+=vCoeff[k]*vC11[k]*(vu[:,0]**vD11[k][0])*(vu[:,1]**vD11[k][1])*(vu[:,2]**vD11[k][2])
        H12[:]+=vCoeff[k]*vC12[k]*(vu[:,0]**vD12[k][0])*(vu[:,1]**vD12[k][1])*(vu[:,2]**vD12[k][2])
        H13[:]+=vCoeff[k]*vC13[k]*(vu[:,0]**vD13[k][0])*(vu[:,1]**vD13[k][1])*(vu[:,2]**vD13[k][2])
        H22[:]+=vCoeff[k]*vC22[k]*(vu[:,0]**vD22[k][0])*(vu[:,1]**vD22[k][1])*(vu[:,2]**vD22[k][2])
        H23[:]+=vCoeff[k]*vC23[k]*(vu[:,0]**vD23[k][0])*(vu[:,1]**vD23[k][1])*(vu[:,2]**vD23[k][2])
        H33[:]+=vCoeff[k]*vC33[k]*(vu[:,0]**vD33[k][0])*(vu[:,1]**vD33[k][1])*(vu[:,2]**vD33[k][2])
    H11s=H22*H33-H23*H23
    H12s=H23*H13-H12*H33
    H13s=H12*H23-H22*H13
    H22s=H11*H33-H13*H13
    H23s=H12*H13-H11*H23
    H33s=H11*H22-H12*H12
    vNorm2=G1**2+G2**2+G3**2
    vNorm=np.sqrt(vNorm2)
    G1[:]/=vNorm;G2[:]/=vNorm;G3[:]/=vNorm
    KG=(H11s*G1*G1+H22s*G2*G2+H33s*G3*G3+2.0*(H12s*G1*G2+H13s*G1*G3+H23s*G2*G3))/vNorm2
    if(allKG):return KG
    return np.min(KG)
        


epsMin,epsMax=1.0e-5,5.0e-4
def cvxOnePt22(degree,vMons,vMises,vSigma,vDir,ffPoints,ffPointsFMT):    
    va=np.zeros(vMises.shape[0])
    lbdA,lbdB,minKG=0.0,10.0,1.0
    #wt11=0
    while(True):
        #print('----wt11 = ',wt11);wt11+=1
        va[:]=vMises[:]+lbdB*vDir[:]
        minKG=PolyN_GaussCheck22(va,vMons,vSigma)
        ffPoints.write(ffPointsFMT.format(*va[:],minKG))
        if(minKG<0.0):break
        if(minKG>2.0*epsMin):lbdA=lbdB
        else:lbdA=0.5*lbdB    
        lbdB+=10.0
    #wt22=0    
    while(True):
        #print('+++++++++wt22 = ',wt22);wt22+=1
        #if(np.abs(minKG-0.0005)<=4.0e-4):
        if((minKG>=epsMin) and (minKG<epsMax)):break
        lbd=0.5*(lbdA+lbdB)
        va[:]=vMises[:]+lbd*vDir[:]
        minKG=PolyN_GaussCheck22(va,vMons,vSigma)
        ffPoints.write(ffPointsFMT.format(*va[:],minKG))
        if(minKG<epsMin):lbdB=lbd
        else:lbdA=lbd
        #print(lbdA,lbdB,minKG)    
    return va, minKG 

def cvxBdry22(degree,vMons,vMises,vSigma,vPoints,fileIdx,fileAct='w'):
    fName=gDrive+'cvxBdry_deg_'+str(degree)+'_batch_'+str(fileIdx)
    #fKG=gDrive+fName+'KG'
    fPoints=fName+'_test.txt'
    #fCurv=gDrive+fKG+'_test.txt'
    fName+='.txt'
    #fKG+='.txt'
    #ff,ffKG,ffPoints,ffCurv=open(fName,'w'),open(fKG,'w'),open(fPoints,'w'),open(fCurv,'w')
    ff,ffPoints=open(fName,fileAct),open(fPoints,fileAct)
    frmat='{},'*(vMises.shape[0]) +'{}\n'
    for kk in range(vPoints.shape[0]):
        print('kk = ',kk)
        va,KG=cvxOnePt22(degree,vMons,vMises,vSigma,vPoints[kk,:],ffPoints,frmat)
        ff.write(frmat.format(*va[:],KG))
        va,KG=cvxOnePt22(degree,vMons,vMises,vSigma,-vPoints[kk,:],ffPoints,frmat)
        ff.write(frmat.format(*va[:],KG))
    ff.close();ffPoints.close()
    print('cvxBdry22: DONE')    
    



def plotPoly(vCoeff,data,vMonoms={},saveFig=True,figDirPlot='./',saveName='myPlot',doublePlot=False,vCoeff2=[]):
    zdata={'thetaR':[],'rValue':[],'thetaS':[0.0],'sigma':[1.0]}
    for item in data:
        if(np.abs(item['q'])<1.0e-7): ###uniaxial: q=0
            if(item['s'] != '*'):
                zdata['thetaS'].append(item['theta'])
                zdata['sigma'].append(item['s'])
            if(item['r'] != '*'):
                zdata['thetaR'].append(item['theta'])
                zdata['rValue'].append(item['r'])    
    degree=vMonoms['nQ'];oDeg=1.0/float(degree)
    NN=vMonoms['nMons']
    Nq=4
    vq=np.linspace(0,-0.25,Nq) ###; print(vq)
    if(doublePlot):
        Nq=1;vq=np.linspace(0,-0.25,Nq)
    vLW=np.ones(Nq);vLW[0]=2
    Nphi=101
    vphi=np.linspace(0,np.pi/2,Nphi);vDegs=(180/np.pi)*vphi
    vcos=np.cos(vphi);vsin=np.sin(vphi)
    vcos2,vsin2,vsc=vcos*vcos,vsin*vsin,vcos*vsin
    vsx,vsy,vsxy=np.zeros(Nphi),np.zeros(Nphi),np.zeros(Nphi)
    vS,vR, mxS,mnS = np.zeros(Nphi),np.zeros(Nphi), -1.0, 100.0
    vQ=vMonoms['vQ']
    vD1=vMonoms['vD1Q'];vD2=vMonoms['vD2Q'];vD3=vMonoms['vD3Q']
    vC1=vMonoms['vC1Q'];vC2=vMonoms['vC2Q'];vC3=vMonoms['vC3Q']
    vDX,vDY,vDXY=np.zeros(Nphi),np.zeros(Nphi),np.zeros(Nphi)
    if(doublePlot):
        vS2,vR2 = np.zeros(Nphi),np.zeros(Nphi)
        vDX2,vDY2,vDXY2=np.zeros(Nphi),np.zeros(Nphi),np.zeros(Nphi)
    fig,fig2=plt.figure(),plt.figure()
    ax,ax2=fig.add_subplot(),fig2.add_subplot()
    for kq in range(Nq):
        q=vq[kq]
        vsx[:]=vcos2+q*vsin2
        vsy[:]=vsin2+q*vcos2
        vsxy[:]=(1.0-q)*vsc
        if(doublePlot):
            for k in range(NN):
                vS[:]+=vCoeff[k]*(vsx**vQ[k][0])*(vsy**vQ[k][1])*(vsxy**vQ[k][2])
                vDX[:]+=vCoeff[k]*vC1[k]*(vsx**vD1[k][0])*(vsy**vD1[k][1])*(vsxy**vD1[k][2])
                vDY[:]+=vCoeff[k]*vC2[k]*(vsx**vD2[k][0])*(vsy**vD2[k][1])*(vsxy**vD2[k][2])
                vDXY[:]+=vCoeff[k]*vC3[k]*(vsx**vD3[k][0])*(vsy**vD3[k][1])*(vsxy**vD3[k][2])
                vS2[:]+=vCoeff2[k]*(vsx**vQ[k][0])*(vsy**vQ[k][1])*(vsxy**vQ[k][2])
                vDX2[:]+=vCoeff2[k]*vC1[k]*(vsx**vD1[k][0])*(vsy**vD1[k][1])*(vsxy**vD1[k][2])
                vDY2[:]+=vCoeff2[k]*vC2[k]*(vsx**vD2[k][0])*(vsy**vD2[k][1])*(vsxy**vD2[k][2])
                vDXY2[:]+=vCoeff2[k]*vC3[k]*(vsx**vD3[k][0])*(vsy**vD3[k][1])*(vsxy**vD3[k][2])
        else:
            for k in range(NN):
                vS[:]+=vCoeff[k]*(vsx**vQ[k][0])*(vsy**vQ[k][1])*(vsxy**vQ[k][2])
                vDX[:]+=vCoeff[k]*vC1[k]*(vsx**vD1[k][0])*(vsy**vD1[k][1])*(vsxy**vD1[k][2])
                vDY[:]+=vCoeff[k]*vC2[k]*(vsx**vD2[k][0])*(vsy**vD2[k][1])*(vsxy**vD2[k][2])
                vDXY[:]+=vCoeff[k]*vC3[k]*(vsx**vD3[k][0])*(vsy**vD3[k][1])*(vsxy**vD3[k][2])                    
        vR[:]=(vDXY*vsc-vDX*vsin2-vDY*vcos2)/(vDX+vDY)
        vS[:]=1.0/vS**oDeg
        mxS=max(max(vS),mxS);mnS=min(min(vS),mnS)
        if(doublePlot):
            vR2[:]=(vDXY2*vsc-vDX2*vsin2-vDY2*vcos2)/(vDX2+vDY2)
            vS2[:]=1.0/vS2**oDeg
            mxS=max(max(vS2),mxS);mnS=min(min(vS2),mnS)
        ax.plot(vDegs,vR,linewidth=vLW[kq],color='k')
        ax2.plot(vDegs,vS,linewidth=vLW[kq],color='k')
        if(doublePlot):
            ax.plot(vDegs,vR2,linewidth=vLW[kq],color='k',linestyle='--')
            ax2.plot(vDegs,vS2,linewidth=vLW[kq],color='k',linestyle='--')
        vS[:],vDX[:],vDY[:],vDXY[:]=0.0,0.0,0.0,0.0
    ax.plot(zdata['thetaR'],zdata['rValue'],linestyle='',marker='o',markersize=6,markerfacecolor='b')
    x0,x1=-0.25,90.25
    y0,y1=ax.get_ylim()
    ax.text(x1-0.05*(x1-x0),y0+0.01*(y1-y0),r'$\theta$',fontsize=16)
    if(doublePlot):
        ax.text(x0+0.01*(x1-x0),y1-0.05*(y1-y0),r'$r(\theta)$',fontsize=16)
    else:
        ax.text(x0+0.01*(x1-x0),y1-0.05*(y1-y0),r'$r(q,\theta)$',fontsize=16) 
    majorXticks=[0,15,30,45,60,75,90]
    minorXticks=np.linspace(0,90,4*6+1)
    ax.set_xticks(minorXticks,minor=True)
    ax.set_xticks(majorXticks,minor=False)
    ax.set_xticklabels([r'$0^o$','$15^o$','$30^o$','$45^o$','$60^o$','$75^o$','$90^o$'])
    #ax.grid()
    yMin=int(min(zdata['rValue']))
    yMax=int(max(zdata['rValue']))+1
    majorYticks=np.linspace(yMin,yMax,yMax-yMin+1)
    #minorYticks=np.linspace(yMin,yMax,4*(yMax-yMin)+1)
    #ax.set_yticks(majorYticks,minor=False)
    #ax.set_yticks(minorYticks,minor=True)
    #ax.set_yticklabels([str(round(k,2)) for k in minorYticks],minor=True,fontsize=10)
    #ax.set_yticklabels([str(round(k,2)) for k in ax.get_yticks(minor=False)],minor=False,fontsize=10)
    #ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.65)
    ax.set_xlim(-0.85,90.85)
    ax.set_ylim(yMin-0.0275,0.9*yMax)
    fName=saveName+'_P{}'.format(str(degree))
    if(saveFig):
        fig.savefig(figDirPlot+fName+'_rValue.png',bbox_inches='tight',dpi=300)
    ax2.plot(zdata['thetaS'],zdata['sigma'],linestyle='',marker='o',markersize=6,markerfacecolor='b')
    x0,x1=-0.25,90.25
    ax2.set_xlim(x0,x1)
    y0,y1=ax2.get_ylim()
    y1+=0.01
    ax2.set_ylim((y0,y1))
    ax2.text(x1-0.05*(x1-x0),y0+0.01*(y1-y0),r'$\theta$',fontsize=16)
    if(doublePlot):
        ax2.text(x0+0.01*(x1-x0),y1-0.05*(y1-y0),r'$\overline{\sigma}(\theta)$',fontsize=16)
    else:
        ax2.text(x0+0.01*(x1-x0),y1-0.05*(y1-y0),r'$\overline{\sigma}(q,\theta)$',fontsize=16)    
    ax2.set_xticks([0,15,30,45,60,75,90])
    ax2.set_xticklabels(['$0^o$','$15^o$','$30^o$','$45^o$','$60^o$','$75^o$','$90^o$'])
    ax2.grid()
    if(saveFig):
        fig2.savefig(figDirPlot+fName+'_Sigma.png',bbox_inches='tight',dpi=300)    
    ###exit()
    ###########################################################################################################
    fig.clf()
    ax=fig.add_subplot()
    NX,NY,maxX,maxY=201,201,1.25,1.25
    vsx=np.linspace(-maxX,maxX,NX)
    vsy=np.linspace(-maxY,maxY,NY)
    maxSXY=1.0/(vCoeff[-1]**oDeg)
    maxSXY2=1.0/(vCoeff2[-1]**oDeg)
    #vsxy=[0.0]
    vsxy=maxSXY*np.array([0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.975,0.99,0.998])
    vsxy2=maxSXY2*np.array([0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.975,0.99,0.998])
    nIdx=int(degree/2+1)
    vTerms=np.zeros((nIdx,NY,NX))
    if(doublePlot):
        vTerms2=np.zeros((nIdx,NY,NX))
        for kidx in range(nIdx):
            vPowers=vMonoms['vQidx'][kidx][1]#; print(vPowers)
            #nPowers=len(vPowers)   
            for ky in range(NY):
                for i in vPowers:
                    vTerms[kidx,ky,:]+=vCoeff[i]*(vsx**(vMonoms['vQ'][i][0]))*(vsy[ky]**(vMonoms['vQ'][i][1]))
                    vTerms2[kidx,ky,:]+=vCoeff2[i]*(vsx**(vMonoms['vQ'][i][0]))*(vsy[ky]**(vMonoms['vQ'][i][1]))
        #for sxy in vsxy:
        for kkxy in range(len(vsxy)):
            sxy=vsxy[kkxy]
            sxy2=sxy*sxy
            sxy22=vsxy2[kkxy]*vsxy2[kkxy]
            vz=vTerms[nIdx-1]
            vz2=vTerms2[nIdx-1]
            for i in range(nIdx-2,-1,-1):
                vz=vTerms[i]+sxy2*vz
                vz2=vTerms[i]+sxy22*vz2
            #vz=vz**(1.0/degree)
            ax.contour(vsx,vsy,vz,levels=[1.0],linewidths=1,colors=['k'])
            ax.contour(vsx,vsy,vz2,levels=[1.0],linewidths=1,colors=['k'],linestyles='--')
            #vz2=vTerms[0]**(1.0/degree);print(np.max(vz-vz2),np.min(vz-vz2))
            #ax.contour(vsx,vsy,vz2,levels=[1.0],linewidths=1,colors=['r'])
    else:        
        for kidx in range(nIdx):
            vPowers=vMonoms['vQidx'][kidx][1]#; print(vPowers)
            #nPowers=len(vPowers)   
            for ky in range(NY):
                for i in vPowers:
                    vTerms[kidx,ky,:]+=vCoeff[i]*(vsx**(vMonoms['vQ'][i][0]))*(vsy[ky]**(vMonoms['vQ'][i][1]))
        for sxy in vsxy:
            sxy2=sxy*sxy
            vz=vTerms[nIdx-1]
            for i in range(nIdx-2,-1,-1):
                vz=vTerms[i]+sxy2*vz#;print(sxy2,i)
            #vz=vz**(1.0/degree)
            ax.contour(vsx,vsy,vz,levels=[1.0],linewidths=1,colors=['k'])
            #vz2=vTerms[0]**(1.0/degree);print(np.max(vz-vz2),np.min(vz-vz2))
            #ax.contour(vsx,vsy,vz2,levels=[1.0],linewidths=1,colors=['r'])
    ax.grid()
    ax.set_aspect('equal')
    x0,x1=ax.get_xlim()
    y0,y1=ax.get_ylim()
    ax.text(x1-0.075*(x1-x0),y0+0.02*(y1-y0),r'$\overline{\sigma}_x$',fontsize=14)
    ax.text(x0+0.01*(x1-x0),y1-0.075*(y1-y0),r'$\overline{\sigma}_y$',fontsize=14)
    if(saveFig):
        fig.savefig(figDirPlot+fName+'_sxySections.png',bbox_inches='tight',dpi=300)
    else:
        plt.show()





def plotKG(vCoeff,vMonoms,vu,saveFig=True,figDirPlot='./',saveName='myPlot'):
    vQ=vMonoms['vQ']
    vS=np.zeros((vu.shape[0],1))
    for k in range(vMonoms['nMons']):
        vS[:,0]+=vCoeff[k]*(vu[:,0]**vQ[k][0])*(vu[:,1]**vQ[k][1])*(vu[:,2]**vQ[k][2])
    vS[:,0]=vS[:,0]**(1.0/vMonoms['nQ'])
    vKG=PolyN_GaussCheck22(vCoeff,vMonoms,vu,allKG=True)
    idx=np.argwhere(vKG<0)
    vKG=np.ones(vKG.shape[0])
    vKG[idx[:,0]]=-1.0
    fg3=plt.figure(figsize=(8,8))
    ax=fg3.add_subplot()
    ax.set_aspect('equal')
    ax.spines[['right', 'top']].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=16)
    spl=ax.scatter(vu[:,0]/vS[:,0],vu[:,1]/vS[:,0],c=vKG, cmap='seismic')
    cbar=fg3.colorbar(spl,ax=ax, shrink=0.6,pad=0.0075,location='right')
    #cbar.ax.set_ylim(-2.0, 10.0)
    #cbar.outline.set_visible(False)
    #vticks=np.linspace(np.min(KG),np.max(KG),6)
    #cbar.set_ticks(vticks)
    #cbar.set_ticklabels([str(round(x,2)) for x in vticks],fontsize=11)
    vticks=[-0.5,0.5]
    cbar.set_ticks(vticks)
    cbar.set_ticklabels([r'$\mathregular{<0}$',r'$\mathregular{\geq 0}$'],fontsize=18)
    ##cbar.ax.set_title(r'$\mathregular{\sigma_{xy}}$',fontsize=15)
    cbar.ax.set_title(r'$\,\,\mathregular{K_G}$',fontsize=18)
    if(saveFig):fg3.savefig(saveName+'_KG.png',bbox_inches='tight',dpi=300,facecolor='w')
    else:plt.show()

        
############################   Optimization section ------------------------------------------ 
'''
xZero = initial guess (assumed to be within the feasible set)
fFunc, gradF = the cost function to be minimized and its gradient
gFunc, gradG = the constraint mapping the boundary of the feasible domain (gFunc(x)<=0.0) 
                  and its gradient
fParam, gParam = parameters of fFunc and gFunct (if any)
dLambda = the scaling factor of the distance along a gradient descent step
'''




def constrOptim33(xZero,xCenter,fFunc,gradF,gFunc,gradG,fParam=[],gParam=[],dLambda=0.001,
               nIter=10**4,epsDX=1.0e-7,epsG=1.0e-7,epsGF=1.0e-3):
    nDim=xZero.shape[0]
    x0=np.zeros(nDim)
    x0[:]=xZero[:]
    xA,xAA,xB,gF=np.zeros(nDim),np.zeros(nDim),np.zeros(nDim),np.zeros(nDim)
    fF,fA,fAA,fB,gA,gZer=0.0,0.0,0.0,0.0,0.0,0.0
    flagExit,flagG=False,False
    xMid,grdG,normG,normF=np.zeros(nDim),np.zeros(nDim),1.0,1.0
    kIter=0
    while(kIter<nIter):
        fF,gF[:]=fFunc(x0,fParam),-gradF(x0,fParam)
        if(flagG):###project gradient 
            grdG[:]=gradG(x0,gParam)
            normG=np.sqrt(np.sum(grdG*grdG))
            grdG[:]/=normG
            gF[:]=gF-np.sum(gF*grdG)*grdG
        normF=np.sqrt(np.sum(gF*gF,axis=0))  ##; print('normF = ',normF)
        if(normF<epsGF):flagExit=True;print('zero increment and gA = ',gA);break
        xA[:]=x0[:]+dLambda*gF[:]
        gA,gZer=gFunc(xA,gParam),gFunc(x0,gParam)
        if(gA>0.0):
            x0[:]=xCenter
            while(True):
                xMid[:]=x0[:]+0.8*(xA[:]-x0[:])
                gA=gFunc(xMid,gParam) #;print(gA,gZer,normF)
                if(gA>=0.0):xA[:]=xMid[:]
                else:x0[:]=xMid[:]
                if(np.sum(np.sqrt((xA-x0)**2))<1.0e-7):xA[:]=x0[:];break
                #kW+=1
            flagG=True
        else: flagG=False                
        x0[:]=xA[:]    
        kIter+=1  
        if(kIter%500==0):
            print(f'iter = {kIter}, funcVal = {fF}, gradNorm = {normF}, gA = {gA}, gZero = {gZer}, flagG = {flagG}')
            print(f'x0 = {x0}')
    if(not flagExit):print('reached max number of iterations in the constrained part of the algorithm')    
    return fF,x0
    
    
def constrOptim33B(xZero,xCenter,fFunc,gradF,gFunc,gradG,fParam=[],gParam=[],dLambda=0.001,
               nIter=10**4,epsDX=1.0e-7,epsG=1.0e-7,epsGF=1.0e-3):
    nDim=xZero.shape[0]
    x0=np.zeros(nDim)
    x0[:]=xZero[:]
    xA,xAA,xB,gF=np.zeros(nDim),np.zeros(nDim),np.zeros(nDim),np.zeros(nDim)
    fF,fA,fAA,fB,gA,gZer=0.0,0.0,0.0,0.0,0.0,0.0
    flagExit,flagG=False,False
    xMid,grdG,normG,normF=np.zeros(nDim),np.zeros(nDim),1.0,1.0
    vLd=np.linspace(0.0,dLambda,10)
    rnld=range(1,vLd.shape[0]-1)
    kIter=0
    while(kIter<nIter):
        gF[:]= -gradF(x0,fParam)
        if(flagG):###project gradient 
            grdG[:]=gradG(x0,gParam)
            normG=np.sqrt(np.sum(grdG*grdG))
            grdG[:]/=normG
            gF[:]=gF-np.sum(gF*grdG)*grdG
        normF=np.sqrt(np.sum(gF*gF,axis=0))  ##; print('normF = ',normF)
        if(normF<epsGF):flagExit=True;print('zero increment and gA = ',gA);break
        for kj in rnld:
            xAA[:]=x0[:]+vLd[kj]*gF[:]
            xA[:]=x0+vLd[kj+1]*gF
            fA,fAA=fFunc(xA,fParam),fFunc(xAA,fParam)
            if(fAA<=fA):xA[:]=xAA[:];break
        gA,gZer=gFunc(xA,gParam),gFunc(x0,gParam)
        if(gA>0.0):
            x0[:]=xCenter
            while(True):
                xMid[:]=x0[:]+0.8*(xA[:]-x0[:])
                gA=gFunc(xMid,gParam) #;print(gA,gZer,normF)
                if(gA>=0.0):xA[:]=xMid[:]
                else:x0[:]=xMid[:]
                if(np.sum(np.sqrt((xA-x0)**2))<1.0e-7):xA[:]=x0[:];break
                #kW+=1
            flagG=True
        else: flagG=False                
        x0[:]=xA[:]    
        kIter+=1  
        if(kIter%500==0):
            print(f'iter = {kIter}, funcVal = {fA}, gradNorm = {normF}, gA = {gA}, gZero = {gZer}, flagG = {flagG}')
            print(f'x0 = {x0}')
    if(not flagExit):print('reached max number of iterations in the constrained part of the algorithm')    
    return fF,x0  
    
def constrOptim33C(xZero,xCenter,fFunc,gradF,gFunc,gradG,fParam=[],gParam=[],dLambda=0.001,
               nIter=10**4,epsDX=1.0e-7,epsG=1.0e-7,epsGF=5.0e-3):
    nDim=xZero.shape[0]
    x0=np.zeros(nDim)
    x0[:]=xZero[:]
    xA,xAA,xB,gF=np.zeros(nDim),np.zeros(nDim),np.zeros(nDim),np.zeros(nDim)
    fF,fA,fAA,fB,gA,gZer=0.0,0.0,0.0,0.0,0.0,0.0
    flagExit,flagG=False,False
    xMid,grdG,normG,normF=np.zeros(nDim),np.zeros(nDim),1.0,1.0
    vLd=np.linspace(0.0,dLambda,10)
    rnld=range(1,vLd.shape[0]-1)
    kIter=0
    while(kIter<nIter):
        gF[:]= -gradF(x0,fParam)
        if(flagG):###project gradient 
            grdG[:]=gradG(x0,gParam)
            normG=np.sqrt(np.sum(grdG*grdG))
            grdG[:]/=normG
            gF[:]=gF-np.sum(gF*grdG)*grdG
        for kj in rnld:
            xAA[:]=x0[:]+vLd[kj]*gF[:]
            xA[:]=x0+vLd[kj+1]*gF
            fA,fAA=fFunc(xA,fParam),fFunc(xAA,fParam)
            if(fAA<=fA):xA[:]=xAA[:];break
        gA=gFunc(xA,gParam)
        normF=np.sqrt(np.sum(gF*gF,axis=0))  ##; print('normF = ',normF)
        if((normF<epsGF) and (gA<=0.0)):flagExit=True;print('zero increment and gA = ',gA);break
        if(gA>0.0):
            x0[:]=xCenter
            while(True):
                xMid[:]=x0[:]+0.8*(xA[:]-x0[:])
                gA=gFunc(xMid,gParam) #;print(gA,gZer,normF)
                if(gA>=0.0):xA[:]=xMid[:]
                else:x0[:]=xMid[:]
                if(np.sum(np.sqrt((xA-x0)**2))<1.0e-7):xA[:]=x0[:];break
                #kW+=1
            flagG=True
        else: flagG=False                
        x0[:]=xA[:]    
        kIter+=1  
        if(kIter%500==0):
            print(f'iter = {kIter}, funcVal = {fA}, gradNorm = {normF}, gA = {gA},  flagG = {flagG}')
            print(f'x0 = {x0}')
    if(not flagExit):print('reached max number of iterations in the constrained part of the algorithm')    
    return fF,x0      


def constrOptim33CD(xZero,xCenter,fFunc,gradF,gFunc,gradG,fParam=[],gParam=[],dLambda=0.001,
               nIter=10**4,epsDX=1.0e-7,epsG=1.0e-7,epsGF=1.0e-3):
    nDim=xZero.shape[0]
    x0=np.zeros(nDim)
    x0[:]=xZero[:]
    xA,xAA,xB,gF=np.zeros(nDim),np.zeros(nDim),np.zeros(nDim),np.zeros(nDim)
    fF,fA,fAA,fB,gA,gZer=0.0,0.0,0.0,0.0,0.0,0.0
    flagExit,flagG=False,False
    xMid,grdG,normG,normF=np.zeros(nDim),np.zeros(nDim),1.0,1.0
    solMin={'xMin':[],'fMin':9999.0,'iter':100,'gFnorm':9999.0,'gA':111.0}
    vLd=np.linspace(0.0,dLambda,4)
    rnld=range(0,vLd.shape[0]-1)
    kIter=0
    while(kIter<nIter):
        fF,gF[:]= fFunc(x0,fParam), -gradF(x0,fParam)
        if(flagG):###project gradient 
            grdG[:]=gradG(x0,gParam)
            normG=np.sqrt(np.sum(grdG*grdG))
            grdG[:]/=normG
            gF[:]=gF-np.sum(gF*grdG)*grdG
        xA[:]=x0+dLambda*gF
        #for kj in rnld:
        #    xAA[:]=x0[:]+vLd[kj]*gF[:]
        #    xA[:]=x0+vLd[kj+1]*gF
        #    fA,fAA=fFunc(xA,fParam),fFunc(xAA,fParam)
        #    if(fAA<=fA):xA[:]=xAA[:];break
        gA=gFunc(xA,gParam)
        normF=np.sqrt(np.sum(gF*gF,axis=0))  ##; print('normF = ',normF)
        if((normF<epsGF) and (gA<=0.0)):
            flagExit=True
            print(f'iter = {kIter}, funcVal = {fA}, gradNorm = {normF}, gA = {gA}, flagG = {flagG}')
            print('zero increment and gA = ',gA);break
        #if((np.sqrt(np.sum((xA-x0)**2))<1.0e-7) and (gA<=0.0)):
        #    flagExit=True
        #    print(f'iter = {kIter}, funcVal = {fA}, gradNorm = {normF}, gA = {gA}, flagG = {flagG}')
        #    print('cost is no longer decreasing along search direction');break
        if(gA>0.0):
            xB[:]=xCenter
            while(True):
                xMid[:]=xB[:]+0.8*(xA[:]-xB[:])
                gA=gFunc(xMid,gParam) #;print(gA,gZer,normF)
                if(gA>=0.0):xA[:]=xMid[:]
                else:xB[:]=xMid[:]
                if(np.sum(np.sqrt((xA-xB)**2))<1.0e-7):xA[:]=xB[:];break
                #kW+=1
            flagG=True
        else: flagG=False
        #fA=fFunc(xA,fParam)
        #if((fA>fF) and (gA<=0.0)):
        #    flagExit=True
        #    print(f'iter = {kIter}, funcVal = {fA}, gradNorm = {normF}, gA = {gA}, flagG = {flagG}')
        #    print('cost is no longer decreasing');break     
        fA=fFunc(xA,fParam)
        if((fA<solMin['fMin']) and (gA<=0)):
            solMin['fMin'],solMin['xMin'],solMin['iter'],solMin['gFnorm'],solMin['gA']=fA,xA,kIter,normF,gA
        if(kIter>solMin['iter']+1001):
            flagExit=True
            print('Solution no longer decreasing');break               
        x0[:]=xA[:]    
        kIter+=1  
        if(kIter%100==0):
            print(f'iter = {kIter}, funcVal = {fA}, gradNorm = {normF}, gA = {gA}, flagG = {flagG}')
            print(f'x0 = {x0}')
    if(not flagExit):print('reached max number of iterations in the constrained part of the algorithm') 
    print(f'Final number of iter = {kIter}')   
    return fF,x0,solMin


if(0):
    ### test functions: Both are quadratics 
    def testF(x,vParam):
        return np.sum(vParam[1]*(x-vParam[0])**2,axis=0)
    def testGradF(x,vParam):
        return 2.0*vParam[1]*(x-vParam[0])
    def testG(x,vParam):
        return np.sum(vParam[1]*(x-vParam[0])**2,axis=0)-1.0
    def testGradG(x,vParam):
        return 2.0*vParam[1]*(x-vParam[0])
    ### vectorized forms
    def vtestF(x,vParam):
        vx=x[:,0]-vParam[0][0]
        vy=x[:,1]-vParam[0][1]
        return vParam[1][0]*vx*vx+vParam[1][1]*vy*vy    
    def vtestG(x,vParam):
        vx=x[:,0]-vParam[0][0]
        vy=x[:,1]-vParam[0][1]
        return vParam[1][0]*vx*vx+vParam[1][1]*vy*vy - 1.0 
          

    if(0):### test case: min is within the feasible domain (passed)
        vfParam=[np.array([0.9,0.2]),np.array([1.0,2.0])]
        vgParam=[np.array([0.0,0.0]),np.array([1.0,0.5])]

    if(0):### test case: min is on the boundary of the feasible domain (passed)
        vfParam=[np.array([1.00005,0.0]),np.array([1.0,2.0])]
        vgParam=[np.array([0.0,0.0]),np.array([1.0,0.5])]

    if(1):### test case: min is outside the feasible domain ()
        vfParam=[np.array([1.2,1.2]),np.array([1.0,2.0])]
        vgParam=[np.array([0.0,0.0]),np.array([1.0,0.5])]
        
        


    xZro=np.array([-0.9,0.1])
    xCenter=np.zeros(2)
    fMin,xSol=constrOptim33(xZro,xCenter,testF,testGradF,testG,testGradG,vfParam,vgParam)
    print('minF = ',fMin)
    print('xSol = ', xSol)

    if(1):### test case: min is outside the feasible domain ()
        ###test set
        ndx,sq2=1001,np.sqrt(2.0)
        vx=np.linspace(-1.0,1.0,ndx)
        vy=np.linspace(-sq2,sq2,ndx)
        X,Y=np.meshgrid(vx,vy)
        X,Y=X.reshape(ndx*ndx),Y.reshape(ndx*ndx)
        vPoints=np.zeros((ndx*ndx,2))
        vPoints[:,0],vPoints[:,1]=X,Y
        idx=np.argwhere(vtestG(vPoints,vgParam)<=0.0)
        fVals=vtestF(vPoints[idx[:,0]],vfParam)
        print('minF = ',np.min(fVals))
        ##print('maxF = ',np.max(fVals))
        idmin=np.argmin(fVals)
        print('idmin = ',idmin)
        print(vPoints[idx[idmin,0],:],  'g(xMin) = ', testG(vPoints[idx[idmin,0]],vgParam),',  ', vtestG(vPoints[idx[idmin,0]].reshape((1,2)),vgParam))
        print('minF = ',fVals[idmin])
        fg=plt.figure();ax=fg.add_subplot()
        ax.scatter(vPoints[idx[:,0],0],vPoints[idx[:,0],1],alpha=0.5)
        ax.plot(vPoints[idx[idmin,0],0],vPoints[idx[idmin,0],1],color='r',marker='o',markersize=9,alpha=0.5)
        ax.plot(xSol[0],xSol[1],color='k',marker='o',markersize=7)
        ax.grid()
        ax.set_aspect('equal')
        fg.savefig('feasibleDomain.png',format='png',dpi=300,bbox_inches='tight',facecolor='w')

    
        

import numpy as np 

zro=0.0
zroTol=1.0e-10
gTol=1.0e-9

## Hill'48 yield function
## Note: a0=1 (normalized)
def ysHill(sx,sy,sxy,a1,a2,a3):
    return np.sqrt(sx*sx+a1*sx*sy+a2*sy*sy+a3*sxy*sxy)

def fGYF_Hill(sx,sy,sxy,vcf):
    a1,a2,a3=vcf[:]
    yf=np.sqrt(sx*sx+a1*sx*sy+a2*sy*sy+a3*sxy*sxy)
    dx=(2.0*sx+a1*sy)
    dy=(a1*sx+2.0*a2*sy)
    dxy=2.0*a3*sxy
    return yf,[dx,dy,dxy]
    
## Calculate Hill'48 parameters a_i from data
## Data: (normalized) yield stresses y_45/y_0 and y_90/y_0; r-value r_0
## Data: or the three r-values: r_0, r_45, r_90
def paramHill(vData, rval=True):
    if(rval):
        r0,r45,r90=vData
        a1=-2*r0/(1+r0)
        a2=-a1*(1+r90)/(2.0*r90)
        a3=(1.0+a1+a2)*(1.0+2.0*r45)
    else:
        y45,y90,r0=vData
        a1=-2*r0/(1+r0)
        a2=1/(y90*y90)
        a3=4/(y45*y45)-(1+a1+a2)
    return a1,a2,a3    




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
    

def nMonoms(degQ):
    deg=degQ//2    
    return (deg+1)**2,deg*(deg+1)

def vPoly(degree):
    nQ=int(degree)
    if(nQ%2):print("'degree' must be an even integer\nCalculations aborted");exit()
    dd={'nQ':nQ}
    vP=[];vPidx=[];jj=0
    while(jj<=nQ-1):
        lv=len(vP)
        vPidx.append([jj,[k for k in range(lv,lv+nQ-jj)]])
        vP+=[(nQ-1-jj-k,k,jj) for k in range(nQ-jj)]
        jj+=2
    dd['vP']=vP ;dd['vPidx']=vPidx   
    dd['vD1P']=[(mon[0]-1,mon[1],mon[2]) if(mon[0]>0) else (0,0,0) for mon in vP]
    dd['vC1P']=np.array([mon[0] if(mon[0]>0) else 0 for mon in vP])
    dd['vD2P']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in vP]
    dd['vC2P']=np.array([mon[1] if(mon[1]>0) else 0 for mon in vP])
    dd['vD3P']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in vP]
    dd['vC3P']=np.array([mon[2] if(mon[2]>0) else 0 for mon in vP]) 
    dd['vH11P']=[(mon[0]-1,mon[1],mon[2]) if(mon[0]>0) else (0,0,0) for mon in dd['vD1P']]
    dd['vCH11P']=np.array([cf*mon[0] for (cf,mon) in zip(dd['vC1P'],dd['vD1P'])])
    dd['vH12P']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in dd['vD1P']]
    dd['vCH12P']=np.array([cf*mon[1] for (cf,mon) in zip(dd['vC1P'],dd['vD1P'])])
    dd['vH13P']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD1P']]
    dd['vCH13P']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC1P'],dd['vD1P'])])
    dd['vH22P']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in dd['vD2P']]
    dd['vCH22P']=np.array([cf*mon[1] for (cf,mon) in zip(dd['vC2P'],dd['vD2P'])])
    dd['vH23P']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD2P']]
    dd['vCH23P']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC2P'],dd['vD2P'])])
    dd['vH33P']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD3P']]
    dd['vCH33P']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC3P'],dd['vD3P'])])    
    vQ=[];vQidx=[];jj=0
    while(jj<=nQ):
        lv=len(vQ)
        vQidx.append([jj,[k for k in range(lv,lv+nQ+1-jj)]])
        vQ+=[(nQ-jj-k,k,jj) for k in range(nQ+1-jj)]
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


##nQ,nP=nMonoms(degQ)
##ddMon=vPoly(degQ)
#
def fSHYqp(vsx,vsy,vsxy,vCoeff,ddMon,nQ,nP):
    sq3=np.sqrt(3.0);sq32=np.sqrt(1.5);sq6=np.sqrt(6.0);sq2=np.sqrt(2.0)
    vP=ddMon['vP'];vPidx=ddMon['vPidx'];nPidx=len(vPidx)##;print("vPidx\n",vPidx)
    vQ=ddMon['vQ'];vQidx=ddMon['vQidx'];nQidx=len(vQidx)##;print("vQidx\n",vQidx)
    nSamples=vsx.shape[0]
    vx,vy,vz=np.zeros(nSamples),np.zeros(nSamples),np.zeros(nSamples)
    vx[:]=vsx[:];vy[:]=vsy[:];vz[:]=vsxy[:]
    #zz=np.zeros(nSamples)
    vx=(2.0*vx-vy)/sq6;vy=vy/sq2;vx2y2=vx**2+vy**2       
    vMod=np.zeros(nSamples)
    vux=np.zeros(nSamples);vuy=np.zeros(nSamples);vuz=np.zeros(nSamples)
    vMonP=np.zeros((nSamples,nP));vTermsP=np.zeros((nSamples,nPidx))
    vMonQ=np.zeros((nSamples,nQ));vTermsQ=np.zeros((nSamples,nQidx))
    vMod[:]=np.sqrt(vx2y2+2*vz**2)
    vux[:]=vx/vMod;vuy[:]=vy/vMod;vuz[:]=(2*vz**2)/vMod**2
    for k in range(nP):
        vMonP[:,k]=vux**vP[k][0]*vuy**vP[k][1]
    for k in range(nPidx):
        vTermsP[:,k]=np.dot(vMonP[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
    vsP=vTermsP[:,-1]
    for k in range(nPidx-2,-1,-1):
        vsP[:]=vTermsP[:,k]+vuz*vsP
    for k in range(nQ):
        vMonQ[:,k]=vux**vQ[k][0]*vuy**vQ[k][1]
    for k in range(nQidx):
        vTermsQ[:,k]=np.dot(vMonQ[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])         
    vsQ=vTermsQ[:,-1]
    for k in range(nQidx-2,-1,-1):
        vsQ[:]=vTermsQ[:,k]+vuz*vsQ
    return sq32*vMod*(1.0+vsP+vsQ)
    

def fGYF_SHYqp(vsx,vsy,vsxy,vCoeff,ddMon,nQ,nP):
    degQ=ddMon['nQ']
    degQm1=degQ-1;degPm1=degQ-2
    sq3=np.sqrt(3.0);sq6=np.sqrt(6.0);sq2=np.sqrt(2.0)
    vP=ddMon['vP'];vPidx=ddMon['vPidx'];nPidx=len(vPidx)
    vD1P=ddMon['vD1P'];vD2P=ddMon['vD2P'];vD3P=ddMon['vD3P']
    vC1P=ddMon['vC1P'];vC2P=ddMon['vC2P'];vC3P=ddMon['vC3P']
    vQ=ddMon['vQ'];vQidx=ddMon['vQidx'];nQidx=len(vQidx)
    vD1Q=ddMon['vD1Q'];vD2Q=ddMon['vD2Q'];vD3Q=ddMon['vD3Q']
    vC1Q=ddMon['vC1Q'];vC2Q=ddMon['vC2Q'];vC3Q=ddMon['vC3Q']
    try:
        nSamples=vsx.shape[0]
    except:
        nSamples=1    
    vx,vy,vz=np.zeros(nSamples),np.zeros(nSamples),np.zeros(nSamples)
    vx[:]=vsx[:];vy[:]=vsy[:];vz[:]=vsxy[:]
    #zz=np.zeros(nSamples)
    vx=(2.0*vx-vy)/sq6;vy=vy/sq2;vx2y2=vx**2+vy**2       
    vMod=np.zeros(nSamples)
    vMod[:]=np.sqrt(vx2y2+2*vz**2)
    vu=np.zeros((nSamples,3))
    vu[:,0]=vx/vMod
    vu[:,1]=vy/vMod
    vu[:,2]=vz/vMod
    vu3=2.0*vu[:,2]**2
    vMon=np.zeros((nSamples,nP))
    vMonD1=np.zeros((nSamples,nP));vMonD2=np.zeros((nSamples,nP))
    vTerms=np.zeros((nSamples,nPidx))
    vTermsD1=np.zeros((nSamples,nPidx));vTermsD2=np.zeros((nSamples,nPidx));vTermsD3=np.zeros((nSamples,nPidx))
    for k in range(nP):
        vMon[:,k]=vu[:,0]**vP[k][0]*vu[:,1]**vP[k][1]
        vMonD1[:,k]=vC1P[k]*vu[:,0]**vD1P[k][0]*vu[:,1]**vD1P[k][1]
        vMonD2[:,k]=vC2P[k]*vu[:,0]**vD2P[k][0]*vu[:,1]**vD2P[k][1]
    for k in range(nPidx):
        vTerms[:,k]=np.dot(vMon[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3P[vPidx[k][1][0]]*vTerms[:,k]    
    vvP=vTerms[:,-1];vPD1=vTermsD1[:,-1];vPD2=vTermsD2[:,-1];vPD3=vTermsD3[:,-1]
    for k in range(nPidx-2,-1,-1):
        vvP[:]=vTerms[:,k]+vu3*vvP
        vPD1[:]=vTermsD1[:,k]+vu3*vPD1
        vPD2[:]=vTermsD2[:,k]+vu3*vPD2 
    for k in range(nPidx-2,0,-1):        
        vPD3[:]=vTermsD3[:,k]+vu3*vPD3
    vPD3[:]*=vu[:,2]   
    vMon=np.zeros((nSamples,nQ))
    vMonD1=np.zeros((nSamples,nQ));vMonD2=np.zeros((nSamples,nQ))
    vTerms=np.zeros((nSamples,nQidx))
    vTermsD1=np.zeros((nSamples,nQidx));vTermsD2=np.zeros((nSamples,nQidx));vTermsD3=np.zeros((nSamples,nQidx))
    for k in range(nQ):
        vMon[:,k]=vu[:,0]**vQ[k][0]*vu[:,1]**vQ[k][1]
        vMonD1[:,k]=vC1Q[k]*vu[:,0]**vD1Q[k][0]*vu[:,1]**vD1Q[k][1]
        vMonD2[:,k]=vC2Q[k]*vu[:,0]**vD2Q[k][0]*vu[:,1]**vD2Q[k][1]
    for k in range(nQidx):
        vTerms[:,k]=np.dot(vMon[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3Q[vQidx[k][1][0]]*vTerms[:,k]        
    vvQ=vTerms[:,-1];vQD1=vTermsD1[:,-1];vQD2=vTermsD2[:,-1];vQD3=vTermsD3[:,-1]
    for k in range(nQidx-2,-1,-1):
        vvQ[:]=vTerms[:,k]+vu3*vvQ
        vQD1[:]=vTermsD1[:,k]+vu3*vQD1
        vQD2[:]=vTermsD2[:,k]+vu3*vQD2
    for k in range(nQidx-2,0,-1):        
        vQD3[:]=vTermsD3[:,k]+vu3*vQD3
    vQD3[:]*=vu[:,2]      
    fY=np.sqrt(1.5)*vMod*(1.0+vvP+vvQ)
    vPQ=1.0-degPm1*vvP-degQm1*vvQ
    D1=vu[:,0]*vPQ+vPD1+vQD1
    D2=vu[:,1]*vPQ+vPD2+vQD2
    DZ=(vu[:,2]*vPQ+vPD3+vQD3)
    DX=np.sqrt(2.0/3.0)*D1
    DY=D2/sq2-D1/sq6
    return fY,[DX,DY,DZ]
    




### function 'getCoeff': utility function
### use it to extract back the SHYqp coefficients from a '*_Err_and_Coeff.txt' file    
def getCoeff(fName):
    ff=open(fName,'r')
    line=ff.readline()
    line=ff.readline()
    deg=int(line.strip().split(':')[1])
    print('degQ = ',deg)
    while('P-Coeffs' not in line):
        line=ff.readline()
    line=ff.readline()    
    vCoeff=[]    
    while('Q-Coeffs' not in line):
        vCoeff.append(float(line.strip()))
        line=ff.readline()
    #print('P-coeff: Done')        
    while(True):
        line=ff.readline().strip()
        if(line):
            vCoeff.append(float(line))
        else:
            break
    #print('Q-coeff: Done')        
    return deg,np.array(vCoeff) 
##nQ,nP=nMonoms(degQ)
##ddMon=vPoly(degQ)



## generate a pseudo-uniform distribution of points on the unit sphere
def genSpherePoints(nPoints,moreBiax=False):
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
    nPoints=vP.shape[0]
    vP2=np.zeros((nPoints,3))
    vP2[:,:2]=vP[:,:2];vP2[:,2]=-vP[:,2]
    vP=np.concatenate([vP,vP2],axis=0)
    vP=np.concatenate([vP,-vP],axis=0)
    if(moreBiax):
        vt1=np.linspace(0,2*np.pi,600)
        vc1=np.cos(vt1);vs1=np.sin(vt1)
        vP2=np.column_stack([vc1,vs1,np.zeros(600)])
        vP3=np.column_stack([vc1,vs1,0.075*np.ones(600)])
        vP4=np.column_stack([vc1,vs1,-0.075*np.ones(600)])
        vP=np.concatenate((vP,vP2,vP3,vP4),axis=0)        
    return vP 


## generate a pseudo-uniform distribution of points on the TOP HALF of the unit sphere
def genSpherePointsTop_old(nPoints):
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
    nPoints=vP.shape[0]
    vP2=np.zeros((nPoints,3))
    vP2[:,1]=-vP[:,1];vP2[:,0]=vP[:,0];vP2[:,2]=vP[:,2]
    vP=np.concatenate([vP,vP2],axis=0)       
    return vP      

def genSpherePointsTop(nPoints):
    N1=int(0.25*nPoints)
    vt1=np.linspace(0,np.pi/2,N1+1)
    vc1=np.cos(vt1[1:]);vs1=np.sin(vt1[1:])
    vN2=np.int_(nPoints*vs1)     
    vP=np.zeros((vt1.shape[0]+np.sum(vN2),3))
    vN2[:]+=1
    vP[0,2]=1.0
    jj=1
    for kk in range(0,len(vt1)-1):
        ct1=vc1[kk];st1=vs1[kk]
        N2=vN2[kk];jN2=jj+N2
        vt2=np.linspace(0,np.pi,N2)
        vc2=np.cos(vt2[0:N2]);vs2=np.sin(vt2[0:N2])
        vP[jj:jN2,0]=st1*vc2
        vP[jj:jN2,1]=st1*vs2
        vP[jj:jN2,2]=ct1
        jj=jN2 
    nPoints=vP.shape[0]
    vP2=np.zeros((nPoints,3))
    vP2[:,1]=-vP[:,1];vP2[:,0]=vP[:,0];vP2[:,2]=vP[:,2]
    vP=np.concatenate([vP,vP2],axis=0)       
    return vP       
    
    
def genPointsCvxCheckOpt(nRandom=10**4):
    np.random.seed(99)     
    vuRandom=np.random.normal(0,1,(nRandom,3)) ##nPoints on 2D unit sphere of 3D
    vNorm=np.sqrt(np.sum(vuRandom**2,axis=1))
    vuRandom[:,0]/=vNorm;vuRandom[:,1]/=vNorm;vuRandom[:,2]/=vNorm
    vu=genSpherePoints(nPoints=300)
    vu=np.concatenate((vu[:,0:3],vuRandom),axis=0)
    #nPoints=vu.shape[0]
    return vu    
    
def HomN_GaussCheck(vf,vG,vH,Nhom):
    ncp=vf.shape[0]
    vyf=vf**(1.0/Nhom)
    vyf=vyf.reshape((ncp,1))
    tt=vyf/(Nhom*vf)
    tt=tt.reshape((ncp,1))
    G1=tt*vG[:,0].reshape((ncp,1))
    G2=tt*vG[:,1].reshape((ncp,1))
    G3=tt*vG[:,2].reshape((ncp,1))
    H11=tt*vH[:,0].reshape((ncp,1))-(Nhom-1)*(G1*G1)/vyf
    H12=tt*vH[:,1].reshape((ncp,1))-(Nhom-1)*(G1*G2)/vyf
    H13=tt*vH[:,2].reshape((ncp,1))-(Nhom-1)*(G1*G3)/vyf
    H22=tt*vH[:,3].reshape((ncp,1))-(Nhom-1)*(G2*G2)/vyf
    H23=tt*vH[:,4].reshape((ncp,1))-(Nhom-1)*(G2*G3)/vyf
    H33=tt*vH[:,5].reshape((ncp,1))-(Nhom-1)*(G3*G3)/vyf
    H11s=H22*H33-H23*H23
    H12s=H23*H13-H12*H33
    H13s=H12*H23-H22*H13
    H22s=H11*H33-H13*H13
    H23s=H12*H13-H11*H23
    H33s=H11*H22-H12*H12
    #print('here: ',H33s.shape)
    vNorm2=G1**2+G2**2+G3**2
    vNorm=np.sqrt(vNorm2)
    G1/=vNorm;G2/=vNorm;G3/=vNorm
    KG=(H11s*G1*G1+H22s*G2*G2+H33s*G3*G3+2.0*(H12s*G1*G2+H13s*G1*G3+H23s*G2*G3))/vNorm2
    return  KG
    




### Given a distribution 'vPoints' of points on the unit sphere,
### 'JensenPoints' calculates for each point 'p' three points 
### in the tangent plane at 'p' with center of mass='p' 
def JensenPoints(vPoints):
    nPoints=vPoints.shape[0]
    testPoints=np.zeros((nPoints,9))
    ee=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0]])
    sp,idx,vv,tau,tauA,tauB=0.0,0,np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3)
    norm=1.0
    c1,s1=np.cos(2.0*np.pi/3.0),np.sin(2.0*np.pi/3.0)
    c2,s2=np.cos(4.0*np.pi/3.0),np.sin(4.0*np.pi/3.0)
    for kk in range(nPoints):
        tau[:]=vPoints[kk,:]
        idx=0
        sp=np.sum(tau*ee[idx,:])
        if(sp>0.9):idx=1;sp=np.sum(tau*ee[idx,:])
        vv[:]=ee[idx,:]-sp*tau
        norm=np.sqrt(np.sum(vv*vv))
        tauA[:]=vv[:]/norm
        tauB[0]=tau[1]*tauA[2]-tau[2]*tauA[1]
        tauB[1]=tau[2]*tauA[0]-tau[0]*tauA[2]
        tauB[2]=tau[0]*tauA[1]-tau[1]*tauA[0]
        norm=np.sqrt(np.sum(tauB*tauB))
        tauB[:]/=norm
        testPoints[kk,0:3]=tauA[:]
        testPoints[kk,3:6]=c1*tauA[:]+s1*tauB[:]
        testPoints[kk,6:9]=c2*tauA[:]+s2*tauB[:]
    return testPoints


        
        
        




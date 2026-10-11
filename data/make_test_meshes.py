import numpy as np
def write(fn,V,F,UV=None):
    with open(fn,'w') as f:
        for v in V: f.write('v %.10f %.10f %.10f\n'%tuple(v))
        if UV is not None:
            for t in UV: f.write('vt %.10f %.10f\n'%tuple(t))
            for t in F: f.write('f %d/%d %d/%d %d/%d\n'%(t[0]+1,t[0]+1,t[1]+1,t[1]+1,t[2]+1,t[2]+1))
        else:
            for t in F: f.write('f %d %d %d\n'%tuple(np.array(t)+1))
# cap (disk), with a naive projection UV (x,y) as reference map
def cap(n=12,m=24,h=0.8):
    V=[[0,0,1]]
    for i in range(1,n+1):
        th=h*np.pi/2*i/n
        for j in range(m):
            ph=2*np.pi*j/m+(0.5*np.pi/m if i%2 else 0)
            V.append([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
    F=[[0,1+j,1+(j+1)%m] for j in range(m)]
    for i in range(1,n):
        a=1+(i-1)*m;b=1+i*m
        for j in range(m): F+= [[a+j,b+j,b+(j+1)%m],[a+j,b+(j+1)%m,a+(j+1)%m]]
    return np.array(V),np.array(F)
V,F=cap(); write('cap.obj',V,F,UV=V[:,:2])
# icosphere
def icosphere(sub=3):
    t=(1+5**.5)/2
    V=[[-1,t,0],[1,t,0],[-1,-t,0],[1,-t,0],[0,-1,t],[0,1,t],[0,-1,-t],[0,1,-t],[t,0,-1],[t,0,1],[-t,0,-1],[-t,0,1]]
    F=[[0,11,5],[0,5,1],[0,1,7],[0,7,10],[0,10,11],[1,5,9],[5,11,4],[11,10,2],[10,7,6],[7,1,8],[3,9,4],[3,4,2],[3,2,6],[3,6,8],[3,8,9],[4,9,5],[2,4,11],[6,2,10],[8,6,7],[9,8,1]]
    V=[np.array(v,float)/np.linalg.norm(v) for v in V]
    for _ in range(sub):
        cache={};NF=[]
        def mid(a,b):
            k=(min(a,b),max(a,b))
            if k not in cache:
                p=V[a]+V[b];V.append(p/np.linalg.norm(p));cache[k]=len(V)-1
            return cache[k]
        for a,b,c in F:
            ab,bc,ca=mid(a,b),mid(b,c),mid(c,a); NF+=[[a,ab,ca],[b,bc,ab],[c,ca,bc],[ab,bc,ca]]
        F=NF
    return np.array(V),np.array(F)
V,F=icosphere(3); V=V*np.array([1.0,0.7,0.5]); write('ellipsoid.obj',V,F)
# torus
def torus(n=40,m=20,R=1,r=0.4):
    V=[];F=[]
    for i in range(n):
        for j in range(m):
            u=2*np.pi*i/n;v=2*np.pi*j/m
            V.append([(R+r*np.cos(v))*np.cos(u),(R+r*np.cos(v))*np.sin(u),r*np.sin(v)])
    for i in range(n):
        for j in range(m):
            a=i*m+j;b=((i+1)%n)*m+j;c=((i+1)%n)*m+(j+1)%m;d=i*m+(j+1)%m
            F+=[[a,b,c],[a,c,d]]
    return np.array(V),np.array(F)
V,F=torus(); write('torus.obj',V,F)

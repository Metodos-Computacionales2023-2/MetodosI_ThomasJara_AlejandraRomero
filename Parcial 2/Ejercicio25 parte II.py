import sympy as sym
import numpy as np
x= sym.Symbol("x",Real=True)
y= sym.Symbol("y",Real=True)
n=2

def GetLaguerre(n,x):
    if n==0:
        poly = sym.Number(1)
    elif n==1:
        poly = 1-x
    else:
        poly = ((2*(n-1)+1-x)*GetLaguerre(n-1,x)-(n-1)*GetLaguerre(n-2,x))/(n)
   
    return sym.expand(poly,x)
def GetDLaguerre(n,x):
    Pn = GetLaguerre(n,x)
    return sym.diff(Pn,x,1)
def GetNewton(f,df,xn,itmax=10000,precision=1e-14):
    
    error = 1.
    it = 0
    
    while error >= precision and it < itmax:
        
        try:
            
            xn1 = xn - f(xn)/df(xn)
            
            error = np.abs(f(xn)/df(xn))
            
        except ZeroDivisionError:
            print('Zero Division')
            
        xn = xn1
        it += 1
        
    if it == itmax:
        return False
    else:
        return xn
def GetRoots(f,df,x,tolerancia = 10):
    
    Roots = np.array([])
    
    for i in x:
        
        root = GetNewton(f,df,i)

        if  type(root)!=bool:
            croot = np.round( root, tolerancia )
            
            if croot not in Roots:
                Roots = np.append(Roots, croot)
                
    Roots.sort()
    
    return Roots
def GetAllRootsGLag(n):

    xn = np.linspace(0,n+(n-1)*(n**(1/2)),100)
    
    Laguerre = []
    DLaguerre = []
    
    for i in range(n+1):
        Laguerre.append(GetLaguerre(i,x))
        DLaguerre.append(GetDLaguerre(i,x))
    
    poly = sym.lambdify([x],Laguerre[n],'numpy')
    Dpoly = sym.lambdify([x],DLaguerre[n],'numpy')
    Roots = GetRoots(poly,Dpoly,xn)

    if len(Roots) != n:
        ValueError('El número de raíces debe ser igual al n del polinomio.')
    
    return Roots
def GetWeightsGLag(n):

    Roots = GetAllRootsGLag(n)

    DLaguerre = []
    
    for i in range(n+1):
        DLaguerre.append(GetLaguerre(i+1,x))
    
    Dpoly = sym.lambdify([x],DLaguerre[n],'numpy')
    Weights = (Roots)/((n+1)**2*(Dpoly(Roots))**2)
    
    return Weights
print(GetWeightsGLag(n))
#Los resultados son los siguiente: [0.85355339 0.14644661]


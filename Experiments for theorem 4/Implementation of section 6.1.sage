
#Implementation of section 6.1



		
p=5748807539017702360253319010999088144012726967851985367500146318267957074534925160103194116945235074554875586145392358191684453970325928931959667840219964464748401277481704415790487893532578153076323436880582050234300403800137
q=3519016870771270648543366986935010687561815857359691693435895779291629268139843310188884623813310251155732262846042084327493598849030523012760650888680833888192311370027521414306561221313207362880019974195239882715761458208993
print('p=',p) 
print('\n')
print('the length of p is ',p.nbits())
print('\n')
print('q=',q)
print('\n')
print('the length of q is ',q.nbits())
print('\n')


# To compute \tilde{p}
List=p.bits()
# let nlsb be the number of least significant bits of p
nlsb=600
print('the number of known most signifant bits of p is',p.nbits()-nlsb)
print('\n')
ss=0
for i in srange(nlsb):
    ss=ss+2^i
p1=ss
ss=0
for i in srange(nlsb,p.nbits()):
    ss=ss+List[i]*2^i
p2=ss
p0=p1+p2
print('p0=',p0)
print('\n')
N=p*q
print('N=',N)
print('\n')
phi=(p-1)*(q-1)
print('phi=',phi)
print('\n')
print('the length of N is ',N.nbits())
print('\n')
# To compute \tilde{q}
q0=floor(N/p0)
print('q0=',q0)
print('\n')
u=p+q-p0-q0
print('u=',u)
print('\n')


	
r=66648350391960798602407896493517617079851395986506093699268036723928579696998993737014907386100496261948392581
print('r=',r)
print('\n')
s=1957387984051876164535698583906365043842863609584676870413702683133840359463826244243075276509507478919826621
print('s=',s)
print('\n')
e=Mod(s/r,phi).lift()
print('e=',e)
print('\n')
print('the length of e is ',e.nbits())
print('\n')
x=(e*r-s)/phi
print('x=',x)
print('\n')
b0=p0+q0-N-1
print('b0=',b0) 
print('\n')


# To verify that e*r-x*phi=s, the difference e*r-x*phi-s must be zero
print('e*r-x*phi-s=',e*r-x*phi-s)
print('\n')
# To compute the exponents: alpha,delta,epsilon,gamma

alpha= 0.99999
delta= 0.244
epsilon= 0.4
gamma= 0.24
print('\n')

print('alpha=',alpha)
print('delta=',delta)
print('epsilon=',epsilon)
print('gamma=',gamma)
print('\n')

#To check the conditions of Theorem 4


print('To check that 1-0.5*gamma-2*epsilon<delta<1-0.5*gamma-sqrt(alpha*epsilon), the answer must be true')
print('\n')
print(1-0.5*gamma-2*epsilon<delta<1-0.5*gamma-sqrt(alpha*epsilon))
print('\n')
print('1-0.5*gamma-2*epsilon=',1-0.5*gamma-2*epsilon)
print('\n')
print('delta=',delta)
print('\n')
print('1-0.5*gamma-sqrt(alpha*epsilon)=',1-0.5*gamma-sqrt(alpha*epsilon))
print('\n')

print('To check that epsilon<alpha<4*epsilon the answer must be true')
print('\n')
print(epsilon<alpha<4*epsilon)
print('\n')
print('epsilon=',epsilon)
print('\n')
print('alpha=',alpha)
print('\n')
print('4*epsilon=',4*epsilon)
print('\n')

print('To check that |s|<inf(|er|,xu), the answer must be true')
print('\n')
print(abs(s)<min(abs(e*r),abs(x*u)))
print('\n')
print('|s|=',abs(s))
print('\n')
print('inf(|er|,xu)=',min(abs(e*r),abs(x*u)))
print('\n')





# To compute the tuples (c,a,b) for (a,b) in E_c \cup F_c; for all c in the range [0,v]
	
# The list from E_c
def L1(v,kappa) :
    List =[]
    for c in srange (v+1) :
        for a in srange (c , c+1) :
            for b in srange (c,c+floor(kappa*c)+2) :
                List.extend([(c,a,b)])
    return sorted(List)
	
# The list from F_c
def L2(v,kappa) :
    List =[]
    for c in srange (v+1) :
        for a in srange (c+1 ,v+1) :
            for b in srange(c ,c+1):
                List.extend([(c,a,b)])
    return sorted(List)
	
# To compute the union from E_c and F_c 
def L(v,kappa):
    Ic=L1(v,kappa);
    Jc=L2(v,kappa);
    Ic.extend(Jc)
    return sorted(Ic)
	
	
# To define the ring Z[z_1,z_2,z_3,z_4] with z_4=z_1z_2+z_3
R.<z1,z2,z3,z4,w>=ZZ['z1,z2,z3,z4,w']

# We define an auxiliary variable $w$, which will later represent z_4 - z_3
I=(z1*z2-w)*R
# To satisfy  w=z1*z2, we define the quotient R/I
Q=R.quo(I)

#To compute the polynomial T(z1,z2,z3)
T=z1*z2+b0*z1+z3
#To compute the polynomial theta(z1,z4)
theta=z4+b0*z1
v=5;
kappa=1;
L=L(v,kappa)


# we recall that (z_1,z_2,z_3)=(−x, u, s) and z4=z_1z_2+z_3=-xu+s


print('T at the root (z1,z2,z3) modulo e must be zero')
print('\n')
print('Mod(T(z1,z2,z3),e)=',Mod(T(-x,u,s,z4,w),e))
print('\n')


print('theta at the root (z1,z4) modulo e must be zero')
print('\n')
print('Mod(theta(z1,z4),e)=',Mod(theta(-x,z2,z3,-x*u+s,w),e))
print('\n')



# To  replace z1*z2 with z4-z3, we replace w with z4-z3 in the computations

print( 'To verify that  S_{c,a,b}(z1,z2,z3,z4)=0 (mod e^v)')
for (c,a,b) in L:
    S=Q(z1^(a-c)*z2^(b-c)*z3^(v-a)*theta^c*e^(v-c)).lift().subs(w=z4-z3)	
    print([c,a,b],Mod(S(-x,u,s,-x*u+s,w),e^v))
print('\n')

print('Bounds of roots are')
Z1=floor(4*N^(alpha+delta-1))
print('Z1=',Z1)
Z2=floor(2*N^(epsilon))
print('Z2=',Z2)
Z3=floor(N^(gamma))
print('Z3=',Z3)
Z4=floor(16*N^(alpha+delta+epsilon-1))
print('Z4=',Z4)
print('\n')

print('To verify that Z1,Z2,Z3,Z4 are bounds of z1,z2,z3,z4, the answer must be true')
print('\n')
print(abs(-x)<Z1 and abs(u)<Z2 and abs(s)<Z3 and abs(-x*u+s)<Z4)
print('\n')

#The construction of the matrix of the lattice 

M=Matrix(0,len(L))  
for (c,a,b) in L:
    S=Q(z1^(a-c)*z2^(b-c)*z3^(v-a)*theta^c*e^(v-c)).lift().subs(w=z4-z3)	
    h=S(z1*Z1,z2*Z2,z3*Z3,z4*Z4,w)
    L1=vector(h[z1^(a-c)*z2^(b-c)*z3^(v-a)*z4^c] for (c,a,b) in L)
    M = M.stack(L1)
	
print('(v,kappa)=',(v,kappa))
print('The dimension of the lattice  is', M.nrows(), '\n')




#The polynomials of the rows of matrix before LLL

#To find the polynomials from the matrix
L2=[z1^(a-c)*z2^(b-c)*z3^(v-a)*z4^c for (c,a,b) in L]
r=0
LM=[]
for i in srange(len(L)):
        D=M[i,:]*transpose(matrix([L2]))
        r=D[0,0]
        LM.extend([r])
	
	
print('To verify that  S_{c,a,b}(z1,z2,z3,z4)=0 mod e^v from the matrix')
for i in srange(len(L)):
    h=LM[i]((-x)/Z1,(u)/Z2,(s)/Z3,(-x*u+s)/Z4,w)
    print([i,Mod(h,e^v)])

print('\n')	

# Performing the LLL reduction

import time
from datetime import timedelta
start_time = time.monotonic()
MR=M.LLL() 
end_time = time.monotonic()
print('\n')
print('time of LLL in seconds');
print(timedelta(seconds=end_time - start_time));
print('\n')

#To collect the polynomials h(z1,z2,z3,z4) after the LLL algorithm

#To find the polynomials from  LLL-matrix

L3=[(z1/Z1)^(a-c)*(z2/Z2)^(b-c)*(z3/Z3)^(v-a)*(z4/Z4)^c for (c,a,b) in L]
r=0
LL=[]
for i in srange(len(L)):
        D=MR[i,:]*transpose(matrix([L3]))
        r=D[0,0]
        LL.extend([r])

print('To check that polynomilas after LLL are equal to 0 mod e^v')
for i in srange(len(L)):
    S=LL[i]
    h=S(-x,u,s,-x*u+s,w)
    print([i,Mod(h,e^v)])
print('\n')	

print('To collect the polynomials such that h(x0,y0,z0,w0)=0 over the integers (Theorem 2 of Howgrave Graham)')
Lis=[]
print('The list in which h(z1,z2,z3,z4)=0 ')
for i in srange(len(L)):
    S=LL[i]
    h=S(-x,u,s,-x*u+s,w)
    if h==0:
        Lis.extend([i])
print(Lis)
print('\n')
print('verification that h(z1,z2,z3,z4)=0 ')
for i in Lis:
    S=LL[i]
    h=S(-x,u,s,-x*u+s,w)
    print([i,h])

# Performing Grobner basis computations 

print('\n')

start_time = time.monotonic()
Gg=[]
for i in Lis:
    Gg.extend([LL[i]])
I=Ideal(Gg)
E=I.groebner_basis()
z1,z2,z3,z4=var('z1 z2 z3 z4')
p1=E[0](z1,z2,z3,z4,1)
p2=E[1](z1,z2,z3,z4,1)
p3=E[2](z1,z2,z3,z4,1)
p4=E[3](z1,z2,z3,z4,1)
from sympy import solve
print('\n')
print('The list of solutions is')
print('\n')
print(solve([p1, p2, p3, p4], z1, z2, z3, z4, dict=True))
end_time = time.monotonic()
print('\n')
print('time of Grobner basis computations ');
print(timedelta(seconds=end_time - start_time));
print('\n')
print('To check that the solution z1;z2;z2;z4 is in the list')
print('\n')
print('z1=',-x)
print('z2=',u)	
print('z3=',s)	
print('z4=',-x*u+s)	
print('\n')
print('To check that the private exponent d is of magnitude N')
print('\n')
d=Mod(1/e,phi).lift();print('d=',d)
print('\n')
print('d=N^(delta0) with delta0=',numerical_approx(log(d)/log(N),digits=4))



import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt("matrice_Vandermonde.txt", dtype= float, delimiter=',')
b = np.loadtxt("cazuri.txt", dtype = float, delimiter=',')
c = np.loadtxt("decese.txt", dtype = float, delimiter=',')

print('------------Matricea A de tip Vandermonde-------------\n', A)
print('------------Vectorul nr de cazuri covid-------------\n', b)
print('------------Vectorul nr de decese covid-------------\n', c)

solx = np.linalg.pinv(np.copy(A))@np.copy(b)
soly = np.linalg.pinv(np.copy(A))@np.copy(c)
m,n  = np.shape(A)

# TORT - triangularizare ortogonala cu reflectori 
def  TORT (A):

    U = np.zeros((m,n))
    p = min(m-1,n)
    beta = np.zeros((p,1))
    for k in range(p):
        # Calcularea lui sigma
        s=0
        for i in range(k,m):
         s = s + A[i, k]**2 
        sigma = np.sign(A[k,k])*np.sqrt(s)
        
        if sigma != 0:
            # Calcularea vectorilor householder
            U[k,k] = A[k,k] +sigma
            for i in range(k+1,m):
             U[i,k] = A[i,k]

            beta[k] = sigma*U[k,k]

            #Aplicarea reflectorului Householder pe A
            for j in range(k,n):
                tau = 0
                for q in range(k,m):
                 tau = tau + U[q,k]*A[q,j]
                tau =tau/beta[k]
                
                for i in range(k,m):
                 A[i,j] = A[i, j] - tau *U[i,k]
                    
    return U, A, beta

# Algoritmul UTRIS
def Utris(U, b):
    n = len(U)
    x = np.zeros((n,1))
    for i in range(n-1,-1, -1):
        s = b[i]
        for j in range(i+1,n):
            s = s -U[i][j]*x[j]
        x[i] = s/U[i][i]
    return x

def Utris(U, c):
    n = len(U)
    y = np.zeros((n,1))
    for i in range(n-1,-1, -1):
        s = c[i]
        for j in range(i+1,n):
            s = s -U[i][j]*y[j]
        y[i] = s/U[i][i]
    return y

#Triangularizarea ortogonala a lui A
U,R,beta  = TORT(A)

# Aplicarea Reflectorilor asupra lui b
for k in range(n):
    tau=0
    for i in range(k,m):
     tau = U[i,k]*b[i] + tau 
    tau =tau/beta[k]

    for i in range(k,m):
     b[i] = b[i] - tau*U[i,k]

for k in range(n):
    tau=0
    for i in range(k,m):
     tau = U[i,k]*c[i] + tau 
    tau =tau/beta[k]

    for i in range(k,m):
     c[i] = c[i] - tau*U[i,k]

#Calcularea solutiei CMMP
x  = Utris(R[0:n,:], b[0:n])
y  = Utris(R[0:n,:], c[0:n])

print('-----------x------------\n',x)
print('-------solutia x a sistemului Ax=b---------\n', solx)

print('-----------y------------\n',y)
print('-------solutia y a sistemului Ay=c--------\n', soly)

caz_2024 = np.array([])
decese_2024 = np.array([])

print("Pentru anul 2024 avem urmatoarele predictii:")
for i in range (12):
    print("\nNumarul de cazuri noi pentru luna ", end = " ")
    if i == 0:
       print("ianuarie va fi ", end = " ")
       caz = x[0] + x[1]*7 + x[2]*49
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna ianurie va fi ", end = " ")
       deces = y[0] + y[1]*7 + y[2]*49
       print(deces)
       decese_2024 = np.append(decese_2024,deces)
    if i == 1:
       print("februarie va fi ",end = " ")
       caz = x[0] + x[1]*8 + x[2]*64
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna februarie va fi ",end = " ")
       deces = y[0] + y[1]*8 + y[2]*64
       print(deces)
       decese_2024 = np.append(decese_2024,deces)
    if i == 2:
       print("martie va fi ",end = " ")
       caz = x[0] + x[1]*9 + x[2]*81
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna martie va fi ",end = " ")
       deces = y[0] + y[1]*9 + y[2]*81
       print(deces)
       decese_2024 = np.append(decese_2024,deces)
    if i == 3:
       print("aprilie va fi ",end = " ")
       caz = x[0] + x[1]*10 + x[2]*100
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna aprilie va fi ",end = " ")
       deces = y[0] + y[1]*10 + y[2]*100
       print(deces)
       decese_2024 = np.append(decese_2024,deces)
    if i == 4:
       print("mai va fi ",end = " ")
       caz = x[0] + x[1]*11 + x[2]*121
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna mai va fi ",end = " ")
       deces = y[0] + y[1]*11 + y[2]*121
       print(deces)
       decese_2024 = np.append(decese_2024,deces)
    if i == 5:
       print("iunie va fi ",end = " ")
       caz = x[0] + x[1]*12 + x[2]*144
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna iunie va fi ",end = " ")
       deces = y[0] + y[1]*12 + y[2]*144
       print(deces)
       decese_2024 = np.append(decese_2024,deces)
    if i == 6:
       print("iulie va fi ",end = " ")
       caz = x[0] + x[1]*13 + x[2]*169
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna iulie va fi ",end = " ")
       deces = y[0] + y[1]*13 + y[2]*169
       print(deces)
       decese_2024 = np.append(decese_2024,deces)
    if i == 7:
       print("august va fi ",end = " ")
       caz = x[0] + x[1]*14 + x[2]*196
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna august va fi ",end = " ")
       deces = y[0] + y[1]*14 + y[2]*196
       print(deces)
       decese_2024 = np.append(decese_2024,deces)
    if i == 8:
       print("septembrie va fi ",end = " ")
       caz = x[0] + x[1]*15 + x[2]*225
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna septembrie va fi ",end = " ")
       deces = y[0] + y[1]*15 + y[2]*225
       print(deces)
       decese_2024 = np.append(decese_2024,deces)
    if i == 9:
       print("octombrie va fi ",end = " ")
       caz = x[0] + x[1]*16 + x[2]*256
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna octombrie va fi ",end = " ")
       deces = y[0] + y[1]*16 + y[2]*256
       print(deces)
       decese_2024 = np.append(decese_2024,deces)
    if i == 10:
       print("noiembrie va fi ",end = " ")
       caz = x[0] + x[1]*17 + x[2]*289
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna noiembrie va fi ",end = " ")
       deces = y[0] + y[1]*17 + y[2]*289
       print(deces)
       decese_2024 = np.append(decese_2024,deces)
    if i == 11:
       print("decembrie va fi ",end = " ")
       caz = x[0] + x[1]*18 + x[2]*324
       print(caz)
       caz_2024 = np.append(caz_2024, caz)
       print("\nNumarul de decese pentru luna decembrie va fi ",end = " ")
       deces = y[0] + y[1]*18 + y[2]*324
       print(deces)
       decese_2024 = np.append(decese_2024,deces)


lunile = ["Ianuarie", "Februarie", "Martie", "Aprilie", "Mai", "Iunie", "Iulie", "August", "Septembrie", "Octombrie", "Noiembrie", "Decembrie"]
caz_2023 = np.array([3.325066, 3.340342, 3.374825, 3.393902, 3.404197, 3.407995, 3.409997, 3.417177, 3.468372, 3.497226, 3.505080, 3.512369])
decese_2023 = np.array([0.067576, 0.067704, 0.067917, 0.068089, 0.068184, 0.068240, 0.068253, 0.068270, 0.068372, 0.068525, 0.068653, 0.068728])

# Creăm o serie de valori numerice pentru axa X
luni_numerice = range(len(lunile))

plt.scatter(luni_numerice, caz_2023, color='blue', label='Date reale', marker='o')
plt.xticks(luni_numerice, lunile, rotation=25)  # Setăm etichetele axei X cu numele lunilor
plt.xlabel('Lunile anului 2023')
plt.ylabel('Cazuri COVID')
plt.title('Grafic COVID 2023')
plt.legend()
plt.show()

plt.scatter(luni_numerice, decese_2023, color='pink', label='Date reale', marker='o')
plt.xticks(luni_numerice, lunile, rotation=45) 
plt.xlabel('Lunile anului 2023')
plt.ylabel('Numarul de decese')
plt.title('Grafic COVID 2023')
plt.legend()
plt.show()

plt.scatter(luni_numerice, caz_2024, color='red', label='Date prezise', marker='o')
plt.xticks(luni_numerice, lunile, rotation=45) 
plt.xlabel('Lunile anului 2024')
plt.ylabel('Cazuri covid')
plt.title('Predictia COVID pentru anul 2024 ')
plt.legend()
plt.show()

plt.scatter(luni_numerice, decese_2024, color='green', label='Date prezise', marker='o')
plt.xticks(luni_numerice, lunile, rotation=45) 
plt.xlabel('Lunile anului 2024')
plt.ylabel('Numarul de decese')
plt.title('Predictia COVID pentru anul 2024')
plt.legend()
plt.show()
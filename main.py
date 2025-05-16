import numpy as np
import matplotlib.pyplot as plt

#Le code risque d'etre peu compréhensible sans avoir lu les slides


#Variables globales

#Constantes ******************************
alpha = 0.1 
beta = 0.02 
gamma = 0.03 
delta = 0.01 
K = 100 
#Constantes ******************************

#Autres variables ************************
max_iter = 10000
eps = 1e-6
Tf = 450 #temps final 


    #Conditions initiales
CA0 = 80
CT0 = 20
CS0 = 10

C0 = [CA0, CT0, CS0]
#Plus tard, le tableau de tableau C contiendra à chaque colonne n le tableau Cn
t0 = 0 #temps initial
h = 0.1 #pas entre chaque temps
#Autres variables ************************

#fonctions********************************

def f(Cn): 
    """
        fonction f définie dans la présentation.
        Args :
            Cn : numpy array ou python array (représente l'état du système au temps tn)
        Returns :
            numpy array représente f(Cn)
            
    """
    y1 = Cn[1] * (-alpha * (1 - (Cn[1])/K) + beta) + delta * Cn[2]
    y2 = Cn[1] * (alpha * (1 - (Cn[1])/K) - beta - delta - gamma)
    y3 = Cn[1] * (gamma + delta) - delta * Cn[2]
    return np.array([y1,y2,y3])


def F1(Cnk_suiv, Cnk_prec):
    """
        fonction F1 qui est l'itération de Euler Implicite
        Args:
            Cnk_suiv : numpy array ou python array (terme actuel de la suite approchant Cn, c'est Cn,k)
            Cnk_prec : numpy array ou python array (terme précédent de la suite approchant Cn)
        Returns:
            numpy array : représente le terme suivant de la suite approchant Cn, c'est Cn,k+1
    """
    Cnk_prec = np.array(Cnk_prec)
    Cnk_suiv = np.array(Cnk_suiv)
    return Cnk_prec + h*f(Cnk_suiv) #Cn,k+1 = Cn + hf(Cn,k)


def F2(Cnk_suiv, Cnk_prec):
    """
        fonction F2, même que F1 mais avec la méthode des trapèzes
        Args:
            Cnk_suiv : numpy array ou python array
            Cnk_prec : numpy array ou python array
        Returns:
            numpy array
    """
    Cnk_prec = np.array(Cnk_prec)
    Cnk_suiv = np.array(Cnk_suiv)
    return Cnk_prec + (h/2)*(f(Cnk_suiv) + f(Cnk_prec))

def f_newton(Cnk_suiv, Cnk_prec):
    """
        fonciton f_newton qui implémente celle décrite dans le rapport, on cherche à l'annuler
        Args:
            Cnk_suiv : numpy array ou python array
            Cnk_prec : numpy array ou python array
        Returns:
            numpy array 
    """
    return F2(Cnk_suiv, Cnk_prec) - Cnk_suiv

def df_newton(Cnk):
    """
        fonction df_newton qui représente la jacobienne de f_newton évaluée en le veccteur Cnk
        Args:
            Cnk : numpy array ou python array
        Returns:
            numpy array : jacobienne de f_newton évaluée en Cnk
    """
    A = np.zeros((3,3))
    A[:,1] = [-alpha + (2*Cnk[1])/K + beta, alpha - (2*Cnk[1])/K - beta - delta - gamma, gamma + delta]
    A[:,2] = [delta, 0, -delta]
    return (h/2) * A - np.eye(3)

#fonctions********************************



#methodes numeriques**********************

def point_fixe(X0, _F, _eps=eps, _max_iter = max_iter):
    """
        fonction point_fixe : résout X = _F(X) par la méthode du point fixe
        Args:
            X0 : numpy array ou python array (premier terme de la suite Cnk, (représente Cn-1))
            _F : fonction de point fixe (on pourra utiliser F1 ou F2 définits plus haut)
            _eps : float (tolérance d'erreur)
            _max_iter : int (nombre maximum d'itérations de la méthode)
        Returns:
            numpy array : c'est le terme suivant (donc Cn+1)
    """
    X0 = np.array(X0)
    Xk = X0 #Xk = Cn,k
    Xk_1 = X0 #Xk_1 représente Xk+1 soit Cn,k+1
    for i in range(_max_iter):
        Xk_1 = _F(Xk, X0)
        if(np.linalg.norm(Xk - Xk_1) < eps):
            break
        Xk = Xk_1
    return Xk_1

def newton(X0, f, df, _eps=eps, _max_iter= max_iter):
    """
        fonction newton : résout f(X) = 0 par la méthode de newton
        Args:
            X0 : numpy array ou python array (premier terme de la suite de newton Xnk, (représente Cn-1))
            f : fonction de newton choisie 
            _eps : float (tolérance d'erreur)
            _max_iter : int (nombre maximum d'itérations de la méthode)
        Returns:
            numpy array : c'est le point où f s'annule (ici c'est Cn+1)
    """
    X0 = np.array(X0)
    Xk = X0
    for i in range(_max_iter):
        if(np.linalg.norm(f(Xk, X0)) < _eps):
            break
        Xk = Xk - np.linalg.inv(df(Xk)) @ f(Xk, X0)
    return Xk

def eulerImplicite(_C0):
    """
        applique la méthode d'euler implicite pour résoudre le probleme, et utilise le point fixe pour approcher chaque terme Cn
        Args:
            _C0 : numpy array ou python array (condition initiale de l'EDO)
        Returns :
            C : numpy array contenant à la colonne n le vecteur Cn
            T : python array qui contient tous les temps tn
    """
    t = t0
    C = np.zeros((3,1))
    T = [t0]
    C[:,0] = _C0
    k = 1
    while(t < Tf):
        t = t+h
        C = np.append(C, np.transpose([point_fixe(C[:,k-1],F1)]), axis=1) #ajoute Cn au tableau C tel que à la colonne i on ait Ci
        T.append(t)
        k = k+1
    return C, T


def eulerExplicite(_C0):
    """
        applique la méthode d'euler explicite pour résoudre le probleme
        Args:
            _C0 : numpy array ou python array (condition initiale de l'EDO)
        Returns :
            C : numpy array contenant à la colonne n le vecteur Cn
            T : python array qui contient tous les temps tn
    """
    t = t0
    C = np.zeros((3,1))
    T = [t0]
    C[:,0] = _C0
    k = 1
    while(t < Tf):
        t = t+h
        C = np.append(C, np.transpose( [C[:,k-1] + h*f(C[:,k-1])] ), axis = 1) #ajoute Cn au tableau C tel que à la colonne i on ait Ci
        T.append(t)
        k = k+1
    return C, T

def trapeze_newton(_C0):
    """
        applique la méthode de newton et des trapèzes combinés (comme décrit dans les slides) pour résoudre le problème
        Args:
            _C0 : numpy array ou python array (condition initiale de l'EDO)
        Returns :
            C : numpy array contenant à la colonne n le vecteur Cn
            T : python array qui contient tous les temps tn
    """
    t = t0
    C = np.zeros((3,1))
    T = [t0]
    C[:,0] = _C0
    k = 1
    while(t < Tf):
        t = t+h
        C = np.append(C, np.transpose( [newton(C[:,k-1],f_newton, df_newton)] ), axis = 1) #ajoute Cn au tableau C tel que à la colonne i on ait Ci
        T.append(t)
        k = k+1
    return C, T

#methodes numeriques**********************

C1, T1 = eulerImplicite(C0)
C2, T2 = eulerExplicite(C0)
C3, T3 = trapeze_newton(C0)

plt.subplot(3,1,1)
plt.plot(T1, C1[0], label='CA(t)')
plt.plot(T1, C1[1], label='CT(t)')
plt.plot(T1, C1[2], label='Cs(t)')
plt.legend()
plt.title("Point fixe à partir d'Euler Implicite")

plt.subplot(3,1,2)
plt.plot(T2, C2[0], label='CA(t)')
plt.plot(T2, C2[1], label='CT(t)')
plt.plot(T2, C2[2], label='CS(t)')
plt.legend()
plt.title("Euler explicite")

plt.subplot(3,1,3)
plt.plot(T3, C3[0], label='CA(t)')
plt.plot(T3, C3[1], label='CT(t)')
plt.plot(T3, C3[2], label='CS(t)')
plt.legend()
plt.xlabel("t")
plt.title("trapèze newton")

plt.show()

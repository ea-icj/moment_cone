from cone.root import *
from cone.weight import *
from cone.inequality import *

import itertools

class Monomial:
    """
    A monomial is a list of roots that represent de product of the corresponding variables in the ring QV.
    We allow toral roots [k,i,i] that correspondant to 1's in the monomial.
    A monomial has a degree.
    """
    def __init__(self, roots):
        """
        Initialize a Monomial object.

        :param roots: List of Root objects.
        """
        # Validate that roots is a list; otherwise, raise an error
        if not isinstance(roots, list):
            raise ValueError("The parameter 'roots' must be a list of Root objects.")
        
        # Validate that all elements are instances of Root
        if not all(isinstance(beta, Root) for beta in roots):
            raise ValueError("All elements of 'roots' must be instances of the Root class.")
        
        # Assign the roots
        self.roots = roots
        
        # Calculate the default degree, excluding toral roots [k, i, i]
        self._degree = sum(1 for beta in self.roots if beta.i != beta.j)

    #def __init__(self, roots):
        """
        Initialize a Monomial object.

        :param roots: List of Root objects.
        """
     #   self.roots = roots if isinstance(roots, list) else []
        # Calculate the default degree as the length of the list of roots
      #  self._degree = sum([1 for beta in self.roots if beta.i!=beta.j])

    @property
    def degree(self):
        """
        Getter for degree.
        """
        return self._degree

    @degree.setter
    def degree(self, value):
        """
        Setter for degree. Allows modification with `monomial.degree = new_value`.
        """
        if not isinstance(value, int):
            raise ValueError("Degree must be an integer.")
        self._degree = value

    def __repr__(self):
        """
        String representation of the Monomial.
        """
        roots_repr = ", ".join(str(root) for root in self.roots)
        return f"Monomial(roots=[{roots_repr}], degree={self._degree})"

def Is_Max_un(l : list[int]) -> bool:
    """
    Check if the maximum of a list of integers is at most 1 or not.
    """
   for x in l :
       if x>1 :
           return(False)
   return(True)          

def List_var_1(l : list[Monomial]) -> list[Root]:
    """
    Given a list of monomials, return the list of roots that give linear equalities. 
    """
    ll=[]
    for m in l :
        if m.degree==1 : 
           for beta in m.roots:
               if beta.i!=beta.j and beta not in ll :
                   ll.append(beta)
                   
    return(ll)

def Is_Lin_Triangular(equa) -> bool: #equa is dictionnary (index of a Weight) -> list of Monomials
    """
    Check if a given list of equalities is Linear Trinagular or not.
    Each equality is actually given by its list of monomials (without coefficients)
    """
    if len(equa.keys())==0 : return(True)
    List_eq_lin=[]
    list_var_deg_1=[]
    for id_p in equa.keys():
        if Is_Max_un([m.degree for m in equa[id_p]]) :
            List_eq_lin.append(id_p)
            for beta in List_var_1(equa[id_p]):
                if beta not in list_var_deg_1 :
                    list_var_deg_1.append(beta)
    if len(List_eq_lin)==0 :
        return(False) 
    if len(List_eq_lin)>len(list_var_deg_1) : return('Trop ineg lin')   
    if len(List_eq_lin)!=len(list_var_deg_1) :
        return(False)
       
    # We put the linear variables at 0
    newequa={}
    for id_p in equa.keys():
        if id_p not in List_eq_lin :
            lmon=[]
            for x in equa[id_p]:
                if not(set(x.roots) & set(list_var_deg_1)) :
                    #nx=x.copy()
                    lmon.append(x)               
            newequa[id_p]=lmon
    return(Is_Lin_Triangular(newequa))        

def Fiber_LT(ineq : Inequality) -> bool:
    """
    Input : an inequality. 
    Return True iff the equations of the fibers of pi are linear triangular.
    """
    tau=ineq.tau
    d=ineq.tau.d
    # Inversions of U(w)xT with a dictionnary pos --> list of roots
    inv_w_dic={}
    for k in range(len(d)):
        inv_w_dic[k]=[Root(k,i,i) for i in range(d[k])]                
    for beta in ineq.inversions:
        inv_w_dic[beta.k].append(beta)    

    # Weights with positive scalar with ineq.tau    
    Poids_pos=[]
    for ll in list(tau.positive_weights.values()):
        Poids_pos+=ll

    # Weights with non positive scalar with ineq.tau       
    Poids_neg0=[]
    for ll in list(tau.non_positive_weights.values()):
        Poids_neg0+=ll

    #Initialization of equa as the dictionnary positive weight --> empty list    
    equa={}
    for p in Poids_pos:
        equa[p.index_in(d)]=[]
        
    
    for betas in itertools.product(*inv_w_dic.values()):
        weight_new=Weight([beta.i for beta in betas])
        weight_old=Weight([beta.j for beta in betas])
        if (weight_new in Poids_pos) and (weight_old in Poids_neg0) :
            equa[weight_new.index_in(d)].append(Monomial(list(betas)))
    return(Is_Lin_Triangular(equa))
       

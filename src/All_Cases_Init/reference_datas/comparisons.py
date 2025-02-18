from .Klyachko import Klyachko_3_7, Klyachko_3_8, Klyachko_4_8
from .Vergne_Walter import Vergne_Walter_444
from ..tau import *

#TODO: remove if not used
def Klyachkoineqs_to_ineqs(lK):
    """ converts a list of flattened inequalities (like the ones in Klyachko.py) to our format Inequality
    """
    n=len(lK[0])
    list_ineq=[]
    for v in lK:
        list_ineq.append(Inequality.from_tau( Tau((v,), LinGroup([n]))))
    return list_ineq

#TODO check unitary tests
def compare(list1,list2, comment0="elements", comment1="1", comment2="2"):
    """ Compares two lists. returns elements from the first one that don't belong to the second one and conversely.
    elements have to be hashable (e.g. tuples) so that set(list1) works
    optional arguments comments are here to caracterise the nature of elements (comment0) and particularities of the compared lists (cf compareK_ineq and compareVW_ineq below)
    >>>compare([(1,3,6),(2,4,9)],[(1,3,6),(12,43,24),(1,3,6)])
    2 elements from list 1  vs  2 elements from list 2
    1 elements in both lists
    1 elements appearing only in list 1 {(2, 4, 9)}
    1 elements appearing only in list 2 {(12, 43, 24)}
    [{(2, 4, 9)}, {(12, 43, 24)}]
    """
    set1=set(list1)
    set2=set(list2)
    print(len(list(set1)), comment0, "from list", comment1, " vs ", len(list(set2)), comment0, "from list", comment2)

    #if not((k,n) in [(3,7),(3,8),(4,8)]:
    #   print "case not supported by Klyachko"
    inter=set1.intersection(set2)
    print(len(list(inter)),"elements in both lists")
    only=[set1-set2, set2-set1]
    print(len(only[0]), comment0, "appearing only in list",comment1, only[0])
    print(len(only[1]), comment0, "appearing only in list",comment2, only[1])
    return only

def compareK_ineq(Klyach,list_ineq):
    """
    Klyach is a list of flattened inequalities (like the ones in Klyachko.py)
    list_ineq is a list of computed object of type Inequalities.
    the function compares the underlying inequalities (fermionic case) 
    """
    opp_Klyach=[tuple([-x for x in v]) for v in Klyach]
    list2=[ineq.wtau.components[0] for ineq in list_ineq]
    return compare(opp_Klyach,list2,"inequalities", "of reference (Klyachko)"," computed")
    
def compareVW_ineq(VW,list_ineq):
    """
    VW is a list of flattened inequalities (like the ones in Vergne_Walter.py)
    list_ineq is a list of computed object of type Inequalities.
    the function compares the underlying inequalities (kron case)"""
    G=list_ineq[0].tau.G
    opp_VW_tau=[Tau.from_flatten(tuple([-x for x in v]),G) for v in VW]
    opp_VW=[tau.sort_mod_sym_dim.flattened for tau in opp_VW_tau]
    list2=[ineq.wtau.sl_representative.sort_mod_sym_dim.flattened for ineq in list_ineq]
    return compare(opp_VW,list2,"inequalities", "of reference (Vergne Walter)"," computed")
    
def compare_to_reference(list_ineq,V):
    """
    list_ineq is a list of Inequalities computed for a certain representation V.
    If exists, it will be compared to a reference list of inequalities (currently only the cases of Klyachko.py for fermions and Vergne_Walter.py for kronecker)
    Other references can be added later
    """
    if V.type=="fermion":
        n=V.G[0]
        k=V.nb_part
        if (k,n)==(3,7):
            reference=Klyachko_3_7
        elif (k,n)==(3,8):
            reference=Klyachko_3_8
        elif (k,n)==(4,8):
            reference=Klyachko_4_8
        else:
            print("no reference for", V,"included")
            return None
        return compareK_ineq(reference,list_ineq)
    elif V.type=="kron":
        if tuple(V.G)==(4,4,4,1):
            return compareVW_ineq(Vergne_Walter_444,list_ineq)
        else:
            print("no reference for", V,"included")
            return None
    else:
        print("no reference for", V,"included")
        return None


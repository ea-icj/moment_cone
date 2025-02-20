
from .ineq_Python_4_4_4_1_kron_Vergne_Walter import *
from .ineq_Python_7_fermion_3_Klyachko import *
from .ineq_Python_8_fermion_3_Klyachko import *
from .ineq_Python_8_fermion_4_Klyachko import *

from ..tau import *


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
    print(len(only[0]), comment0, "appearing only in list",comment1)
    print(len(only[1]), comment0, "appearing only in list",comment2)
    return only


def compare_to_reference(list_ineq,V):
    """
    list_ineq is a list of Inequalities computed for a certain representation V.
    If exists, it will be compared to a reference list of inequalities (currently only the cases of Klyachko.py for fermions and Vergne_Walter.py for kronecker)
    Other references can be added later
    """
    if isinstance(V,FermionRepresentation):
        n=V.G[0]
        k=V.particle_cnt
        if (k,n)==(3,7):
            reference=inequalities_K37
        elif (k,n)==(3,8):
            reference=inequalities_K38
        elif (k,n)==(4,8):
            reference=inequalities_K48
        else:
            print("no reference for", V,"included")
            return None
        return compare(list_ineq,reference, "inequalities", " computed","of reference (Klyachko)")
    elif isinstance(V,KroneckerRepresentation):
        reference_sym=[Inequality.from_tau(ineq.wtau.sort_mod_sym_dim) for ineq in inequalities_VW444]
        list_ineq_sym= [Inequality.from_tau(ineq.wtau.sl_representative.sort_mod_sym_dim) for ineq in list_ineq]
        if tuple(V.G)==(4,4,4,1):
            return compare(list_ineq_sym, reference_sym,"inequalities (up to S3-sym)"," computed","of reference (Vergne-Walter)")
        else:
            print("no reference for", V,"included")
            return None
    else:
        print("no reference for", V,"included")
        return None





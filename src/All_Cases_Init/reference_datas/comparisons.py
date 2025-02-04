from reference_datas.Klyachko import Klyachko_3_7, Klyachko_3_8, Klyachko_4_8
from reference_datas.Vergne_Walter import Vergne_Walter_444

def Klyachkoineqs_to_ineqs(lK):
    n=len(lK[0])
    list_ineq=[]
    for v in lK:
        list_ineq.append(Inequality.from_tau( Tau((v,), LinGroup([n]))))
    return list_ineq

def compare(list1,list2, comment1="1", comment2="2"):
    print(len(list1), "inequalities from list", comment1, " vs ", len(list2), "inequalities from list", comment2)
    set1=set(list1)
    set2=set(list2)
    #if not((k,n) in [(3,7),(3,8),(4,8)]:
    #   print "case not supported by Klyachko"
    inter=set1.intersection(set2)
    print(len(list(inter)),"elements in both lists")
    only=[set1-set2, set2-set1]
    print("inequalities appearing only in list",comment1, only[0])
    print("inequalities appearing only in list",comment2, only[1])
    return only

def compareK_ineq(Klyach,list_ineq):
    opp_Klyach=[tuple([-x for x in v]) for v in Klyach]
    list2=[ineq.wtau.components[0] for ineq in list_ineq]
    return compare(opp_Klyach,list2,"of reference (Klyachko)"," computed")
    
def compareVW_ineq(VW,list_ineq):
    G=list_ineq[0].tau.G
    opp_VW_tau=[Tau.from_flatten(tuple([-x for x in v]),G) for v in VW]
    opp_VW=[tau.sort_mod_sym_dim.flattened for tau in opp_VW_tau]
    list2=[ineq.wtau.sl_representative.sort_mod_sym_dim.flattened for ineq in list_ineq]
    return compare(opp_VW,list2,"of reference (Vergne Walter)"," computed")
    
def compare_to_reference(list_ineq,V):
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


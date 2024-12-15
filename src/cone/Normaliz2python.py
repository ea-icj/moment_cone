from .typing import *
from re import split as re_split
from .tau import Tau
from .inequality import *
from .utils import dictionary_list_lengths

def convert_lines_Nout2pyth(lines:Sequence[str],d:"Dimension") -> Sequence[Tau]:
    """ lines is a list of lines read from a Normaliz output file, e.g. via
    fichier = open("/home/bm29130h/Documents/Recherche/Ressources_autres/GDT/Machine Learning/calculs Kron/2 oct/ineq_Normaliz"+st+".out","r")
    lines = fichier.readlines()
    fichier.close() 
    
    They are converted to 1 parameter compatible with our implementation. 
    """
    res=[]
    concerned_line=False
    for l in lines[16:]:
        if l=='\n':
            concerned_line=False
        if concerned_line:
            mp=[]
            m=re_split(r" ",l)
            m[-1]=m[-1][:-1]
            mp=[]
            for mm in m:
                if mm!='':
                   mp.append(int(mm))
            res.append(Tau.from_flatten(mp,d))
        if l[-21:]=="support hyperplanes:\n":
            concerned_line=True
    return(res)

def convert_file_Nout2pyth(rep_path:str, d:"Dimension") -> Sequence[Tau]:
    """rep_path is the path to the repository containing a Normaliz output, d is the dimension under concern (typically, a given repository can contain many Normaliz outputs for differents dimension d. e.g. path="/home/bm29130h/Documents/Recherche/Ressources_autres/GDT/Machine Learning/calculs Kron/2 oct/"
    
    the function then applies convert_lines_Nout2pyth
    """
    st=""
    for di in d:
       st+="-"+str(di)
    fichier = open(rep_path+"ineq_Normaliz"+st+".out","r")
    lines = fichier.readlines()
    fichier.close() 
    return(convert_lines_Nout2pyth(lines,d))

#path="/home/bm29130h/Documents/Recherche/Ressources_autres/GDT/Machine Learning/calculs Kron/2 oct/"
#reference=convert_file_Nout2pyth(path,d)

    
def compare_tau_candidates_reference_mod_sym_dim(Candidates:Sequence[Inequality],Reference:Sequence[Inequality])->dict[str,Sequence[Inequality]]:
    l1=[ineq.wtau.end0_representative.opposite for ineq in unique_modulo_symmetry_list_of_ineq(Candidates)]
    l2=[ineq.wtau for ineq in unique_modulo_symmetry_list_of_ineq(Reference)]
    res = {key: [] for key in ["candidates_only","reference_only","both"]}
    for wtau in list(set(l1+l2)):
        wtau_opp=wtau.opposite
        if wtau in l1:
            if wtau in l2: 
                res["both"].append(Inequality.from_tau(wtau_opp))
            else:
                res["candidates_only"].append(Inequality.from_tau(wtau_opp))
        else: 
            res["reference_only"].append(Inequality.from_tau(wtau_opp))
    print(dictionary_list_lengths(res))
    #print(len(res["cand_only"]), "irredundant inequalities found in given list")
    #print(len(res["both"]), "redundant inequalities found in given list")
    return(res)



#ineq_normaliz=convertre(lines)

#ineq_normaliz_SL=[tau2SL(ineq,d0) for ineq in ineq_normaliz] 

#ineq








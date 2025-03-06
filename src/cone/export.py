""" Python to normaliz and Python to latex """

__all__ = (
    "ExportFormat",
    "export_normaliz",
    "export_latex",
    "export_python",
    "export_many",
)

from .typing import *
from .tau import Tau
from .inequality import Inequality
from .representation import *

# Available export formats
ExportFormat = Literal[
    "Normaliz",
    "LaTeX",
    "Python",
    "Terminal",
    "None",
]


#######################################################################
#I. Normaliz export
#######################################################################

def Liste_to_Normaliz_string(liste: list[list[int]], sgn: int = 1) -> str: #sgn=1 ou -1 allows to change the sign, exchange >=0 and <=0
    """ converts a list of list of numbers to a string with
    *newline characters at the end of each sublist
    *a space between each number of a given sublist
    """
    chaine=''
    for l in liste:
        for x in l[:-1]: # no space at the end
            chaine+=str(sgn*x)+' '
        chaine+=str(sgn*l[-1])+'\n'
    return chaine
    

def info_from_GV(V: Representation) -> str:
    """ V a representation of a LinGroup G. It returns a nomalized string encoding the main information on G and V. To be used in the names of the input files given to Normaliz"""
    G = V.G
    info = '-'.join(str(d) for d in G)
    info += "_" + type(V).__name__
    if isinstance(V, ParticleRepresentation):
        info += "_" + str(V.particle_cnt)
    return info

def export_normaliz(
        V: Representation,
        inequations: Sequence[Inequality | list[int]],
        r: Optional[int] = None,
        #equations=[],
        extra_info: str = "",
        add_dominance: str = "",
        add_equations: str = ""
        ) -> None:
    """ V is a representation of a LinGroup G, 
    r is the rank of G aka the dimension of the ambient space of the equations (except possible customization, if our inequations contain a second member)
    inequations is a list of inequations (either of Inequality type or a list of coefficients)
    equations is an optional list of equations
    extra_info is a string that one may add to the standard name of the output file, e.g. to indicate the origin of our list of inequations 
    add_dominance can take 2 non-trivial arguments: "all" and "sym". "all" argument adds all the inequalities expressing dominance. "sym" arguments adds the same inequalities
    add_equations can take 2 non-trivial arguments: "all" and "sym". if V is KronerckerRepresentation, it adds the equations determining the subspace in which lie the equations
    """
    if r is None:
        r = V.G.rank
    
    True_inequations: list[list[int]] = []
    for ineq in inequations:
        if isinstance(ineq, Inequality):
            True_inequations.append(list(ineq.wtau.flattened))
        else:
            True_inequations.append(ineq)

    if add_dominance in ("all", "sym"):
        sym = add_dominance != "all"
        True_inequations += [list(ineq.wtau.flattened) for ineq in Inequality.dominance(V,sym)]

    info=info_from_GV(V)+extra_info
    fileO = open('ineq_Normaliz-'+info+'.in','w')
    fileO.write('amb_space '+str(r)+'\n\n')
    new_equations=[] #copy of the list
    if add_equations in ["all","sym"] and isinstance(V, KroneckerRepresentation):
        for k,dk in enumerate(V.G[:-1]):
            if add_equations=="all" or dk!=V.G[k+1]:
                new_equations.append(sum(V.G[:k])*[0]+dk*[1]+(sum(V.G[k+1:])-1)*[0]+[-1])
    if len(new_equations)!=0:
        fileO.write('equations '+str(len(new_equations))+'\n\n')
        List_Eq=Liste_to_Normaliz_string(new_equations)
        fileO.write(List_Eq)
    fileO.write('\n'+'inequalities '+str(len(True_inequations))+'\n\n')
    fileO.write(Liste_to_Normaliz_string(True_inequations,-1)) #Our conventions so far work with inequalities of the form \sum a_i\lambda_i<=0 whil Normaliz standard is ">=0"
    fileO.close()


#######################################################################
#II. Latex export
#######################################################################

def Latex_string_of_tau(
        tau: Tau,
        lambda_notation: bool = False,
        sgn: int = 1
        ) -> str:
    chaine=''
    if not lambda_notation:
        for taui in tau.components:
            if len(taui)>1:
                chaine+='('
            for x in taui:
                chaine+=str(x*sgn)+' ,'
            chaine=chaine[:-2]
            if len(taui)>1:
                chaine+=') \\;'
    else:
        started=False #stores whether we have already met a non-zero coefficient 
        for i,x in enumerate(tau.flattened):
            y=sgn*x
            if y!=0:
                if y>0:
                    if started: #(no sign + for first coeeficient)
                        chaine+=' + '
                if not(y in {1,-1}):
                    chaine+=str(y)
                elif y==-1:
                    chaine+='-'
                chaine+='\\lambda_{'+str(i+1)+'}'
                started=True
        chaine+='\\geq 0'
    return chaine 

def Latex_string_of_cluster_dom1PS(
        inequations: Sequence[Inequality],
        lambda_notation: bool = False,
        sgn: int = 1
        ) -> str: #sgn=1 ou -1 allows to change the sign, exchange >=0 and <=0
    """ converts a list of Inequalities (associated to a given dominant 1PS  taudom) to a string describing part of a latex tabular
    """
    n=len(inequations)
    chaine='\\multirow{'+str(n)+'}{*}{'
    chaine+=Latex_string_of_tau(inequations[0].tau, False, sgn)+' } ' 
    for ineq in inequations:
        chaine+=' & '+Latex_string_of_tau(ineq.wtau, lambda_notation, sgn)+' & '
        s=len(ineq.w)
        if s>1: #kron type
           s-=1
        chaine+=Latex_string_of_tau(Tau([tuple(wi) for wi in ineq.w[:s]]))
        chaine+='\\\\ \n \\cline{2-3} \n'
    chaine=chaine[:-15]
    chaine+='\\hline'
    return(chaine)
    
def group_by_dom1PS(inequations: Sequence[Inequality]) -> list[list[Inequality]]:
    """ from a list of Inequalities: groups them by same dominant 1 parameter subgroup
    """
    remaining_indices=[i for i in range(len(inequations))]
    grouped_ineqs: list[list[Inequality]] = []
    while remaining_indices != []:
        j=0
        taudom=inequations[remaining_indices[0]].tau
        taudom_ineqs: list[Inequality] = []
        while j<len(remaining_indices):
            if inequations[remaining_indices[j]].tau==taudom:
                taudom_ineqs.append(inequations[remaining_indices[j]])
                del(remaining_indices[j])
            else:
                j+=1
        grouped_ineqs.append(taudom_ineqs)
    return grouped_ineqs


#TODO make a caption

def export_latex(
        V: Representation,
        inequations: Sequence[Inequality],
        sgn: int = 1,
        extra_info: str = ""
        ) -> None: #sgn=1 ou -1 allows to change the sign, exchange >=0 and <=0
    """ converts a list of Inequalities associated to a given dominant 1PS to a string describing part of a latex tabular
    """
    caption=""
    if isinstance(V, KroneckerRepresentation):
        lambda_notation=False
        caption="Kronecker case in dimension $("
        for dk in V.G[:-1]:
            caption+=str(dk)+","
        caption=caption[:-1]
        caption+=")$"
    else:
        assert(isinstance(V, ParticleRepresentation))
        lambda_notation=True
        if isinstance(V, FermionRepresentation):
            caption="Ferminonic case $\\wedge^{"+str(V.particle_cnt)+"}\\mathbb{C}^{"+str(V.G[0])+"}$"
        if isinstance(V, BosonRepresentation):
            caption="Bosonic case $S^{"+str(V.particle_cnt)+"}\\mathbb{C}^{"+str(V.G[0])+"}$"
    grouped_ineqs=group_by_dom1PS(inequations)
    chaine='$\\begin{array}{|c| c |c|} \n \\hline \n \\textrm{dominant 1-PS} & \\textrm{Inequality} & w \\\\ \n \\hline'
    for taudom_list in grouped_ineqs:
        chaine+=Latex_string_of_cluster_dom1PS(taudom_list,lambda_notation,sgn) 
    chaine+='\\end{array}$'
    info=info_from_GV(V)+extra_info
    file0 = open('ineq_Latex-'+info+'.tex','w')
    
    file0.write("\\documentclass[12pt]{article} \n \\usepackage{amsmath,amssymb} \n \\usepackage{multirow} \n \\usepackage{graphicx} \n \n  \\begin{document}\n \n \\begin{table}\n \\caption{"+ caption+"}\n")
    file0.write(chaine)
    file0.write("\n \n \\end{table} \n \\end{document}")
    file0.close()


#######################################################################
#III. Python export
#######################################################################

def export_python(
        V: Representation,
        inequations: Sequence[Inequality],
        extra_info: str = ""
        ) -> None:
    info=info_from_GV(V)+extra_info
    file0 = open('ineq_Python-'+info+'.py','w')
    file0.write("#Inequalities selected for V of " + type(V).__name__ + " type with dimensions "+str(list(V.G)))
    if isinstance(V, ParticleRepresentation):
       file0.write("with number of particules = "+str(V.particle_cnt))
    file0.write("\n \nG=LinearGroup("+str(list(V.G))+") \n")
    chain = f"V = {type(V).__name__}(G"
    if isinstance(V, ParticleRepresentation):
       chain+=", "+str(V.particle_cnt)
    file0.write(chain+" )")
    file0.write("\n \nbrut_inequations=[")
    chain=""
    for ineq in inequations:
        chain+=str(ineq.wtau.flattened)+", \n"
    chain=chain[:-2]+" \n ] \n\n"   
    file0.write(chain)
    file0.write("#inequalities in our formated type Inequality \n")
    file0.write("inequalities=[Inequality.from_tau(Tau.from_flatten(brut_ineq,G)) for brut_ineq in brut_inequations] \n \n ")
    file0.close()


def export_terminal(
        V: Representation,
        inequations: Sequence[Inequality],
        extra_info: str = ""
        ) -> None:
    print(f"Computed {len(inequations)} inequalities{' (' + extra_info + ')' if extra_info else ''}:")
    for ineq in inequations:
        print(ineq)
    

def export_none(
        V: Representation,
        inequations: Sequence[Inequality],
        extra_info: str = ""
        ) -> None:
    pass


def export_many(
        formats: ExportFormat | Sequence[ExportFormat],
        V: Representation,
        inequations: Sequence[Inequality],
        extra_info: str = "",
        ) -> None:
    """ Export inequalities in given (possible multiple) format(s) """
    if isinstance(formats, str):
        formats = [formats]

    from .utils import to_literal
    for format in formats:
        match to_literal(ExportFormat, format):
            case "Normaliz":
                export_normaliz(V, inequations, extra_info=extra_info)
            case "LaTeX":
                export_latex(V, inequations, extra_info=extra_info)
            case "Python":
                export_python(V, inequations, extra_info=extra_info)
            case "Terminal":
                export_terminal(V, inequations, extra_info=extra_info)
            case "None":
                export_none(V, inequations, extra_info=extra_info)
                
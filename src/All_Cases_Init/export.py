#Python to normaliz and Python to latex
#TODO: integrate in the project


def Liste_to_Normaliz_string(liste,sgn=1): #sgn=1 ou -1 allows to change the sign, exchange >=0 and <=0
    """ converts a list of list of numbers to a string with
    *newline characters at the end of each sublist
    *a space between each number of a given sublist
    """
    chaine=''
    for l in liste:
        for x in l[:-1]: # no space at the end
            chaine+=str(sgn*x)+' '
        chaine+=str(sgn*l[-1])+'\n'
    return(chaine)

def info_from_GV(V):
    """ V a representation of a LinGroup G. It returns a nomalized string encoding the main information on G and V. To be used in the names of the input files given to Normaliz"""
    G=V.G
    info=''
    for i in G[:-1]:
       info+=str(i)+'-'
    info+=str(G[-1])+' '
    info+=V.type+' '
    if V.type in ["boson","fermion"]:
        info+=str(V.nb_part)+' '
    return info

def export_normaliz(V, inequations, r=None, equations=[], extra_info=""):
    """ V is a representation of a LinGroup G, 
    r is the rank of G aka the dimension of the ambient space of the equations (except possible customization, if our inequations contain a second member)
    inequations is a list of inequations (either of Inequality type or a list of coefficients)
    equations is an optional list of equations
    extra_info is a string that one may add to the standard name of the output file, e.g. to indicate the origin of our list of inequations """
    if r==None:
        r=V.G.rank
    if hasattr(inequations[0],"wtau"):
        True_inequations=[ineq.wtau.flattened for ineq in inequations]
    else:
        True_inequations=inequations
    info=info_from_GV(V)+extra_info
    print(True_inequations[0], info)
    fileO = open('ineq_Normaliz-'+info+'.in','w')
    fileO.write('amb_space '+str(r)+'\n\n')
    if len(equations)!=0:
        fileO.write('equations '+str(len(equations))+'\n\n')
        List_Eq=Liste_to_Normaliz_string(equations)
        fileO.write(List_Eq)
    fileO.write('\n'+'inequalities '+str(len(True_inequations))+'\n\n')
    fileO.write(Liste_to_Normaliz_string(True_inequations,-1)) #Our conventions so far work with inequalities of the form \sum a_i\lambda_i<=0 whil Normaliz standard is ">=0"
    fileO.close()


def Latex_string_of_tau(tau,lambda_notation=False, sgn=1):
    chaine=''
    if not(lambda_notation):
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

def Latex_string_of_cluster_dom1PS(inequations, lambda_notation=False, sgn=1): #sgn=1 ou -1 allows to change the sign, exchange >=0 and <=0
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
        chaine+=Latex_string_of_tau(Tau([tuple(wi) for wi in ineq.w[:s]],ineq.tau.G))
        chaine+='\\\\ \n \\cline{2-3} \n'
    chaine=chaine[:-15]
    chaine+='\\hline'
    return(chaine)
    
def group_by_dom1PS(inequations):
    """ from a list of Inequalities: groups them by same dominant 1 parameter subgroup
    """
    remaining_indices=[i for i in range(len(inequations))]
    grouped_ineqs=[]
    while remaining_indices!=[]:
        j=0
        taudom=inequations[remaining_indices[0]].tau
        taudom_ineqs=[]
        while j<len(remaining_indices):
            if inequations[remaining_indices[j]].tau==taudom:
                taudom_ineqs.append(inequations[remaining_indices[j]])
                del(remaining_indices[j])
            else:
                j+=1
        grouped_ineqs.append(taudom_ineqs)
    return grouped_ineqs

    
#TODO check fermion
#TODO create a file
#TODO make a caption

def export_latex(V, inequations, sgn=1, extra_info=""): #sgn=1 ou -1 allows to change the sign, exchange >=0 and <=0
    """ converts a list of Inequalities associated to a given dominant 1PS to a string describing part of a latex tabular
    """
    if V.type=='kron':
        lambda_notation=False
    else:
        lambda_notation=True
    grouped_ineqs=group_by_dom1PS(inequations)
    chaine='$\\begin{array}{|c| c |c|} \n \\hline \n \\textrm{dominant 1-PS} & \\textrm{Inequality} & w \\\\ \n \\hline'
    for taudom_list in grouped_ineqs:
        chaine+=Latex_string_of_cluster_dom1PS(taudom_list,lambda_notation) 
    chaine+='\\end{array}$'
    info=info_from_GV(V)+extra_info
    file0 = open('ineq_Latex-'+info+'.tex','w')
    file0.write("\\documentclass[12pt]{article} \n \\usepackage{amsmath} \n \\usepackage{multirow} \n \\usepackage{graphicx} \n \n  \\begin{document}\n \n")
    file0.write(chaine)
    file0.write("\n \n \\end{document}")
    file0.close()

    


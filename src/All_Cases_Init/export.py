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
    return chaine 

def Latex_string_of_cluster_1PS(taudom, inequations, lambda_notation=False, sgn=1): #sgn=1 ou -1 allows to change the sign, exchange >=0 and <=0
    """ converts a list of Inequalities (associated to a given dominant 1PS  taudom) to a string describing part of a latex tabular
    """
    n=len(inequations)
    chaine='\\multirow{'+str(n)+'}{*}{'
    chaine+=Latex_string_of_tau(taudom, False, sgn)+' } ' 
    for ineq in inequations:
        chaine+=' & '+Latex_string_of_tau(ineq.wtau, lambda_notation, sgn)+' & '
        chaine+='\\\\ \n \\cline{2-3} \n'
    chaine=chaine[:-15]
    chaine+='\\hline'
    return(chaine)
    
def export_latex(V, inequations, lambda_notation=False, sgn=1): #sgn=1 ou -1 allows to change the sign, exchange >=0 and <=0
    """ converts a list of Inequalities associated to a given dominant 1PS to a string describing part of a latex tabular
    """
    if V.type=='kron':
        lambda_notation=False
    else:
        lambda_notation=True
    chaine='$\\begin{array}{|c| c |c|} \n \\hline \n \\textrm{dominant 1-PS} & \\textrm{Inequality} & w \\\\ \n \\hline'
    for ineq in inequations:
        chaine+=Latex_string_of_tau(ineq.tau) 
        chaine+='& '+Latex_string_of_tau(ineq.wtau, lambda_notation, sgn)
        chaine+='&'
        chaine+='\\\\ \n \\hline \n'
    chaine+='\\end{array}$'
    return(chaine)

    



#example of input
#info=info_from_GV(G,V)
#export_normaliz(info,G.rank,[ineq.wtau.flattened for ineq in Birational_Ineq48])

#TODO be more consistent with real objects in Normaliz:
#if amb_space is 8, an inequality 2a+3b<=8 translates to 
"""
amb_space 2

inhom_inequalities 14            
-2 -3 8

#What x1 a_1+...xn a_n +k >=0 translates to (x1,...,xn, k)
#same for equations

#Another possible writing:
"""
amb_space 2

constaints 1
2 3 <= 8
"""


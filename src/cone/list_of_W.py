__all__ = (
    #"ListW_Mod",
    "List_Inv_Ws_Mod",
    "Check_Rank_Tpi",
)

import operator
import itertools
import numpy as np
from typing import Tuple
from collections import defaultdict


from .typing import *
from .permutation import *
from .root import *
from .vector_chooser import *
from .rings import matrix
from .tau import Tau
from .representation import *
from .inequality import Inequality
from .utils import *
from sage.all import Partition,Partitions
from copy import copy, deepcopy

def List_Inv_Ws_Mod(tau : Tau,V: Representation) -> Iterable[dict[int,list[Root]]]:
    """
    Return the list (as iterable) of inversion sets compatible W^{P(tau) and with the C^*-module V^{tau>0}.
    The output is an iterable of dictionnaries int -> list(Root)
    This function initializes the contraints and start the recursive part.
    """
    print("llllllllllllllllllllllll",tau,"ttt",tau.reduced)
    lG=list(tau.G)
    while lG and lG[-1] == 1:
        lG.pop()
    s = len(lG) # Number of non-trivial symmetric groups  
    nbs_blocks=[x-1 for x in tau.reduced.small_d] # list of the number of blocks on the first row 
    max_nb_blocks = max(nbs_blocks) # max number of blocks on the first row (or equivalently last column)
    
    weights_grid = np.empty((s,max_nb_blocks, max_nb_blocks), dtype=int)
    size_grid = np.empty((s,max_nb_blocks, max_nb_blocks), dtype=object)
    #size_grid = np.empty((s,max_nb_blocks, max_nb_blocks), dtype=[('x', int), ('y', int)])
    inner_grid = np.empty((s,max_nb_blocks, max_nb_blocks), dtype=object)
    outer_grid = np.empty((s,max_nb_blocks, max_nb_blocks), dtype=object)
    init_inv = np.empty((s,max_nb_blocks, max_nb_blocks), dtype=object)
    # List of blocks that will still contribute to the weight p
    List_redroots_next,Dic_tau_redroots = Init_pos_redroots(tau)
    current_pos=List_redroots_next[0]   # copy of GL,i,j
    #Initialization of target_weights to len(grading_positive_weights)
    target_weights = {key: len(value) for key, value in tau.positive_weights(V).items()}
    for k in range(s):
        print(nbs_blocks[k])
        for i in range(nbs_blocks[k]):
            for j in range(i,nbs_blocks[k]):
                size_grid[k,i,j] = [tau.reduced.mult[k][i],tau.reduced.mult[k][j+1]]
                weights_grid[k,i,j] = tau.reduced.values[k][i]-tau.reduced.values[k][j+1]
                inner_grid[k,i,j] = [0]*(tau.reduced.mult[k][i])
                outer_grid[k,i,j] = [tau.reduced.mult[k][j+1]]*(tau.reduced.mult[k][i])
    #print("tttt",inner_grid,outer_grid,target_weights)
    result= List_Inv_W_Mod_rec(nbs_blocks,init_inv,weights_grid,size_grid,inner_grid,outer_grid,target_weights,List_redroots_next,Dic_tau_redroots)
    print("ssssssssssssssssssss",result)
    return result

def Init_pos_redroots(tau:Tau)->(Iterable[Tuple[int]],dict[int,Iterable[Tuple[int]]]):
    taured=tau.reduced
    list_result=[]
    dic_result = defaultdict(list)
    for k,tau1 in enumerate(taured.values):
        for i in range(len(tau1)-1):
            for j in range(i,len(tau1)-1):
                dic_result[tau1[j-i]-tau1[j+1]].append((k,j-i,j))
                list_result.append((k,j-i,j))
                #print(list_result)
    return (list_result,dic_result)

def adjust_inner_outer_ijk(inner_ik,outer_ik,entry_ij,entry_jk):
    inner_ik_new=[]
    outer_ik_new=[]
    mi=len(entry_ij)
    mj=len(entry_jk)
    for i1 in range(mi):
         if entry_ij[i1]==0:
             inner_ik_new.append(inner_ik[i1])
         else:
             inner_ik_new.append(max(inner_ik[i1],entry_jk[mj-entry_ij[i1]]))
         if entry_ij[i1]==mj:
             outer_ik_new.append(outer_ik[i1])
         else:
             outer_ik_new.append(min(outer_ik[i1],entry_jk[mj-entry_ij[i1]-1]))
    return (inner_ik_new,outer_ik_new)

def List_Inv_W_Mod_rec(nbs_blocks : list[int],current_inv,weights_grid,size_grid,inner_grid,outer_grid,target_weights,List_redroots_next,Dic_tau_redroots):
    ## We first determine the possible partitions to take in position current_pos
    # Possible lengths
    if len(List_redroots_next)==0:
        return [Table_part_2_inv_list(nbs_blocks,size_grid,weights_grid,current_inv)]
    current_pos=List_redroots_next[0]
    k,i,j=current_pos
    result=[]
    p = weights_grid[*current_pos]
    #print("ggg", current_pos, weights_grid)
    Dic_tau_redroots[p]=Dic_tau_redroots[p][1:]
    List_redroots_next=List_redroots_next[1:]
    if len(Dic_tau_redroots[p]) == 0: # On peut améliorer en tenant compte de inner et outer
        if p in target_weights.keys():
            target_lengths = [target_weights[p]]
        else:
            target_lengths=[0]
    else :
        target_lengths = [l for l in range(target_weights[p]+1)]
    # We run over these lengths
    if nbs_blocks[k]==0:
        return  List_Inv_W_Mod_rec(nbs_blocks, current_inv, weights_grid, size_grid, inner_grid, outer_grid, target_weights, List_redroots_next,Dic_tau_redroots)
    else:
        for l in target_lengths :
        # Je crois que la ligne suivante ne marche pas car notre partition est concurrente de celle de Sage
        #if inner_grid[*current_pos]==None: #case of a single block in tau[k]
        #    partit=[[]]
            partit=[list(part) for part in Partitions(l, inner=inner_grid[*current_pos], outer=outer_grid[*current_pos])]
            current_part_list=[part+[0]*(len(inner_grid[*current_pos])-len(part)) for part in partit]
            #current_part_list=[list(part)+[0]*(len(inner_grid[*current_pos])-len(part)) for part in Partitions(l, inner=inner_grid[*current_pos], outer=outer_grid[*current_pos])]
            s= len(weights_grid) 
            for mu in current_part_list :
                next_inv = current_inv.copy() # copy to avoid confusion
                next_inv[*current_pos]=mu
                # Ajust inner and outer
                inner_grid_next=inner_grid.copy()    #deepcopy?
                outer_grid_next=outer_grid.copy()
                ## above current_pos
                for a in range(i):
                    #print("wwww",next_inv,inner_grid,outer_grid,k,a,i,j)
                    inner_grid_next[k,a,j],outer_grid_next[k,a,j] = adjust_inner_outer_ijk(inner_grid[k,a,j],outer_grid[k,a,j],next_inv[k,a,i-1],next_inv[k,i,j])
                """
                    new_inner=[]
                    new_outer=[]
                    for i in range(size_grid[current_pos[0],a,current_pos[2]][0]): # Run over the rows of the block
                        L_in_i=[inner_grid[current_pos[0],a,current_pos[2]][i]]
                        L_out_i=[outer_grid[current_pos[0],a,current_pos[2]][i]]
                        for b in range(current_pos[2]):
                            lam_i=next_inv[current_pos[0],a,b][i]
                            L_in_i.append(mu[m-lam_i])
                            L_out_i.append(mu[m-lam_i-1])
                        new_inner.append(max(L_in_i))
                        new_outer.append(min(L_out_i))
                        if max(L_in_i)>min(L_out_i):
                            print('Incompatible inner and outer')
                    inner_grid_next[size_grid[current_pos[0],a,current_pos[2]]] = new_inner
                    outer_grid_next[size_grid[current_pos[0],a,current_pos[2]]] = new_outer
                    ## right of current_pos"""
                """for b in range(j+1,nbs_blocks[k]):
                    inner_grid_next[k,i,b],outer_grid_next[k,i,b] = adjust_inner_outer_ijk( inner_grid[k,i,b], outer_grid[k,i,b], next_inv[k,i,j], next_inv[k,j,b])
                    new_inner=[]
                    new_outer=[]
                    print(size_grid,current_pos,b)
                    print(size_grid[current_pos[0],current_pos[1],b])
                    for i in range(size_grid[current_pos[0],current_pos[1],b][0]): # Run over the rows of the block
                        print("aabb",inner_grid,i)
                        print(inner_grid[current_pos[0],current_pos[1],b],i)
                        current_inner=inner_grid[current_pos[0],current_pos[1],b]
                        if i>=len(current_inner):
                            L_in_i=[0]
                        else:
                           L_in_i=[inner_grid[current_pos[0],current_pos[1],b][i]]
                           L_out_i=[outer_grid[current_pos[0],current_pos[1],b][i]]
                        for a in range(current_pos[1]+1,nbs_blocks[current_pos[0]]):
                            lam=next_inv[current_pos[0],a,b]
                            L_in_i.append(lam[m-mu[i]])
                            L_out_i.append(lam[m-mu[i]-1])    
                        new_inner.append(max(L_in_i))
                        new_outer.append(min(L_out_i))
                        if max(L_in_i)>min(L_out_i):
                            print('Incompatible inner and outer')
                    inner_grid_next[size_grid[current_pos[0],current_pos[1],b]] = new_inner
                    outer_grid_next[size_grid[current_pos[0],current_pos[2],b]] = new_outer"""
                # Ajust target_weights
                target_weights_next=target_weights.copy() # copy to avoid confusion
                if p in target_weights_next.keys():
                    target_weights_next[p]-=l
                # Exit if not possible : inner, outer incompatible with target_weights
                for p1 in target_weights_next.keys():
                    MAX_mult=sum(sum(outer_grid_next[*free_pos]) for free_pos in Dic_tau_redroots[p1])
                    MIN_mult=sum(sum(inner_grid_next[*free_pos]) for free_pos in Dic_tau_redroots[p1])
                    if MAX_mult < target_weights_next[p1] or MIN_mult > target_weights_next[p1] :
                        print("a")
                        return []
                    # Recursive call
                result+= List_Inv_W_Mod_rec(nbs_blocks, next_inv, weights_grid, size_grid, inner_grid_next, outer_grid_next, target_weights_next, List_redroots_next, Dic_tau_redroots)
    print("rrrrrrrrrrrrrrr",result)
    return result
        
"""        
def Next_pos(current_pos : Tuple[int],nbs_blocks : list[int]) -> Tuple[int]:
    
    return the next position in the set of blocs. 
    In a given pos we strart with the longuest diagonal and so on. 
    
    pos,i,j=current_pos
    if i<nbs_blocks[pos]-1 and j<nbs_blocks[pos]-1: # progress on the diagonal
        return (pos,i+1,j+1)
    elif j-i<nbs_blocks[pos]-1: # go to the next diagonal
        return (pos,0,j-i+1)
    else : # change pos
        return (pos+1,0,0)
"""


def Table_part_2_inv_list(nbs_blocks : list[int],size_grid,weights_grid,T) -> dict[int, list[Root]]:
    """
    T is a 3-dimensional table of partitions encoding a set of inversions. 
    T[pos,i,j] (with j>=i) row i column j
    nbs_blocks[pos] encode the size of the triangle T[pos,*,*]
    size_grid[pos,i,j] is a pair of integers encoding the size of the block
    weights_grid[pos,i,j] is the weight (for the action of tau) on the roots associated to bowes in the block.
    The function return a dictionnary weight -> list of roots.
    """
    result = defaultdict(list)
    s = len(nbs_blocks)
    for pos,a in enumerate(nbs_blocks):
        shift_i=0
        for i in range(a):
            shift_j=size_grid[pos,0,0][0]
            for j in range(i,a):
                for i1,l in enumerate(T[pos,i,j][::-1]):
                    for j1 in range(l):
                        result[weights_grid[pos,i,j]].append(Root(pos,shift_i+i1,shift_j+j1)) 
                shift_j+=size_grid[pos,j,j][1]
            shift_i+=size_grid[pos,i,i][0]
    return result    

def inversion_set_to_permutation(n, inversions):
    """
    Reconstructs a permutation from its inversion set.
    
    Parameters:
        n (int): The size of the permutation (0 to n-1).
        inversions (set of tuples): The set of inversions, where each inversion is a tuple (i, j) with i < j.
    
    Returns:
        list: The reconstructed permutation.
    """
    # Initialize the inversion count for each element
    inv_count = [0] * n
    for i, j in inversions:
        inv_count[j] += 1
    
    # Construct the permutation by placing numbers from largest to smallest
    w = [-1] * n
    available_positions = list(range(n))
    for num in range(n-1, -1, -1):
        position = available_positions.pop(inv_count[num])
        w[position] = num
    
    return Permutation(w)


def Check_Rank_Tpi(ineq : Inequality, V: Representation, method: Method) -> bool :
    """
    Check if Tpi is inversible at a general point of V^\tau. 
    General means randon if the method is probabilist and formal if the method is symbolic.
    The matrix being block trinagular, the function check successively the diagonal blocks.
    """
    tau=ineq.tau
    G = tau.G 
    # Ring depending on the computational method
    if method == "probabilistic":
        ring = V.QI
    elif method == "symbolic":
        ring = V.QV
    else:
        raise ValueError(f"Invalid value {method} of the computation method")
    zero_weights = tau.orthogonal_weights(V)
    v = point_vect(zero_weights, V, ring, bounds=(-100, 100))
    gw = tau.grading_weights(V)
    gr = ineq.gr_inversions 
    for x in sorted(gr.keys(),reverse=True): # Run over the possible values of tau.scalar(root) for root inversion of w
        M=matrix(ring,len(gr[x]))
        for col,root in enumerate(gr[x]): # List of roots such that tau.scalar(root)=x
               uv=V.action_op_el(root,v)
               for row, chi in enumerate(gw[x]): # List of weights such that tau.scalar(chi)=x 
                   M[row,col]=uv[V.index_of_weight(chi)]
        rank_M = M.change_ring(ring.fraction_field()).rank()
        if rank_M<M.nrows():
               return(False)
    return(True)	       



    

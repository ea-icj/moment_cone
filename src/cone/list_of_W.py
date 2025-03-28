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
from sage.all import Partition, Partitions
from copy import copy, deepcopy

def List_Inv_Ws_Mod(tau : Tau,V: Representation) -> Iterable[dict[int,list[Root]]]:
    """
    Returns the list of inversion sets compatible with W^{P(tau) and with the C^*-module V^{tau>0}.
    The output is an iterable of dictionnaries int -> list(Root)
    This function initializes the contraints and start the recursive part.
    """
    lG=list(tau.G)
    while lG and lG[-1] == 1:
        lG.pop()
    s = len(lG) # Number of non-trivial symmetric groups  
    nbs_blocks=[x-1 for x in tau.reduced.small_d] # for each (k in [0,s-1]), stores the number of blocks-1 of the centralizer of tau.component[k]  
    max_nb_blocks = max(nbs_blocks)
    weights_grid = np.empty((s,max_nb_blocks, max_nb_blocks), dtype=int) #note that indices (k,i,j) in this array correspond to bloc with indices (i,j+1) in the k-th component.  
    sizes_blocks = np.empty((s,max_nb_blocks+1), dtype=int) #sizes of the bloc with indices (i,j+1) in the k-th component can be recovered via (sizes_blocs[k,i],sizes_blocs[k,j])
    #size_grid = np.empty((s,max_nb_blocks, max_nb_blocks), dtype=[('x', int), ('y', int)])
    init_inv = np.empty((s,max_nb_blocks, max_nb_blocks), dtype=object) #each entry of init_inv will store a partition (in fine, encoding the roots corresponding to the inversions of a w in W^P(tau) compatible with the C^*-module V^{tau>0}.)
    inner_grid = np.empty((s,max_nb_blocks, max_nb_blocks), dtype=object) # entries of inner_grid and outer_grid are partitions bounding the possible entries of init_inv
    outer_grid = np.empty((s,max_nb_blocks, max_nb_blocks), dtype=object)
    List_pos,Dic_tau_redroots = Init_pos_redroots(tau)  
    current_pos=List_pos[0]   # first position that will be explored
    #Initialization of target_weights to len(grading_positive_weights)
    target_weights = {key: len(value) for key, value in tau.positive_weights(V).items()}
    for k in range(s):
        sizes_blocks[k,nbs_blocks[k]]=tau.reduced.mult[k][nbs_blocks[k]]
        for i in range(nbs_blocks[k]):
            sizes_blocks[k,i] = tau.reduced.mult[k][i]
            for j in range(i,nbs_blocks[k]):
                weights_grid[k,i,j] = tau.reduced.values[k][i]-tau.reduced.values[k][j+1]
                inner_grid[k,i,j] = [0]*(tau.reduced.mult[k][i])
                outer_grid[k,i,j] = [tau.reduced.mult[k][j+1]]*(tau.reduced.mult[k][i])
    #print("tttt",inner_grid,outer_grid,target_weights)
    sym = tau.outer 
    sum_sym=[sum(sym[:i]) for i in range(len(sym)+1)] 
    test_inc= True #in recursive part, when True, we make checks in order to explore modulo symmetries of tau. When turned to false in some recursive instance, check is no longer needed for some component of tau since monotony requirements modulo symmetries of tau are already met. 
    result= List_Inv_W_Mod_rec(nbs_blocks, sizes_blocks, init_inv, weights_grid, inner_grid, outer_grid, target_weights, List_pos, Dic_tau_redroots,sum_sym,test_inc)
    #print("result global function_List_Inv_WS_Mod",result)
    return result

def Init_pos_redroots(tau:Tau)->(Iterable[Tuple[int]],dict[int,Iterable[Tuple[int]]]):
    """
    dic_result lists, for each weight p, the indices corresponding to blocks of U with weight p. list_result is a single list of these indices, ordered in such a way that, when exploring it, inner_grid and outer_outer_grid can take into account all constraints.
    """
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

def adjust_inner_outer_ijk(inner_ik,outer_ik,entry_ij,entry_jk, mi: int, mj: int):
    """
    encodes the constraints in bloc (i,k) (in some component) of inner_grid and outer_grid  from already fixed inversions in blocs (i,j) and (j,k).
    """
    inner_ik_new=[]
    outer_ik_new=[]
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

def List_Inv_W_Mod_rec(nbs_blocks : list[int],sizes_blocks, current_inv,weights_grid,inner_grid,outer_grid,target_weights,List_pos,Dic_tau_redroots,sum_sym : tuple[int],test_inc : bool):
    """
    recursive part for computing the list of inversion sets compatible with W^{P(tau) and with the C^*-module V^{tau>0}.
    
    """
    if len(List_pos)==0: #last position already hit, we thus have a set of inversions compatible with constraints; converted to a format compatible with attribute ._gr_inversions of an element of class Inequality
        return [Table_part_2_inv_list(nbs_blocks,sizes_blocks,weights_grid,current_inv)]
    current_pos=List_pos[0]
    List_pos_next=List_pos[1:]
    k,i,j=current_pos
    test_prev=(current_pos[0] not in sum_sym) #In the current instance of List_Inv_W_Mod_rec, we will consider comparisons with the previous bloc only when it satisfies tau.components[k-1] = tau.components[k] (until test_inc=False) in order to work modulo symmetries of tau.
    result=[]
    p = weights_grid[*current_pos] 
    Dic_tau_redroots_next=Dic_tau_redroots.copy()
    Dic_tau_redroots_next[p]=Dic_tau_redroots_next[p][1:]  #remove current position from positions yet unsettled.
    # Possible lengths
    if p in target_weights.keys():
        if len(Dic_tau_redroots_next[p]) == 0: #TODO On peut améliorer en tenant compte de inner et outer
            target_lengths = [target_weights[p]]
        else :
            target_lengths = [l for l in range(target_weights[p]+1)]
    else:
        target_lengths=[0]
    # We will explore the different partitions of each of these lengths
    for l in target_lengths :
        partit=[list(part) for part in Partitions(l, inner=inner_grid[*current_pos], outer=outer_grid[*current_pos])]
        current_part_list=[part+[0]*(len(inner_grid[*current_pos])-len(part)) for part in partit] #padding the partitions (still useful?)
        for mu in current_part_list :
            # If test_inc keep only mu that are bigger or equal
            if test_inc and test_prev:
                 mu_ref=current_inv[k-1,i,j]
                 if mu < mu_ref:
                     continue # skip this mu
                 if mu > mu_ref:
                     test_inc = False # for some lexicographic order, we have ensured that w[k]>w[k-1]
            next_inv = current_inv.copy() # copy to avoid modifications in current_inv from the exploration of the current branch
            next_inv[*current_pos]=mu  # setting our entry to mu
            # Ajust target_weights
            target_weights_next=target_weights.copy() # copy to avoid confusion
            if p in target_weights_next.keys():
                target_weights_next[p]-=l
            # Ajust inner and outer and exit if incompatibility (inner bigger than outer).
            to_continue=True
            inner_grid_next=inner_grid.copy()    #deepcopy?
            outer_grid_next=outer_grid.copy()
            ## above current_pos
            bound_a=max(2*i-j-1,0)
            for a in range(bound_a,i):
                inner_grid_next[k,a,j],outer_grid_next[k,a,j] = adjust_inner_outer_ijk(inner_grid[k,a,j],outer_grid[k,a,j], next_inv[k,a,i-1], next_inv[k,i,j], sizes_blocks[k][a], sizes_blocks[k][i])
            bound_b=min(2*j-i+1, nbs_blocks[k])
            for b in range(j+1,bound_b):
                inner_grid_next[k,i,b],outer_grid_next[k,i,b] = adjust_inner_outer_ijk(inner_grid[k,i,b],outer_grid[k,i,b], next_inv[k,i,j], next_inv[k,j+1,b], sizes_blocks[k][i], sizes_blocks[k][j+1])
            # Exit if not possible : inner, outer incompatible with target_weights
            for p1 in target_weights_next.keys():
                MAX_mult=sum(sum(outer_grid_next[*free_pos]) for free_pos in Dic_tau_redroots_next[p1])
                MIN_mult=sum(sum(inner_grid_next[*free_pos]) for free_pos in Dic_tau_redroots_next[p1])
                if MAX_mult < target_weights_next[p1] or MIN_mult > target_weights_next[p1] :
                    to_continue=False
            # Recursive call if the previous tests were passed
            if to_continue:
                if len(List_pos_next)!=0 and List_pos_next[0][0]>k: #Reinit test_inc when a new bloc appears 
                    test_inc=True 
                result+= List_Inv_W_Mod_rec(nbs_blocks, sizes_blocks, next_inv, weights_grid, inner_grid_next, outer_grid_next, target_weights_next, List_pos_next, Dic_tau_redroots_next,sum_sym,test_inc)
    #print("partial function result,next_inv,current_pos",result, next_inv,current_pos)
    return result


def Table_part_2_inv_list(nbs_blocks : list[int],sizes_blocks,weights_grid,T) -> dict[int, list[Root]]:
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
            shift_j=sum([sizes_blocks[pos,l] for l in range(i+1)])
            for j in range(i,a):
                for i1,l in enumerate(T[pos,i,j][::-1]):
                    for j1 in range(l):
                        result[weights_grid[pos,i,j]].append(Root(pos,shift_i+i1,shift_j+j1)) 
                shift_j+=sizes_blocks[pos,j+1]
            shift_i+=sizes_blocks[pos,i]
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



    

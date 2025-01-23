from sage.all import QQ,vector,binomial
from functools import cached_property

import itertools
import operator

#from .rings import ZZ
from .typing import *
from .utils import prod, orbit_symmetries
from .group import *
from sage.all import ZZ
#from .rep import *

#from __future__ import annotations


__all__ = (
    "Weight",
)


class Weight:
    """ 
    A weight is a the data of a linear group G and a character of its maximal torus as a vector of length G.rank with coefficients in Z 
    idx: Optional[int] # Index in the lexicographical order
    """  

    def __init__(self, G : LinGroup,as_vector : Optional[vector]=None, mult: Optional[int]=1,idx: Optional[int] = None, as_list: Optional[list[int]]=None,as_list_of_list: Optional[list[list[int]]]=None):
        """
        Initializes an instance of Weight.
        """
        if not isinstance(G, LinGroup):
            raise TypeError("G must be an instance of the LinGroup class.")
        self.G = G
        self._as_vector = as_vector
        self.mult = mult
        self._idx = idx
        self.as_list = as_list
        self.as_list_of_list = as_list_of_list

    def __repr__(self) -> str:
        if self.as_list!=None:
            return 'One epsilon_i by copy of GL: '+ str(self.as_list)
        elif  self.as_list_of_list!=None:
            return 'A list of epsilon_i for each copy of GL: '+ str(self.as_list_of_list)
        else :
            return 'Coordinates of Weight: '+str(self.as_vector)

    def __eq__(self,other: "Weight") -> bool:
        
        # TODO : J'ai supprimmer ces deux lignes car elles renvoyaient toujours faux. Il faudrait les remettre en corrigent, mais je n'ai pas réussi.
        # L'erreur était à l'intérieur de hyperplane_candidats pour kron quand on demande si chi est dans St.indeterminate
        
        #if not isinstance(other, Weight) or self.G != other.G :
        #    return(False)
        
        if self.as_list!=None and other.as_list!=None:
            return self.as_list == other.as_list
        if self.as_list_of_list!=None and other.as_list_of_list!=None:
            return self.as_list_of_list == other.as_list_of_list
        return self.as_vector == other.as_vector
    
    #@staticmethod
    def idx(self,V : "Representation") -> int:
        if self._idx != None :
            return(self._idx)
        if V.type == 'kron':
            res=0
            for i,x in enumerate(self.as_list):
                res+=x*prod(V.G[i+1:])
            return(res)
        if V.type == 'fermion':
            S=self.as_list_of_list[0]
            k=V.nb_part
            id=sum([binomial(V.G.rank-j-1,k-1) for j in range(S[0])])# Subset with smaller first element
            for i,p in enumerate(itertools.pairwise(S)):
                id+=sum([binomial(V.G.rank-j-1,k-i-2) for j in range(p[0]+1,p[1])]) # With j in position i+1 and equal for smaller indices       
            return(id)
        else :
            for i,chi in enumerate(V.all_weights):
                if chi.as_vector==self.as_vector :
                    return i
        
        ## TODO cas bosons
        
    @cached_property
    def as_vector(self) -> vector:
        if self._as_vector != None:
            return self._as_vector
        if self.as_list!=None:
           v=vector(ZZ,self.G.rank)
           shift=0
           for i,x in enumerate(self.as_list) :
                v[shift+x]+=1
                shift+=self.G[i]
           return v
        elif self.as_list_of_list != None:
           v=vector(ZZ,self.G.rank)
           shift=0
           for i,x in enumerate(self.as_list_of_list) :
                for y in x :
                    v[shift+y]+=1
                shift+=self.G[i]
           return v
       

    def __le__(self, other: "Weight") -> bool:
        """ 
        Used in type kron
        Implementation of self <= other (partial ordering)
        """
        if self.as_list != None and other.as_list != None :
            return all(ws >= wo for ws, wo in zip(self.as_list, other.as_list))
        if self.as_list_of_list != None and other.as_list_of_list != None :
            return all(all(ws >= wo for ws, wo in zip(ls, lo)) for ls, lo in zip(self.as_list_of_list, other.as_list_of_list))
        #return self.leq(other)
    
    def leq(self,other,sym : list[int]=None) -> bool:
        """
        Used in type boson and fermion
        Dominance order (partial sums) on each block given by sym. 
        The convention is 
                     chi2.leq(chi,sym) and tau dominant => <chi2,tau> <= <chi2,tau>
        sym precise the sense of tau strictly dominant (decreasing on each block of sizes given by sym).
        TODO : check that is works. Je ne suis pas sûr du sens vu que ça ne marche pas pour Kron.
        """
        
        if sym == None:
            if self.as_list != None and other.as_list != None :
                return all(ws >= wo for ws, wo in zip(self.as_list, other.as_list))
            if self.as_list_of_list != None and other.as_list_of_list != None :
                return all(all(ws >= wo for ws, wo in zip(ls, lo)) for ls, lo in zip(self.as_list_of_list, other.as_list_of_list))
            
        v=other.as_vector-self.as_vector
        if sym == None :
            sym=self.G
        shift=0
        for s in sym:
            partial_sum=0
            for i in range(s):
                partial_sum+=v[shift+i]
                if partial_sum<0 :
                    return False
            if partial_sum!=0 :
                return False
            shift+=s
        return True

    def orbit_symmetries(self,S : list[int]) -> Iterable["Weight"]:
        if self.as_list != None :
            orbit=orbit_symmetries(self.as_list,S)
            for sym_w in orbit:
                yield Weight(self.G,as_list=list(sym_w))
    

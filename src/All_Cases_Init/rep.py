from functools import cached_property
import itertools
from sage.all import Partition, matrix, binomial

from math import prod
from .typing import *
from .rings import Vector,ZZ,QQ, I,vector,Matrix
from .combi import *
from .partition import *
from .group import *
from .weight import *
from .root import *


class Representation:
    def __init__(self, G : LinGroup, type : str, nb_part=None):
        """
        Initializes an instance of Representation.
        
        :param G: An instance of the LinGroup class.
        :param type: A string, must be 'kron', 'fermion', or 'boson'.
        :param nb_part: An integer, required if type is 'fermion' or 'boson'.

        Example : 
        >>> V=Representation(G,'fermion',2)
        >>> V
        Representation(G=GL(5), type='fermion', number of particules=2)
        >>> V.dim
        10
        >>> V.G
        GL(5)
        >>> V.type
        'fermion'
        """
        
        if not isinstance(G, LinGroup):
            raise TypeError("G must be an instance of the Group class.")
        
        if type not in ['kron', 'fermion', 'boson']:
            raise ValueError("type must be 'kron', 'fermion', or 'boson'.")
        
        self.G = G
        self.type = type

        if type in ['fermion', 'boson']:
            if nb_part is None :
                raise ValueError("nb_part must be an integer for types 'fermion' or 'boson'.")
        self.nb_part = nb_part
        
    
    def __repr__(self):
        """
        String representation of the instance.
        """
        if self.type in ['fermion', 'boson']:
            return f"Representation(G={self.G}, type='{self.type}', number of particules={self.nb_part})"
        return f"Representation(G={self.G}, type='{self.type}')"

    @cached_property
    def dim_cone(self) -> int:
        """
        Expected dimension of the cone. To be checked with Stabilizer of K.
        
        >>> V=Representation(G,'fermion',2)
        >>> V.dim_cone
        5
        """
        
        if self.type == 'kron':
            return self.G.rank-len(self.G)+1 
        else :
            return self.G.rank

    @cached_property
    def dim(self) -> int:
        """
        Dimension of V.
        """
        if self.type == 'kron':
            return prod(self.G) 
        elif self.type == 'fermion':
            return binomial(self.G.rank,self.nb_part)
        else :
            return binomial(self.G.rank+self.nb_part-1,self.G.rank-1)
        
    @cached_property    
    def all_weights(self) -> list["Weight"] : # TODO : j'ai du remplacer itérable ar liste car iterable ne marche qu'une fois ?
        """
        Create de the list of weights of T acting on V. 

      >V=Representation(G,'fermion',2)
      >L=[chi for chi in V.all_weights]
      >L
      [A list of epsilon_i for each copy of GL: [[0, 1]],
      A list of epsilon_i for each copy of GL: [[0, 2]],
      A list of epsilon_i for each copy of GL: [[0, 3]],
      A list of epsilon_i for each copy of GL: [[0, 4]],
      A list of epsilon_i for each copy of GL: [[1, 2]],
      A list of epsilon_i for each copy of GL: [[1, 3]],
      A list of epsilon_i for each copy of GL: [[1, 4]],
      A list of epsilon_i for each copy of GL: [[2, 3]],
      A list of epsilon_i for each copy of GL: [[2, 4]],
      A list of epsilon_i for each copy of GL: [[3, 4]]]
        """
        L=[]
        if self.type == 'kron' :
            for i,w in enumerate(itertools.product(*(range(di) for di in self.G))):
                L.append(Weight(self.G,as_list=list(w),idx=i))
        elif self.type == 'fermion' :
            for i,w in enumerate(itertools.combinations(range(self.G[0]),self.nb_part)):
                L.append(Weight(self.G,as_list_of_list=[list(w)],idx=i))
        else :
            for i,w in enumerate(itertools.combinations(range(self.G[0]+self.nb_part-1),self.G[0]-1)) :
                v=vector(ZZ,self.G.rank)
                v[0]=w[0]
                AsL=w[0]*[0]
                for j,p in enumerate(itertools.pairwise(w)):
                    v[j+1]=p[1]-p[0]-1
                    AsL+=v[j+1]*[j+1]
                v[-1]=self.G.rank+self.nb_part-w[-1]-2
                AsL+=v[-1]*[self.G.rank-1]
                L.append(Weight(self.G,as_vector=v,as_list_of_list=[AsL],idx=i))
        return(L)        

    @cached_property 
    def weights_mod_outer(self) -> Iterable["Weight"]:
        """
        Returns all decreasing weights modulo the symmetries of V that is of G.
        
        It assumes that G is arranged so that blocks are contiguous.
        """
        # Using Partition to generate all decreasing weights within a block.
        # The whole weights will be defined as the Cartesian product of the weights for each block.
        from .partition import Partition
        from .utils import group_by_block

        # For fermion and boson outer is trivial
        if self.type != 'kron' :
            return self.all_weights
        
        def pad(p: Partition, l: int) -> tuple[int, ...]:
            return p.pad(l)
        
        block_weights = tuple(
            tuple( # Converting to tuple (instead of a map) seems necessary to keep the right hi (FIXME)
                p.pad(hi) # Adding trailing zeros if necessary
                for p in Partition.all_of_length(hi, di - 1)
            )
            for di, hi in group_by_block(self.G) # Compress returns (value, multiplicity) for each block of G
        )

        for w in itertools.product(*block_weights):
            Lw=[]
            for x in w :
                Lw+=list(x)
            yield Weight(self.G,as_list=list(vector(sum(w, start=())))) # Summing tuples is concatenating them 

    def weights_of_S(self, p : Partition) -> Iterable["Weight"] : # Could be improved
        """
        Create de the list of weights of S\subset T acting on V. With multiplicities.
        S is given by p.
        Only used for V.type='fermion' or 'boson'
        """
        if self.type not in ['fermion','boson']:
            raise TypeError("Weights_of_S is only implemented for fermions")

        GS=LinGroup([len(p)])
        List_indices = []
        for i,x in enumerate(p) :
            List_indices += x*[i]
        
        chiTot=Weight(GS,as_list_of_list=[List_indices])
        for i,l in enumerate(unique_combinations(List_indices, self.nb_part)):
            chi = Weight(GS,as_list_of_list=[l],idx=i)
            chi.mult=prod([binomial(chiTot.as_vector[i],chi.as_vector[i]) for i in range(len(p))])
            yield chi


    #@staticmethod         
    def rhoEij(self,alpha : Root) -> matrix:
        """
        Return the matrix rho_V(E_alpha).
        Has to work when i==j.
        """
        M = matrix(QQ[I],self.dim,self.dim)
        if self.type == 'kron' : 
            Gred=LinGroup(self.G[:alpha.k]+self.G[alpha.k+1:])
            Vred=Representation(Gred,'kron')
            for w in Vred.all_weights:
                wj = Weight(self.G,as_list=list(w.as_list[:alpha.k])+[alpha.j]+list(w.as_list[alpha.k:]))
                idj=wj.idx(self)
                if alpha.i == alpha.j :
                    M[idj,idj]=1
                else :    
                    wi = Weight(self.G,as_list=list(w.as_list[:alpha.k])+[alpha.i]+list(w.as_list[alpha.k:]))
                    idi=wi.idx(self)
                    M[idi,idj]=1
        else : # Case Fermion and Boson
            
            if self.type == 'fermion' :
                shiftrank=1
            else :
                shiftrank=0
            Vred=Representation(LinGroup([self.G[0]-shiftrank]),self.type,self.nb_part-1)
            for w in Vred.all_weights: 
                L1=[s for s in w.as_list_of_list[0] if s<alpha.j]                
                L2=[s+shiftrank for s in w.as_list_of_list[0] if s>=alpha.j]
                lj=L1+[alpha.j]+L2  # we insert j
                if alpha.i == alpha.j:
                        wj  = Weight(self.G,as_list_of_list=[lj])
                        idj = wj.idx(self)
                        M[idj,idj]=1
                elif self.type == 'boson' or alpha.i not in lj : # Otherwise E_ij v =0
                    wj  = Weight(self.G,as_list_of_list=[lj])
                    idj = wj.idx(self)
                    li=L1+[alpha.i]+L2  # we insert i
                    li.sort()
                    wi = Weight(self.G,as_list_of_list=[li])
                    idi=wi.idx(self)
                    if self.type == 'fermion' :
                        M[idi,idj]=(-1)**(len(L1)-li.index(alpha.i))
                    else :
                        M[idi,idj]=lj.count(alpha.j)
                        
        return M

    @cached_property
    def actionK(self) -> dict[Root, Matrix]: #TODO : typer dictionnaire Root -> matrix
        """
        The list of matrices rho_V(xi) for xi in the bases of K.
        """
        L={}
        for beta in Root.all_of_B(self.G) : 
            i,j = beta.i,beta.j
            if i == j :
                L[beta]=I*self.rhoEij(beta)
            else :
                A=self.rhoEij(beta)
                B=self.rhoEij(beta.opposite)
                L[beta]=A-B
                L[beta.opposite]=I*(A+B)
        return(L)

    def action_op_el(self,alpha: Root, v: Vector) -> Vector:
        """
        Action of E_alpha on the Vector v in V
        The basis of V where v is expressed is indexed by V.all_weights

        Examples: TODO
        >>> from cone import *
        >>> d = Dimension((2, 3))
        >>> v = Vector(QQ, d.dimV)
        >>> v[:] = range(1, d.dimV + 1)
        >>> v
        (1, 2, 3, 4, 5, 6)
        >>> alpha = Root(k=1, i=0, j=2)
        >>> V.action_op_el(alpha, v)
        (3, 0, 0, 6, 0, 0)

        >>> d = Dimension((2, 3, 4))
        >>> v = Vector(QQ, d.dimV)
        >>> v[:] = range(1, d.dimV + 1)
        >>> alpha = Root(k=1, i=0, j=2)
        >>> V.action_op_el(alpha, v)
        (9, 10, 11, 12, 0, 0, 0, 0, 0, 0, 0, 0, 21, 22, 23, 24, 0, 0, 0, 0, 0, 0, 0, 0)
        """
        assert len(v) == self.dim
        from .rings import vector
        vp = vector(v.base_ring(),self.dim)
        if self.type == 'kron' :
            # Generate the weights with alpha.i and alpha.j in position alpha.k
            Gred=LinGroup(self.G[:alpha.k]+self.G[alpha.k+1:])
            Vred=Representation(Gred,'kron')
            for w in Vred.all_weights:
                wj = Weight(self.G,as_list=list(w.as_list[:alpha.k])+[alpha.j]+list(w.as_list[alpha.k:]))
                wi = Weight(self.G,as_list=list(w.as_list[:alpha.k])+[alpha.i]+list(w.as_list[alpha.k:]))
                vp[wi.idx(self)] = v[wj.idx(self)]
                
        else :
            # List of weights with j by inserting j from smaller
            Vred=Representation(LinGroup([self.G[0]-1]),self.type,self.nb_part-1)
            for w in Vred.all_weights: 
                L1=[s for s in w.as_list_of_list[0] if s<alpha.j]
                #posj=len(L1)
                L2=[s+1 for s in w.as_list_of_list[0] if s>=alpha.j]
                lj=L1+[alpha.j]+L2  # we insert j
                if self.type == 'boson' or alpha.i not in lj : # Otherwise E_ij v =0
                    wj  = Weight(self.G,as_list_of_list=[lj])
                    idj = wj.idx(self)  
                    if alpha.i==alpha.j:
                        M[idj,idj]=1
                    else :    
                        li=L1+[alpha.i]+L2  # we insert i
                        li.sort()
                        #posi=li.index(i)
                        wi = Weight(self.G,as_list_of_list=[li])
                        idi=wi.idx(self)
                        if self.type == 'fermion' :
                            vp[wi.idx(self)] = (-1)**(len(L1)-li.index(alpha.i))*v[wj.idx(self)]
                        else :
                            vp[wi.idx(self)] = lj.count(alpha.j)*v[wj.idx(self)]
                            
        return vp


    @cached_property
    def Q(self) -> "Ring":
        from .rings import QQ
        return QQ
    
    @cached_property
    def QI(self) -> "Ring":
        from .rings import QQ, I
        return QQ[I]

    @cached_property
    def QZ(self) -> "PolynomialRingForWeights":
        from .rings import PolynomialRingForWeights, QQ
        return PolynomialRingForWeights(QQ, "z")
    
    @cached_property
    def QV(self) -> "PolynomialRingForWeights":
        from .rings import PolynomialRingForWeights, QQ
        from .weight import Weight
        return PolynomialRingForWeights(QQ, weights=self.all_weights)
    
    @cached_property
    def QV2(self) -> "PolynomialRingForWeights":
        from .rings import PolynomialRingForWeights, QQ
        from .weight import Weight
        return PolynomialRingForWeights(QQ,
            weights=self.all_weights,
            seed=('va', 'vb'),
        )

    @cached_property
    def QIV(self) -> "PolynomialRingForWeights":
        from .rings import PolynomialRingForWeights, QQ, I
        from .weight import Weight
        return PolynomialRingForWeights(
            QQ[I], 
            weights=self.all_weights,
            seed=('vr', 'vi')
        )





    

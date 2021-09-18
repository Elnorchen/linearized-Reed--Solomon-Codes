from sage.coding.linear_code import AbstractLinearCode
from sage.coding.relative_finite_field_extension import *
from sage.rings.polynomial.skew_polynomial_element import *
import random
from sage.structure.sage_object import SageObject
from sage.coding.encoder import Encoder
from sage.coding.decoder import Decoder, DecodingError
import numpy as np



class LinearizedRSCode(AbstractLinearCode):
    r"""
    Class of the linearized Reed--Solomon codes $C_{LRS}[n,k]$
    
    INPUT:
         -``ground_field`` -- A finite field $\\GF{q}$ of a prime power order $q$
         
         -``extension_field`` -- A finite field $\\GF{q^m}$ which is an extension field of degree $m$ of $\\GF{q}$.
         
         -``length`` -- The length of the LRS Code, i.e., length $(n)$ should be less than or equal to $(m)$.
         
         -``dimension`` -- The dimesnion of the Gabidulin Code, i.e., dimension $(k)$ should be less than or equal to the length $(n)$.
         
         -``skew_polynomial_ring`` --the skew polynomial ring of the LRS code, i.e., skew_polynomial_ring $Fqm[x;\sigma]$.
         
         -``block_number`` -- number of blocks of the sum-rank metric. In other words, the sum-rank metric divides the code matrix into $l$ blocks
       
    EXAMPLE:
        sage: load('linearized_RS.py')
        sage: Fqm.<aa> = GF(8)
        sage: Fq.<a> = GF(2)
        sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
        sage: n=3
        sage: k=2
        sage: l=1
        sage: q=Fq.order()
        sage: p=Fq.characteristic()
        sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
        sage: L.<x>=Fqm['x',Frob]
        sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
        sage: C
        Linearized Reed--Solomon Code LRS[3,2] over Finite Field in aa of size 2^3
         
    
    
    
    """
    
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, ground_field, extension_field, skew_polynomial_ring, length, dimension, block_number):
        Fq=ground_field
        Fqm=extension_field
        FE=RelativeFiniteFieldExtension(Fqm, Fq)
        m=FE.extension_degree()
        super(LinearizedRSCode, self).__init__(Fqm, length, "GeneratormatrixEncoder", "WelchBerlekampDecoder")
        self._dimension=dimension
        self._length=length
        self._ground_field=Fq
        self._extension_field=Fqm
        self._extension_degree=m
        self._extension=FE
        self._skew_polynomial_ring=skew_polynomial_ring
        self._block_number=block_number
        self._Frob=skew_polynomial_ring.twist_map()
        self._independent_sets=0
        self._relative_extension=FE
        
    def length(self):
        r"""
        Returns the length n of self
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(8)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.length()
            
            3
        
        
        """
        return self._length
    
    def dimension(self):
        r"""
        Return the dimension k of self
        
        EXAMPLES:
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(8)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.dimension()
            
            2
        
        """
        return self._dimension
    
    def minimum_sum_rank_distance(self):
        r"""
        Return the minimum distance d of self
        
        EXAMPLES:
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(8)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.minimum_sum_rank_distance()
            
            2
        
        """
        return (self._length-self._dimension+1)
    
    def skew_polynomial_ring(self):
        r"""
        Return the skew polynomial ring Fqm[x;sigma] of self
        
        EXAMPLES:
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(8)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.skew_polynomial_ring()
            
            Skew Polynomial Ring in x over Finite Field in aa of size 2^3 twisted by aa |--> aa^2
        
        """
        
        return self._skew_polynomial_ring
    
    def extension_field(self):
        r"""
        Returns extension_field of self
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^4)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=4
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.extesion_field()
            
            Finite Field in aa of size 2^4
        """
        return self._extension_field
        
    
    def cardinality(self):
        r"""
        Return the cardinality of self
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.cardinality()
            
            64
        
        """
        m = self._extension_degree
        k = self._dimension
        q = self._ground_field.order()
        C = q**(m*k)
        return C
    
    def number_of_blocks(self):
        r"""
        Return the number of blocks l of self
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.number_of_blocks()
            
            1
        
        """
        
        return (self._block_number)
    
    def block_length(self):
        r"""
        Return the block length of self
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.block_length()
            
            3
        
        """
        return (self._length//self._block_number)
    
    def Frobenuis_endomorphism(self):
        r"""
        Return the Frobenuis endomorphism sigma of self
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.Frobenuis_endomorphism()
            
            Frobenius endomorphism aa |--> aa^2 on Finite Field in aa of size 2^3
        
        """
        return self._Frob
    
    
    def ground_field(self):
        r"""
        Return the ground field Fq of self
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.ground_field()
            
            Finite Field of size 2
        
        """
        return self._ground_field
    
    def extension_field(self):
        r"""
        Return the extension field Fqm of self
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.extension_field()
            
            Finite Field in aa of size 2^3
        
        """
        return self._extension_field
    
    def extension_degree(self):
        r"""
        Return the extension degree m of self
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: C.extension_degree()
            
            3
        
        """
        return self._extension_degree
    
    def map_into_ground_field(self, vector):
        r"""
        Change an input vector v1 defined in extension field into a vector in ground field
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: v=[aa^2,aa^0,aa]  ## Note that here we use aa^0 to replace 1. Otherwise there will be an error.
            sage: C.map_into_ground_field(v)
            
            [0 0 1]
            [1 0 0]
            [0 1 0]
        
        """
        
        FE = self._relative_extension
        vector_length = len(vector)
        for i in range(vector_length):
            if vector[i] not in self.extension_field():
                raise ValueError("input should be an element of %s" % self.extension_field() )
        vector_ground = []
        for i in range(vector_length):
            vector_ground.append(FE.relative_field_representation(vector[i]))
        #print(vector_ground)
        return matrix(vector_ground)
    
    def sum_rank_weight(self, codeword_vector):
        r"""
        Return the sum rank weight of a vector in Fqm
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: v=[aa^2,aa^0,aa]  ## Note that here we use aa^0 to replace 1. Otherwise there will be an error.
            sage: C.map_into_ground_field(v)
            
            3
        
        """
        
        codeword_matrix=self.map_into_ground_field(codeword_vector)
        l=self._block_number
        s=self.block_length()
        sub_blocks=[]
        sub_ranks=[]
        for i in range(l):
            sub_blocks.append(codeword_matrix[i*l:(i+1)*s])
            sub_ranks.append(codeword_matrix[i*l:(i+1)*s].rank())
        sum_rank_wt=sum(sub_ranks)
        return sum_rank_wt
        
    def sum_rank_distance(self, vector1, vector2):
        r"""
        Return the sum rank distance of two vectors in Fqm
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: v1=[aa^2,aa^0,aa]  ## Note that here we use aa^0 to replace 1. Otherwise there will be an error.
            sage: v2=[aa^0,aa^0,aa^0]
            sage: C.map_into_ground_field(v1,v2)
            
            2

        
        """
        codeword_matrix1=self.map_into_ground_field(vector1)
        codeword_matrix2=self.map_into_ground_field(vector2)
        cw_matrix=codeword_matrix1-codeword_matrix2
        #print("vector matrix=")
        #print(cw_matrix)      #these two lines output the difference matrix of two input vectors
        l=self._block_number
        s=self.block_length()
        sub_blocks1=[]
        sub_ranks1=[]
        for i in range(l):
            sub_blocks1.append(cw_matrix[i*l:(i+1)*s])
            sub_ranks1.append(cw_matrix[i*l:(i+1)*s].rank())
        sum_rank_wt=sum(sub_ranks1)
        return sum_rank_wt
       
    def get_independent_sets(self):
        r"""
        Return the independent sets Beta of self. Note that the initial value of independent sets is 0. The self.generator_matrix should be implemented    first to do the assignment
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: E=C.encoder("GeneratormatrixEncoder")
            sage: b=[[aa^0,aa,aa^1]]     #note that here we have two pairs of brackets
            sage: G=C.generator_matrix(b)
            sage: C.get_independent_sets()
            
            [[1, aa, aa^2]]
        
        """
        return self._independent_sets
    
    def Dab(self,element1, element2):
        r"""
        Return the result of linear operator Dab=sigma(b)a
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: b=aa  
            sage: c=aa^2
            sage: C.Dab(b,c)
            
            aa^2 + aa + 1   ## note that aa^5=aa^2+aa+1

        
        """
        dd=self._Frob(element2)*element1 ##this is the linear operator D_a(b)=sigma(b)a in thesis
        
        return dd
    
    def Dab_power(self,element1,element2,power):
        r"""
        Return the ith power of linear operator (Dab)^i=sigma^i(b)sigma^{i-1}(a)...sigma(a)a
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: b=aa  
            sage: c=aa^2
            sage: C.Dab_power(b,c,3)
            
            aa^2  ## note that aa^23=aa^2

        
        """
        
        a=element1     ########this is the power of operator D_a(b). D^n_a(b)=sigma^n(b)sigma^(n-1)(a)....sigma(a)a
        b=element2
        if power> 0:
            
            for i in range(power):
                b=self.Dab(a,b)
        elif power==0:
            b=element2
        else:
            raise ValueError("power should be non-negative" )
        
        return b
            
        
    
    def Independent_Sets(self):
        r"""
        Assign the independent sets of self with random normal basis
        """
        m=self._extension_degree
        FE=self._extension
        sets=[]
        l=self._block_number
        n=self._length
        s=n//l      ###here we must use // instead of / in order to get an integer
        
        FieldBasis=FE.absolute_field_basis()
        if l>1:
            
            for i in range(l):
                sets.append(random.sample(FieldBasis,s))    ##every element in sets is a independent set in Fqm
        else:
            sets=FieldBasis
            
        self._independent_sets=sets
        return sets
    
    def generator_matrix(self, independent_sets=None):
        r"""
        Creat and return the generator matrix of the LRS code
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: E=C.encoder("GeneratormatrixEncoder")
            sage: b=[[aa^0,aa,aa^1]]     #note that here we have two brackets
            sage: G=C.generator_matrix(b)   ## We strongly advice to assign the independent sets
            sage: G
            
            [   1   aa   aa]
            [   1 aa^2 aa^2]

        
        """
        ## independent_sets shall be written as [[1,alpha,alpha^2],[alpha,alpha^3,1]]
        l=self._block_number
        k=self._dimension
        n=self._length
        s=n//l
        Fqm=self._extension_field
        alpha=Fqm.primitive_element()
        if independent_sets is None:
            independent_sets=self.Independent_Sets()
        else:
            self._independent_sets=independent_sets
        beta=[]
        
        #print(independent_sets)
        G=np.empty(shape=[k,n],dtype=sage.rings.finite_rings.element_givaro.FiniteField_givaroElement) ###这里G 的
        
        if s>1 and l>1:
            for i1 in range(l):
                beta=independent_sets[i1]
                
                for i2 in range(k):
                    for i3 in range(s):
                    
                        G[i2,i1*s+i3]=self.Dab_power(alpha**(i1),beta[i3],i2)
                        
        elif l==1:
            beta=independent_sets
            for i1 in range(n):
                for i2 in range(k):
                    G[i2,i1]=self.Dab_power(1,beta[0][i1],i2)
                    
        elif s==1:
             for i1 in range(l):
                beta=independent_sets[i1]
                
                for i2 in range(k):
                    
                    
                        G[i2,i1]=self.Dab_power(alpha**(i1),beta[0],i2)
                        
                        
                        
        else:
            raise ValueError("block number l should be a positive integer and l<=n" )
            
            
        G=matrix(G)    
        return G
    
    def precomputing_polynomials(self,basis,r):## r=[r1,r2,...,rk]
        r"""
        calculate and return the minimum skew polynomial Fk and the Lagrange interpolating polynomial G
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^3)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: E=C.encoder("GeneratormatrixEncoder")
            sage: b=[aa^0,aa,aa^1]     #note that here we have only one pair of brackets
            sage: r=[aa,aa^2,aa]
            sage: C.precomputing_polynomials(b,r)   ## We strongly advice to assign the independent sets
            
            
            [aa*x, x^2 + (aa^2 + aa + 1)*x + aa^2 + aa] ##the first element of this list is G, the second element is Fk
        
        """
        k=self._dimension
        FE=self._extension
        Fqm=self._extension_field
        L=self._skew_polynomial_ring
        Fro=self._Frob
        
        
        
        F=L([basis[0]*(-1),1])
        G=L([r[0]])
        #print("r=",r)
        #print("basis=",basis)
        for i in range(k-1):
            
          
            G=G+(r[i+1]-G.right_mod(x-basis[i+1]))*(F.right_mod(x-basis[i+1]))**(-1)*F
            
            F=(x-Fro(F.right_mod(x-basis[i+1]))*F.right_mod(x-basis[i+1])**(-1)*basis[i+1])*F
            
        P=[]
        P.append(G)
        P.append(F)
        return P   ####P[0] is G，P[1] is F
    
    def code_space(self):
        r"""
        Return the code space of self
        """
        return VectorSpace(self.extension_field(),self.length())
    
    def LRS_TO_SKEW(self,received_word):
        r"""
        Translate the received word into skew Reed--Solomon code case, and translate the independent sets into P-basis. Note that the independent sets of self should be assigned before you run this self.LRS_TO_SKEW method. So far, this method is only used in the decoding algorithm
        
        EXAMPLE::
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(8)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=3
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: b=[[aa^0,aa,aa^1]]     #note that here we have two pair of brackets
            sage: r=[aa,aa^2,aa]
            sage: C.generator_matrix(b)  ##we run the self.generator_marix first in order to assign the independent sets.
            sage: C.LRS_TO_SKEW(r)
            
            [[aa, aa, 1], [1, aa, aa]]  #the first element is received word in skew Reed--Somon code form. The second element is the P-basis
        
        """
        y=np.mat(received_word)
        y=y.tolist()
        y=y[0]
        
        
        n=self.length()
        k=self.dimension()
        Fqm=self.extension_field()
        sets=self.get_independent_sets()## sets is beta
        r=[0]*n  
        a=[0]*n## a is bij
        l=self.number_of_blocks()
        s=self.block_length()
        aa=Fqm.primitive_element()
        #print("sets=",sets)
        for i in range(l):
            for j in range(s):
                #print("i=",i,"j=",j)
                if y[i*s+j]==1:
                    y[i*s+j]=aa**0###in order to avoid the error since "1" is regarded as an integer, not an element in Finite Field
                
                if sets[i][j]==1:
                    sets[i][j]=aa**0
                r[i*s+j]=y[i*s+j]*(sets[i][j])**(-1)
                a[i*s+j]=C.Dab(aa**(i),sets[i][j])*(sets[i][j])**(-1)
                
        
        p=[]
        
        p.append(r)
        p.append(a)
        return  p          
            
    #def Lagrange_interpolation_poly(self,equivX,equivY):##here equivX=(x1,x2,x3,....,xn),equivY=(y1,y2,....,yn) and (x1,y1),(x2,y2),..,(xn,yn) are the interpolating points
        #points=[]
        #n=self._length
        #for i in range(n):
            #points.append((equivX[i],equivY[i]))
            
        #L=self._skew_polynomial_ring
        #pol=L.lagrange_polynomial(points)
        #return pol#
        
    #def minimum_skew_polynomial(self,elements):##elements should be a linearly independent set [alpha1,alpha2,...,alphas] and alpha in Fqm. Here we use the lclm to calculate it
        #e=elements
        #L=self._skew_polynomial_ring
        #k=self._dimension
        #l=L([elements[0]*(-1),1])
        #for i in range(k-1):
            #l=l.left_lcm(x-elements[i+1])
        #return l
        
    def _repr_(self):
        return "Linearized Reed--Solomon Code LRS[%s,%s] over %s" % (self.length(), self.dimension(), self.extension_field())
    
    def _latex_(self):
        return "\\textnormal{Linearized Reed--Solomon Code LRS}[%s,%s] \\textnormal{over }%s" % (self.length(), self.dimension(), self.extension_field())
    
    def __eq__(self,other):
        return isinstance(other, LinearizedRSCode)\
                    and self.length() == other.length()\
                    and self.dimension() == other.dimension()\
            and self.ground_field() == other.ground_field()\
            and self.extension_field() == other.extension_field()\
    
   
#############################################################    Encoder  part     ####################################################################   
    
    
class GeneratormatrixEncoder(Encoder):
        
        def __init__(self,code):
            super(GeneratormatrixEncoder, self).__init__(code) 
         
        
        def __eq__(self):
            return isinstance(other, Encoder1)\
                    and self.code() == other.code()\
        
        def _latex_(self):
            return "\\textnormal{Generator matrix based encoder for }%s" % self.code()._latex_()
        
        def _repr_(self):
            return "Generator matrix based encoder for %s" % self.code()
        
        def encode(self,Generator,information):
            r"""
              encode using generator matrix
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^4)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=4
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: E=C.encoder("GeneratormatrixEncoder")
            sage: b=[[aa^0,aa,aa^3,aa^2]]     #note that here we have two brackets
            sage: G=C.generator_matrix(b)   
            sage: E=C.encoder("GeneratormatrixEncoder")
            sage: c=[aa,aa^2]        ##c is the information vector
            sage: y=E.encode(G,c)
            sage: y
            
            the information vector is [  aa aa^2]
            the generator matrix is 
            [          1          aa        aa^3        aa^2]
            [          1        aa^2 aa^3 + aa^2      aa + 1]

            [    aa^2 + aa aa^2 + aa + 1     aa^2 + aa          aa^2] ##this is the codeword
        
        """
            C=self.code()
            k=C.dimension()
            i=information
            G=Generator
            if len(i)==k:
                i=matrix(i)
                print("the information vector is",i)
                
                print("the generator matrix is ")
                print(G)
                codeword=i*G
                
            else:
                raise ValueError("the length of information vector i should be equal to code dimension k" )
                
            return codeword
            
        
        
######################################################  Decoder part   ##############################################################################        
        
class WelchBerlekampDecoder(Decoder):
    
    def __init__(self,code):
        if not isinstance(code, LinearizedRSCode):
            raise ValueError("The code given as input is not a Linearized Reed--Solomon code")
        super(WelchBerlekampDecoder, self).__init__(code, code.code_space(),"GeneratormatrixEncoder")###这里相对Musab的code有改动。先这么写占位
        
        
    def __eq__(self):
        return (isinstance(other, WelchBerlekampDecoder)
            and self.code() == other.code())
    
    
    def _repr_(self):
        return "decoder for %s" % self.code()
    
    def _latex_(self):
        return "\\textnormal{decoder for }%s" % self.code()._latex_()
    
    def decode_to_message(self,received_word):##received_word and basis in list form
        r"""
        Decode the received word into the information vector. Note that at least an transmitting error with sum-rank weight 1 should be added to the received word. Otherwise, there might be errors in running the decoder.
        
        EXAMPLES::
        
            sage: load('linearized_RS.py')
            sage: Fqm.<aa> = GF(2^4)
            sage: Fq.<a> = GF(2)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: n=4
            sage: k=2
            sage: l=1
            sage: q=Fq.order()
            sage: p=Fq.characteristic()
            sage: Frob = Fqm.frobenius_endomorphism(log(q,p))
            sage: L.<x>=Fqm['x',Frob]
            sage: C=LinearizedRSCode(Fq,Fqm,L,n,k,l)
            sage: E=C.encoder("GeneratormatrixEncoder")
            sage: b=[[aa^0,aa,aa^3,aa^2]]     #note that here we have two brackets
            sage: G=C.generator_matrix(b)   
            sage: E=C.encoder("GeneratormatrixEncoder")
            sage: D=C.decoder("WelchBerlekampDecoder")
            sage: c=[aa,aa^2]        ##c is the information vector
            sage: y=E.encode(G,c)    ## here y is [aa^2 + aa, aa^2 + aa + 1, aa^2 + aa,aa^2] according to the example in the self.encode()
            sage: i=D.decode_to_message([aa^2+aa,aa^2+aa,aa^2+aa,aa^2]) ##put an error [0,1,0,0] to the received word
            sage: i
            
            the information vector is [  aa aa^2]
            the generator matrix is 
            [          1          aa        aa^3        aa^2]
            [          1        aa^2 aa^3 + aa^2      aa + 1]
            Fk= x^2 + (aa^2 + aa + 1)*x + aa^2 + aa
            G= (aa + 1)*x + aa^2 + 1
            Q= (aa^2 + aa + 1)*x^2 + aa^3 + aa + 1
            L= (aa^3 + aa^2)*x + aa^3 + aa^2
            Y= aa^2*x + aa
            R= 0
               
            [aa, aa^2]  ##the information vector getting from decoding

        
        """
        C=self.code()
        n=C.length()
        k=C.dimension()
        Fqm=C.extension_field()
        sets=C.get_independent_sets()## sets is beta
        r1=received_word
        
        l=C.number_of_blocks()
        s=C.block_length()
        aa=Fqm.primitive_element()
        ss=C.LRS_TO_SKEW(r1)
        r=ss[0]
        a=ss[1] ## a is bij
        
        LL=C.skew_polynomial_ring()
        Fro=C.Frobenuis_endomorphism()
        #for j in range(n):
            #r[j]=r[j]*independent_sets1[j]
      
        GG=C.precomputing_polynomials(a,r)
        Gk=GG[0]
        Fk=GG[1]
        L=LL([0])
        #L1=LL([0,1])
        L1=LL([1])
        Q=Fk
        Q1=Gk
        print("Fk=",Fk)
        print("G=",Gk)
        #print("L,l1,Q,Q1=",L,L1,Q,Q1)
        for i in range(k,n):
            s=L.right_mod(x-(Frob(r[i])*(r[i])**(-1)*a[i]))*r[i]-Q.right_mod(x-a[i])
            s1=L1.right_mod(x-(Frob(r[i])*(r[i])**(-1)*a[i]))*r[i]-Q1.right_mod(x-a[i])
            #print("s,s1=",s,s1)
            L2=L
            L=L1
            L1=L2
            Q2=Q
            Q=Q1
            Q1=Q2
            s2=s
            s=s1
            s1=s2
            
                
            L1=s*L1-s1*L
            Q1=s*Q1-s1*Q
                
            
            if s!=0:
                L=(x-Frob(s)*s**(-1)*a[i])*L
                Q=(x-Frob(s)*s**(-1)*a[i])*Q
            #print("第次计算结果Q，L，Q1，L1=",[Q,L,Q1,L1])
        print("Q=",Q1)
        print('L=',L1)
        Y,R=Q1.left_quo_rem(L1)###还差一个Euclidean division
        print("Y=",Y)
        print("R=",R)
        
        I=Y.coefficients(sparse=False)##this should be the information vector
        
        return I
                

           
                
        
            
            
        
       
    
    
    
#################################   registration of default encoders and decoders #############################################################3    
    
    
LinearizedRSCode._registered_encoders["GeneratormatrixEncoder"] = GeneratormatrixEncoder
LinearizedRSCode._registered_decoders["WelchBerlekampDecoder"] = WelchBerlekampDecoder

     
            
            
   

    
        
  

    
    
    
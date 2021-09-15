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
        -``
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
        return self._length
    
    def dimension(self):
        return self._dimension
    
    def minimum_sum_rank_distance(self):
        return (self._length-self._dimension+1)
    
    def skew_polynomial_ring(self):
        return self._skew_polynomial_ring
    
    def extension_field(self):
        r"""
        Returns extension_field of ``self``
        
        EXAMPLES:
        
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
        m = self._extension_degree
        k = self._dimension
        q = self._ground_field.order()
        C = q**(m*k)
        return C
    
    def number_of_blocks(self):
        return (self._block_number)
    
    def block_length(self):
        return (self._length//self._block_number)
    
    def Frobenuis_endomorphism(self):
        return self._Frob
    
    
    def ground_field(self):
        return self._ground_field
    
    def extension_field(self):
        return self._extension_field
    
    def extension_degree(self):
        return self._extension_degree
    
    def map_into_ground_field(self, vector):
        FE = self._relative_extension
        vector_length = len(vector)
        for i in range(vector_length):
            if vector[i] not in self.extension_field():
                raise ValueError("input should be an element of %s" % self.extension_field() )
        vector_ground = []
        for i in range(vector_length):
            vector_ground.append(FE.relative_field_representation(vector[i]))
        print(vector_ground)
        return matrix(vector_ground)
    
    def sum_rank_weight(self, codeword_vector):
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
        codeword_matrix1=self.map_into_ground_field(vector1)
        codeword_matrix2=self.map_into_ground_field(vector2)
        cw_matrix=codeword_matrix1-codeword_matrix2
        print("cwmatrix=",cw_matrix)
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
        return self._independent_sets
    
    def Dab(self,element1, element2):
        dd=self._Frob(element2)*element1 ##this is the linear operator D_a(b)=sigma(b)a in thesis
        
        return dd
    
    def Dab_power(self,element1,element2,power):
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
            ###here there is an error about the type of s. Maybe I can use another function to do this. But I don't know which function can I use.
        self._independent_sets=sets
        return sets
    
    def generator_matrix(self, independent_sets=None):
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
        return P   ####P[0]是G，P[1]是F
    
    def code_space(self):
        return VectorSpace(self.extension_field(),self.length())
    
    def LRS_TO_SKEW(self,received_word):
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

     
            
            
   

    
        
  

    
    
    
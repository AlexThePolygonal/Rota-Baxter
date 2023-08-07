import sage
import sage.all as sg
from timeit import default_timer as timer
from sage.structure.sage_object import SageObject
from sage.all import singular as singular


_ = singular.LIB("primdec.lib")



def get_entries(mat : sg.Matrix) -> list[sg.RingElement]:
    '''
    Helper func
    '''
    res = []
    for row in mat:
        for entry in row:
            res.append(entry)
    return res

def is_absolutely_prime(sage_QQ_ideal) -> bool:
    '''
    Uses singular to check whether an ideal is prime over the algebraic closure
    '''
    S = singular.ideal(sage_QQ_ideal).absPrimdecGTZ()
    singular.setring(S)
    return singular.execute("size(absolute_primes)") == '1'

def singular_locus(sage_ideal):
    pass


class OperatorComponent(SageObject):
    ideal = None
    is_absolutely_prime = True
    is_smooth = True
    
    def __init__(self, ideal):
        self.ideal = ideal
        self.is_absolutely_prime = is_absolutely_prime(ideal)
            # self.is_smooth


class RBClassifier(SageObject):
    '''
        Inherit from this class to classify RB-operators on an algebra
        The child class must implement the following:
          `self.gens` --- the generators of the algebra
          `self.to_vector` --- decompose an element using the generators
          `self.to_matrix` --- the inverse to `self.to_vector`
    '''

    def __init__(self, weight) -> None:
        '''
        Generate the equations defining the RB operator moduli space
        '''
        self.dim = len(self.gens)
        self.weight = weight
        
        self.rota_vars = sg.SR.var(','.join(f'r{i+1}{j+1}' for i in range(self.dim) for j in range(self.dim)))
        self.vars = self.rota_vars + weight.variables()
        
        self.rota_mat = sg.Matrix(sg.SR, self.dim, self.dim)
        for i in range(self.dim):
            for j in range(self.dim):
                self.rota_mat[i,j] = self.rota_vars[i*self.dim + j]

        R = self.rota
        raw_eqns = []
        for x in self.gens:
            for y in self.gens:
                cur_eqns = R(x)*R(y) - (R(R(x)*y + x*R(y)) - weight*R(x*y))
                raw_eqns.append(get_entries(cur_eqns))
        self.raw_eqns = get_entries(raw_eqns)

    def rota(self, alg_elem):
        '''
        The RB operator
        '''
        return self.to_matrix(self.rota_mat.T * self.to_vector(alg_elem))
    
    def pretty_print(self, ideal):
        pretty_vars = sg.SR.var('a,b,c,d,e,f,g,h,e,f,g,h,i,j,k,l')
        gens = sorted(ideal.gens(), key=lambda x: x.degree())
        lins = [gen for gen in gens if gen.degree() == 1]
        substs = dict()
        occurs = dict()

        for gen in gens:
            for var_ in gen.variables():
                if occurs.get(var_) == None:
                    occurs[var_] = 1
                else:
                    occurs[var_] += 1

        for lin in lins:
            to_elim = None
            for v in lin.variables():
                if occurs[v] == 1:
                    to_elim = v
            if to_elim is not None:
                coeff = lin.monomial_coefficient(to_elim)
                lin_ = lin / coeff - to_elim
                substs[to_elim] = lin_
                gens.remove(lin)
        
        raise NotImplementedError()



        

    def classify(self):
        self.ring = sage.rings.polynomial.polynomial_ring_constructor.PolynomialRing(sg.QQ, self.vars)
        self.irrelevant_ideal = self.ring.ideal(self.ring.gens())

        self.ideal = self.ring.ideal(self.raw_eqns)
        self.QQ_comps = self.ideal.minimal_associated_primes()
        sorted(self.QQ_comps, key=sg.dimension)

        for comp in self.QQ_comps:
            if not is_absolutely_prime(comp):
                raise NotImplementedError("Non-rational compoments are not supported")
        
        self.intersection_locus = self.ring.ideal([1])
        for i in range(len(self.QQ_comps)):
            for j in range(i+1, len(self.QQ_comps)):
                a = self.QQ_comps[i]
                b = self.QQ_comps[j]
                self.intersection_locus.intersection((a+b).radical())
        


        

        

    

    
    

class MnClassifier(RBClassifier):
    '''
    Classify RB-operators on the matrix algebra Mₙ(ℂ)
    '''
    def __init__(self, dim, weight) -> None:
        self.alg = sg.MatrixSpace(sg.SR, dim)
        self.gens = self.alg.basis().values()
        super().__init__(weight)

    def to_vector(self, alg_elem):
        return sg.Matrix(get_entries(alg_elem)).T
    
    def to_matrix(self, alg_elem_in_basis):
        return sum([self.gens[i] * alg_elem_in_basis[i][0] for i in range(len(self.gens))])



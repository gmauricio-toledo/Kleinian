import numpy as np

class Transformation:

    def __init__(self,matrix,affine_chart=True):
        self.matrix = matrix
        assert matrix.shape[0] == matrix.shape[1]
        self.dim = matrix.shape[0]
        self.__get_det()
        self.unitary_representation = (1/self.det)*self.matrix # Representación con det=1
        self.trace = np.sum([self.unitary_representation[i,i] for i in range(self.dim)])
        self.__get_type()
        if self.dim == 2:
            self.affine_chart = affine_chart

    def __get_type(self):
        '''
        This method determines whether the element is parabolic, elliptic or loxodromic
        '''
        if self.dim == 2:
            if self.trace**2 == 4:
                self.type = "parabolic"
            elif self.trace.imag == 0:
                if self.trace**2 < 4:
                    self.type = "elliptic"
                if self.trace**2 > 4:
                    self.type = "hyperbolic"
            else:
                self.type = "loxodromic"
        elif self.dim == 3:
            print("To be added...")
        else:
            print("Not available")

    def __get_det(self):
        '''
        Computes the "standard" form of the matrix
        '''
        self.det = np.linalg.det(self.matrix)
        # print(f"El determinante es {self.det}")
        return

    def inverse(self):
        if self.det != 0:
            return Transformation(np.linalg.inv(self.matrix))
        else:
            print("The matrix is singular")

    # ---- Operation overriding ----

    def __add__(self, other):
        return Transformation(self.unitary_representation + other.unitary_representation)
    
    def __sub__(self, other):
        return Transformation(self.unitary_representation - other.unitary_representation)
    
    def __mul__(self, other):
        return Transformation(self.unitary_representation * other.unitary_representation)
    
    # def __truediv__(self, other):
    #     return Transformation(self.unitary_representation - other.unitary_representation)

    # ---- Printing ----

    # ---- Make the class callable ----

    def __call__(self,z):
        if self.dim == 2:
            if self.affine_chart:  # Nos restringimos a la carta afin de CP1, ignoramos los infinitos
                gz = np.dot(self.unitary_representation,np.array([z,1]))
                if gz[1] != 0:
                    return gz[0]/gz[1]
                else:
                    return None
            else:
                pass
            
        elif self.dim == 3:
            print("To be added...")
        else:
            print("Not available")

    def __repr__(self):
        return f"Transformation:\n{self.matrix}"

    def __str__(self):
        return f"Transformation:\n{self.matrix}"

class KleinianGroup:

    def __init__(self, generadores):
        self.generadores = list(generadores)
        # assert todos tienen la misma dimension
        self.dim = self.generadores[0].dim 

    def is_discrete(self):
        if len(self.generadores)==2:
            print("Using Jørgensen\'s inequality...")
            [g1, g2] = self.generadores
            conj = g1*(g2*(g1.inverse()*g2.inverse()))
            test = np.absolute(g1.trace**2-4) + np.absolute(conj.trace - 2)
            if test >= 1:
                print("The group is discrete")

    def plot(self):
        pass

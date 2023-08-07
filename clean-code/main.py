from rblib import *


# classifier = MnClassifier(2, sg.SR(sg.Integer(0)))
# # classifier = MnClassifier(2, sg.SR.var('w'))
# classifier.classify()

# print(f"number of classes is {len(classifier.QQ_comps)}")
x,y,z,w = sg.SR.var('x,y,z,w')
Ring = sage.rings.polynomial.polynomial_ring_constructor.PolynomialRing(sg.QQ, (x,y,z,w))
ideal = Ring.ideal(x*y - w**2, z)

I = sg.macaulay2(ideal)


J = I.minimalPresentation()

singlocJ = J.singularLocus().ideal()

print(J._sage_())





loadPackage "Cremona"

P = QQ[r11,r12,r13,r14,r21,r22,r23,r24,r31,r32,r33,r34,r41,r42,r43,r44]

IInter = ideal(r41 + r44,r34 - r42,r33 + r44,r31 + r42,r24 - r43,r22 + r44,r21 + r43,r14 + r44,r13 + r43,r12 + r42,r11 - r44,r43^2 - r23*r44,r42*r43 + r44^2,r32*r43 + r42*r44,r42^2 - r32*r44,r23*r42 + r43*r44,r23*r32 - r44^2)


IA = ideal(-r42*r43 + r41*r44,-r42^2 + r32*r44,-r41*r42 + r32*r43,-r41*r42 + r31*r44,-r41^2 + r31*r43,-r32*r41 + r31*r42,r24*r42 + r44^2,r24*r41 + r43*r44,r24*r32 + r42*r44,r24*r31 + r42*r43,-r24*r43 + r23*r44,r23*r42 + r43*r44,r23*r41 + r43^2,r23*r32 + r42*r43,r23*r31 + r41*r43,r34 - r42,r33 - r41,r22 + r44,r21 + r43,r14 + r44,r13 + r43,r12 + r42,r11 + r41)
IB = ideal(-r42*r43 + r41*r44,r34*r43 + r44^2,r34*r41 + r42*r44,-r34*r42 + r32*r44,r32*r43 + r42*r44,r32*r41 + r42^2,-r43^2 + r23*r44,r23*r42 - r41*r43,r23*r34 + r43*r44,r23*r32 + r42*r43,-r41*r43 + r21*r44,-r23*r41 + r21*r43,-r41^2 + r21*r42,r21*r34 + r42*r43,r21*r32 + r41*r42,r33 + r44,r31 + r42,r24 - r43,r22 - r41,r14 + r44,r13 + r43,r12 + r42,r11 + r41)
IC1 = ideal(r24*r42 + r44^2,r43,r41,r34,r33 - r44,r32,r31 - r42,r23,r22 + r44,r21,r14,r13 - r24,r12,r11 + r44)
IC2 = ideal(-r41^2 + r21*r42,r44,r43,r34 + r42,r33 + r41,r32,r31,r24,r23,r22 - r41,r14 + r41,r13 + r21,r12,r11)
IC3 = ideal(r34*r43 + r44^2,r42,r41,r33 + r44,r32,r31,r24,r23,r22 - r44,r21 - r43,r14,r13,r12 - r34,r11 + r44)
IC4 = ideal(-r41^2 + r31*r43,r44,r42,r34,r33 - r41,r32,r24 + r43,r23,r22 + r41,r21,r14 + r41,r13,r12 + r31,r11)

-- Components related by transposition
C3 = P / IC3
C1 = P / IC1
isBirational rationalMap map(C1, C3, {r11, r13, r12, r14, r31, r33, r32, r34, r21, r23, r22, r24, r41, r43, r42, r44})

-- Components related by transposition
C2 = P / IC2
C4 = P / IC4
isBirational rationalMap map(C4, C2, {r11, r13, r12, r14, r31, r33, r32, r34, r21, r23, r22, r24, r41, r43, r42, r44})

-- Components related by transposition
A = P / IA
B = P / IB
isBirational rationalMap map(B, A, {r11, r13, r12, r14, r31, r33, r32, r34, r21, r23, r22, r24, r41, r43, r42, r44})


-- fixed point of transposition
Inter = P / IInter
isBirational rationalMap map(Inter, Inter, {r11, r13, r12, r14, r31, r33, r32, r34, r21, r23, r22, r24, r41, r43, r42, r44})



mInter = minimalPresentation (P / IInter)

mA = minimalPresentation (P / IA)
mB = minimalPresentation (P / IB)

-- All smooth
-- radical ideal singularLocus Spec mA
-- radical ideal singularLocus Spec mB
-- radical ideal singularLocus Spec mInter


import numpy as np
import matplotlib.pyplot as plt

def intersection_two_segs(seg1_m, seg1_b, seg2_m, seg2_b):
    a = np.array([
        [-seg1_m, 1],
        [-seg2_m, 1]
    ])
    b = np.array([seg1_b, seg2_b])
    return np.linalg.solve(a, b)

A = np.array([0, 0])
B = np.array([2, -2])
C = np.array([10, 0])
d = 0.5

# segment BA
dist_BA = np.linalg.norm(A - B)
unit_BA = (A - B) / dist_BA

# segment BC
dist_BC = np.linalg.norm(C - B)
unit_BC = (C - B) / dist_BC

# bisecting vector
bisect = (unit_BA + unit_BC)/2
unit_bisect = bisect / np.linalg.norm(bisect)

# new frame axes
y_prime_hat = -1 * unit_bisect
x_prime_hat = np.array([-unit_bisect[1], unit_bisect[0]])
theta = np.arctan2(x_prime_hat[1], x_prime_hat[0])

# origin of (x_prime_hat, y_prime_hat) coordinate frame
O = B + d*unit_bisect

# projection matrix - (x_hat, y_hat) -> (x_prime_hat, y_prime_hat)
P = np.array([
    [np.cos(theta), -np.sin(theta), O[0]],
    [np.sin(theta), np.cos(theta), O[1]],
    [0, 0, 1]
])
P_inv = np.linalg.inv(P)
# P = np.array([
#     [np.cos(theta), -np.sin(theta), O[0] * np.cos(theta) - O[1] * np.sin(theta)],
#     [np.sin(theta), np.cos(theta), O[0] * np.sin(theta) + O[1] * np.cos(theta)],
#     [0, 0, 1]
# ])

# find intersection point of x-prime-hat and line segment AB
AB_m = (B[1]-A[1])/(B[0]-A[0])
AB_b = (-AB_m * A[0]) + A[1]

x_prime_m = x_prime_hat[1] / x_prime_hat[0]
x_prime_b = (-x_prime_m * O[0]) + O[1]

I1 = intersection_two_segs(AB_m, AB_b, x_prime_m, x_prime_b)

# transform I1 into the x_prime_hat, y_prime_hat space
I1_homo = np.array([I1[0], I1[1], 1])
I1_prime = np.matmul(P_inv, I1_homo)[:2]
print(f"this is O: {O}")
print(f"this is P: {P}")

print(f"this is i1: {I1}")
print(f"this is i1_prime: {I1_prime}")

# transform the unit vector BA into x_prime_hat, y_prime_hat space and use it to calculate I1_prime_m
# unit_BA_homo = np.array([unit_BA[0], unit_BA[1], 1])
# unit_BA_prime = np.matmul(P_inv, unit_BA_homo)[:2]
# I1_prime_m = unit_BA_prime[1] / unit_BA_prime[0]
# print(f"this is i1_prime_m: {I1_prime_m}")
B_homo = np.array([B[0], B[1], 1])
B_prime = np.matmul(P_inv, B_homo)[:2]

I1_prime_m = (B_prime[1] - I1_prime[1]) / (B_prime[0] - I1_prime[0])
print(f"this is i1_prime_m: {I1_prime_m}")

# find a, c parameters of parabola in x_prime_hat, y_prime_hat space
a = I1_prime_m / (2 * I1_prime[0])
c = I1_prime[1] - (a * (I1_prime[0] ** 2))

# plot line segments AB, BC
plt.plot([A[0], B[0]], [A[1], B[1]])
plt.plot([B[0], C[0]], [B[1], C[1]])

# plot origin O of x_prime_hat, y_prime_hat space
plt.scatter(O[0], O[1], s=4)

# plot x_prime_hat, y_prime_hat axes
y_prime_hat_pt = O + y_prime_hat
plt.plot([O[0], y_prime_hat_pt[0]], [O[1], y_prime_hat_pt[1]])

x_prime_hat_pt = O + 3 * x_prime_hat
x_prime_hat_pt2 = O - x_prime_hat
plt.plot([O[0], x_prime_hat_pt[0]], [O[1], x_prime_hat_pt[1]])
plt.plot([O[0], x_prime_hat_pt2[0]], [O[1], x_prime_hat_pt2[1]])

# plot parabola
x_upper = np.sqrt(-c / a) + 0.1
x_lower = -x_upper
print(x_lower, x_upper)
xs_prime = np.linspace(x_lower, x_upper, 500)
print(xs_prime)
ys_prime = [a*(x**2)+c for x in xs_prime]

print(P_inv)
xs = []
ys = []
for i in range(len(xs_prime)):
    pt_homo = np.array([xs_prime[i], ys_prime[i], 1])
    pt = np.matmul(P, pt_homo)[:2]
    xs.append(pt[0])
    ys.append(pt[1])
plt.scatter(xs, ys, s=1)

plt.scatter(I1[0], I1[1], s=25)
plt.show()




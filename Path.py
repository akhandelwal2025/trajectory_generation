import math
import numpy as np
import constants
import matplotlib.pyplot as plt

class LinearJointSegment:
    def __init__(self, start, end):
        # start = (s1, s_dot1), end = (s2, s_dot2)
        self.start = start
        self.end = end
        self.m, self.offset, self.b = self.calc_m_offset_b()
        self.start_s = self.start[0]

    def calc_m_offset_b(self):
        m = (self.end[1] - self.start[1]) / (self.end[0] - self.start[0])
        offset = self.start[0]
        b = self.start[1]
        return m, offset, b

    def f(self, s):
        return self.m * (s - self.offset) + self.b

    def f_prime(self, s):
        return self.m

    def f_prime2(self, s):
        return 0

class ParabolicJointSegment:
    def __init__(self, start, intersection, end):
        self.start = np.array(list(start))
        self.intersection = np.array(list(intersection))
        self.end = np.array(list(end))
        self.a, self.c = self.generate_parabola()

        # record some important info for use in f(s), f_prime(s), f_prime2(s)
        self.start_s = self.I1[0]
        self.start_s_prime = self.I1_prime[0]
        self.end_s = self.I2[0]
        self.end_s_prime = self.I2_prime[0]
    
    def intersection_two_segs(self, seg1_m, seg1_b, seg2_m, seg2_b):
        a = np.array([
            [-seg1_m, 1],
            [-seg2_m, 1]
        ])
        b = np.array([seg1_b, seg2_b])
        return np.linalg.solve(a, b)
        
    def generate_parabola(self):
        # segment BA
        dist_BA = np.linalg.norm(self.start - self.intersection)
        unit_BA = (self.start - self.intersection) / dist_BA

        # segment BC
        dist_BC = np.linalg.norm(self.end - self.intersection)
        unit_BC = (self.end - self.intersection) / dist_BC

        # bisecting vector
        bisect = (unit_BA + unit_BC)/2
        unit_bisect = bisect / np.linalg.norm(bisect)

        # new frame axes
        y_prime_hat = -1 * unit_bisect
        x_prime_hat = np.array([-unit_bisect[1], unit_bisect[0]])
        theta = np.arctan2(x_prime_hat[1], x_prime_hat[0])

        # find intersection point of x-prime-hat and line segment AB
        AB_m = (self.intersection[1]-self.start[1])/(self.intersection[0]-self.start[0])
        AB_b = (-AB_m * self.start[0]) + self.start[1]

        BC_m = (self.end[1]-self.intersection[1])/(self.end[0]-self.intersection[0])
        BC_b = (-BC_m * self.intersection[0]) + self.intersection[1]

        # origin of (x_prime_hat, y_prime_hat) coordinate frame
        print(f"slope: {BC_m - AB_m}")
        slope_diff = BC_m - AB_m
        d = constants.blend_radius
        if slope_diff == 120:
            d = 0.05
        O = self.intersection + d*unit_bisect

        # projection matrix - (x_hat, y_hat) -> (x_prime_hat, y_prime_hat)
        self.P = np.array([
            [np.cos(theta), -np.sin(theta), O[0]],
            [np.sin(theta), np.cos(theta), O[1]],
            [0, 0, 1]
        ])
        self.P_inv = np.linalg.inv(self.P)

        
       
        x_prime_m = x_prime_hat[1] / x_prime_hat[0]
        x_prime_b = (-x_prime_m * O[0]) + O[1]

        self.I1 = self.intersection_two_segs(AB_m, AB_b, x_prime_m, x_prime_b)
        self.I2 = self.intersection_two_segs(BC_m, BC_b, x_prime_m, x_prime_b)

        # transform I1 into the x_prime_hat, y_prime_hat space
        I1_homo = np.array([self.I1[0], self.I1[1], 1])
        self.I1_prime = np.matmul(self.P_inv, I1_homo)[:2]

        # transform I2 into the x_prime_hat, y_prime_hat space
        I2_homo = np.array([self.I2[0], self.I2[1], 1])
        self.I2_prime = np.matmul(self.P_inv, I2_homo)[:2]

        # transform the unit vector BA into x_prime_hat, y_prime_hat space and use it to calculate I1_prime_m
        B_homo = np.array([self.intersection[0], self.intersection[1], 1])
        B_prime = np.matmul(self.P_inv, B_homo)[:2]

        I1_prime_m = (B_prime[1] - self.I1_prime[1]) / (B_prime[0] - self.I1_prime[0])

        # find a, c parameters of parabola in x_prime_hat, y_prime_hat space
        a = I1_prime_m / (2 * self.I1_prime[0])
        c = self.I1_prime[1] - (a * (self.I1_prime[0] ** 2))

        # """
        #     PLOTTING STUFF
        # """
        # A = self.start
        # B = self.intersection
        # C = self.end
        # plt.plot([A[0], B[0]], [A[1], B[1]])
        # plt.plot([B[0], C[0]], [B[1], C[1]])
        # print(A, B, C)
        # # plot origin O of x_prime_hat, y_prime_hat space
        # plt.scatter(O[0], O[1], s=4)

        # # plot x_prime_hat, y_prime_hat axes
        # y_prime_hat_pt = O + y_prime_hat
        # # plt.plot([O[0], y_prime_hat_pt[0]], [O[1], y_prime_hat_pt[1]])

        # x_prime_hat_pt = O + x_prime_hat
        # # x_prime_hat_pt2 = O - x_prime_hat
        # # plt.plot([O[0], x_prime_hat_pt[0]], [O[1], x_prime_hat_pt[1]])
        # # plt.plot([O[0], x_prime_hat_pt2[0]], [O[1], x_prime_hat_pt2[1]])

        # # plot parabola
        # x_upper = np.sqrt(-c / a) + 0.1
        # x_lower = -x_upper
        # print(x_lower, x_upper)
        # xs_prime = np.linspace(x_lower, x_upper, 500)
        # print(xs_prime)
        # ys_prime = [a*(x**2)+c for x in xs_prime]

        # print(self.P_inv)
        # xs = []
        # ys = []
        # for i in range(len(xs_prime)):
        #     pt_homo = np.array([xs_prime[i], ys_prime[i], 1])
        #     pt = np.matmul(self.P, pt_homo)[:2]
        #     xs.append(pt[0])
        #     ys.append(pt[1])
        # plt.scatter(xs, ys, s=1)

        # plt.scatter(self.I1[0], self.I1[1], s=25)
        # plt.show()
        return a, c

    """
        given: s -> [s1, s2] | output: x -> [x1, x2]
        linear mapping formula: x1 + ((s - s1)/(s2 - s1) * (x2 - x1)) 
    """
    def map_s_to_x_prime(self, s, x1, x2, s1, s2):
        return x1 + ((s-s1)/(s2-s1) * (x2-x1))
    
    def f(self, s):
        x_prime = self.map_s_to_x_prime(s, self.start_s_prime, self.end_s_prime, self.start_s, self.end_s)
        y_prime = self.a * (x_prime ** 2) + self.c
        y_prime_homo = np.array([x_prime, y_prime, 1])
        return np.matmul(self.P, y_prime_homo)[1]
    
    def f_prime(self, s):
        x_prime = self.map_s_to_x_prime(s, self.start_s_prime, self.end_s_prime, self.start_s, self.end_s)
        y_prime = 2 * self.a * x_prime
        y_prime_homo = np.array([x_prime, y_prime, 1]) 
        return np.matmul(self.P, y_prime_homo)[1]
    
    def f_prime2(self, s):
        return 2 * self.a
    
class JointPath:
    def __init__(self, joint_positions, vel_constraint, accel_constraint):
        # joint_positions = list of thetas for a single joint
        # ex. joint1 = [30, 90, 60, 60]
        # vel_constraint = q-dot-max for the joint
        self.joint_positions = joint_positions
        self.s_interval = 1/len(self.joint_positions)
        self.q_dot_max = vel_constraint
        self.q_dot2_max = accel_constraint

        self.linear_segments = self.construct_linear_segments()
        self.segments, self.starting_s = self.construct_parabolas()

        """
            plot segments
        """
        # all_s = []
        # all_theta = []
        # print(self.starting_s)
        # for i in range(len(self.segments)-1):
        #     start_s = self.starting_s[i]
        #     end_s = self.starting_s[i+1]
        #     seg = self.segments[i]
        #     s = np.linspace(start_s, end_s, 100)
        #     theta = [seg.f(s_i) for s_i in s]
        #     all_s.extend(s)
        #     all_theta.extend(theta)
        #     plt.scatter(all_s, all_theta, s=2)
        #     # plt.scatter(s, theta, s=2)
        #     plt.show()

    def construct_parabolas(self):
        starting_s = []
        segments = []
        import pdb
        for i in range(len(self.linear_segments)-1):
            seg1 = self.linear_segments[i]
            seg2 = self.linear_segments[i+1]

            # if consecutive segments have the same slope, no need for parabolic blend
            if seg1.m == seg2.m:
                segments.append(seg1)
                starting_s.append(seg1.start_s)
                continue
            
            parabolic_seg = ParabolicJointSegment(seg1.start, seg1.end, seg2.end)

            # add first linear segment and new parabolic blend. second linear segment will get added at the end
            segments.append(seg1)
            segments.append(parabolic_seg)

            # record the starting positions of each of the segments
            starting_s.append(seg1.start_s)
            starting_s.append(parabolic_seg.start_s)

            # set the starting point of the next linear segment to the intersection point with the parabola
            seg2.start = parabolic_seg.I2
            seg2.start_s = seg2.start[0]

        # add last remaining segment
        segments.append(self.linear_segments[-1])
        starting_s.append(self.linear_segments[-1].start_s)
        return segments, starting_s

    def construct_linear_segments(self):
        start = (0, 0) # (s, theta)
        segments = []
        for i, end_theta in enumerate(self.joint_positions):
            end_s = (i+1) * self.s_interval
            end = (end_s, end_theta)
            segments.append(LinearJointSegment(start, end))
            start = end
        return segments
    
    def f(self, s):
        seg_idx = 0
        for i, segment in enumerate(self.segments):
            if s >= segment.start_s:
                seg_idx = i 
        return self.segments[seg_idx].f(s)

    def f_prime(self, s):
        seg_idx = 0
        for i, segment in enumerate(self.segments):
            if s >= segment.start_s:
                seg_idx = i
        return self.segments[seg_idx].f_prime(s)

    def f_prime2(self, s):
        seg_idx = 0
        for i, segment in enumerate(self.segments):
            if s >= segment.start_s:
                seg_idx = i
        return self.segments[seg_idx].f_prime2(s)
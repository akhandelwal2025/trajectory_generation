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
        self.start_s, self.start_theta = self.I1
        self.start_x, self.start_y = self.I1_prime
        self.end_s, self.end_theta = self.I2
        self.end_x, self.end_y = self.I2_prime
        
        # sample a bunch of x-y values and project back to s-theta
        self.s_theta_pts = self.sample_s_theta_pts()
    
    def sample_s_theta_pts(self):
        x = np.linspace(self.start_x, self.end_x, constants.s_theta_sampling_frequency)
        y = [self.a * (x_i ** 2) + self.c for x_i in x]
        x_y_pts = list(zip(x, y))
        s_theta_pts = []
        for x, y in x_y_pts:
            pt_homo = np.array([x, y, 1])
            s_theta_pt = np.matmul(self.P, pt_homo)
            s_theta_pts.append((s_theta_pt[0], s_theta_pt[1]))
        # print(len(s_theta_pts))
        return s_theta_pts
    
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

        """
            find O - origin of the new frames
            1. find points at center of segment AB, BC (center_BA, center_BC)
            2. construct line segment between center_BA, center_BC
            3. find where this intersects bisecting vector segment - this is O
        """
        center_BA = (self.start + self.intersection)/2
        center_BC = (self.intersection + self.end)/2

        center_m = (center_BC[1] - center_BA[1]) / (center_BC[0] - center_BA[0])
        center_b = (-center_m * center_BA[0]) + center_BA[1]

        bisect_m = unit_bisect[1]/unit_bisect[0]
        bisect_b = (-bisect_m * self.intersection[0]) + self.intersection[1]

        if bisect_m == float('inf') or bisect_m == float('-inf'):
            x = self.intersection[0]
            y = center_m * x + center_b
            O = (x, y)
        else:
            O = self.intersection_two_segs(center_m, center_b, bisect_m, bisect_b)
        # plt.plot([self.start[0], self.intersection[0]], [self.start[1], self.intersection[1]])
        # plt.plot([self.intersection[0], self.end[0]], [self.intersection[1], self.end[1]])
        # plt.plot([self.intersection[0], self.intersection[0] + unit_bisect[0]], [self.intersection[1], self.intersection[1] + unit_bisect[1]])
        # plt.scatter(center_BA[0], center_BA[1])
        # plt.scatter(center_BC[0], center_BC[1])
        # plt.scatter(O[0], O[1])
        # plt.xlim(0, 1)
        # plt.ylim(0, 180)
        # plt.show()
        # new frame axes
        y_prime_hat = -1 * unit_bisect
        x_prime_hat = np.array([-unit_bisect[1], unit_bisect[0]])
        theta = np.arctan2(x_prime_hat[1], x_prime_hat[0])

        # find intersection point of x-prime-hat and line segment AB
        AB_m = (self.intersection[1]-self.start[1])/(self.intersection[0]-self.start[0])
        AB_b = (-AB_m * self.start[0]) + self.start[1]

        BC_m = (self.end[1]-self.intersection[1])/(self.end[0]-self.intersection[0])
        BC_b = (-BC_m * self.intersection[0]) + self.intersection[1]

        x_prime_m = x_prime_hat[1] / x_prime_hat[0]
        x_prime_b = (-x_prime_m * O[0]) + O[1]

        self.I1 = self.intersection_two_segs(AB_m, AB_b, x_prime_m, x_prime_b)
        self.I2 = self.intersection_two_segs(BC_m, BC_b, x_prime_m, x_prime_b)

        # projection matrix - (x_hat, y_hat) -> (x_prime_hat, y_prime_hat)
        self.P = np.array([
            [np.cos(theta), -np.sin(theta), O[0]],
            [np.sin(theta), np.cos(theta), O[1]],
            [0, 0, 1]
        ])
        self.P_inv = np.linalg.inv(self.P)

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

        return a, c

    """
        given: s -> [s1, s2] | output: x -> [x1, x2]
        linear mapping formula: x1 + ((s - s1)/(s2 - s1) * (x2 - x1)) 
    """
    def map_s_to_x_prime(self, s, x1, x2, s1, s2):
        return x1 + ((s-s1)/(s2-s1) * (x2-x1))
    
    def calc_proj_s_theta(self, x_proj):
        y_proj = self.a * (x_proj ** 2) + self.c 
        proj_pt = np.matmul(self.P, np.array([x_proj, y_proj, 1]))
        return proj_pt[0], proj_pt[1]
    
    def f(self, s):
        # --------- IMPLEMENT RAYS HACKY SOLUTION ---------
        # if s == 0:
        #     return 0
        # delta_x = math.sqrt(abs(s/self.a))
        # x = self.start_x + delta_x
        # y = self.a * (x ** 2) + self.c
        # pt_homo = np.array([x, y, 1])
        # theta = np.matmul(self.P, pt_homo)[1]
        # return theta
        # --------- IMPLEMENT BINARY SEARCH THROUGH ALR SAMPLED POINTS ---------
        # print(s, self.start_s, self.end_s)
        s_sampled, theta_sampled = -1, -1
        high = len(self.s_theta_pts)-1
        low = 0
        while low <= high:
            mid = int((high+low)/2)
            s_sampled, theta_sampled = self.s_theta_pts[mid]  

            if s < s_sampled:
                high = mid-1
            elif s > s_sampled:
                low = mid+1
            else:
                return theta_sampled
        # print(s, self.s_theta_pts[-1])
        
        
        s_low, theta_low = self.s_theta_pts[high]
        if high == len(self.s_theta_pts)-1:
            s_high, theta_high = self.end_s, self.end_theta
        else:
            s_high, theta_high = self.s_theta_pts[high+1]
        return self.map_s_to_x_prime(s, theta_low, theta_high, s_low, s_high)

        # --------- IMPLEMENT BINARY SEARCH : KINDA WORKS, BUT NOT REALLY BECAUSE ITS NOT ABLE TO SOLVE AT HIGH FREQUENCY ---------
        # high = self.end_x
        # low = self.start_x
        # mid = (high + low)/2 
        # s_proj, theta_proj = self.calc_proj_s_theta(mid)
        # while abs(s - s_proj) > constants.s_epsilon:
        #     if s > s_proj:
        #         low = mid
        #     else:
        #         high = mid
        #     mid = (high + low)/2
        #     s_proj, theta_proj = self.calc_proj_s_theta(mid)
        # return theta_proj

        # --------- TRANSFORMING S INTO X WITH PROJECTION MATRIX : DOESN'T WORK BECAUSE DON'T KNOW WHAT THE THETA-VALUE SHOULD BE OF THE (S, THETA) PAIR --------- 
        # x = np.matmul(self.P_inv, np.array([s, 0, 1]))[0]
        # y = self.a * (x ** 2) + self.c
        # y_homo = np.array([x, y, 1])
        # return np.matmul(self.P, y_homo)[1]
    
        # --------- ESTABLISHING LINEAR MAPPING FROM S -> X : DOESN'T WORK BECAUSE S -> X IS NOT A LINEAR MAPPING --------- 
        # x_prime = self.map_s_to_x_prime(s, self.start_s_prime, self.end_s_prime, self.start_s, self.end_s)
        # y_prime = self.a * (x_prime ** 2) + self.c
        # y_prime_homo = np.array([x_prime, y_prime, 1])
        # return np.matmul(self.P, y_prime_homo)[1]
    
    def f_prime(self, s):
        # find (s - epsilon, theta1) and (s + epsilon, theta2)
        # print(s, self.start_s, self.end_s)
        s1 = max(self.start_s, s - constants.parabolic_epsilon)
        theta1 = self.f(s1)
        s2 = min(self.end_s, s + constants.parabolic_epsilon)
        theta2 = self.f(s2)

        return (theta2-theta1)/(s2-s1) 

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
        self.segments = self.construct_parabolas()

    def construct_parabolas(self):
        segments = []
        import pdb
        for i in range(len(self.linear_segments)-1):
            seg1 = self.linear_segments[i]
            seg2 = self.linear_segments[i+1]

            # if consecutive segments have the same slope, no need for parabolic blend
            if seg1.m == seg2.m:
                segments.append(seg1)
                continue
            
            parabolic_seg = ParabolicJointSegment(seg1.start, seg1.end, seg2.end)

            # add first linear segment and new parabolic blend. second linear segment will get added at the end
            segments.append(seg1)
            segments.append(parabolic_seg)

            # set the starting point of the next linear segment to the intersection point with the parabola
            seg2.start = parabolic_seg.I2
            seg2.start_s = seg2.start[0]

        # add last remaining segment
        segments.append(self.linear_segments[-1])
        return segments

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
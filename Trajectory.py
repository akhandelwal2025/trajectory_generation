from Path import *
from Limit import *
import constants
import matplotlib.pyplot as plt
import time
import numpy as np

class Trajectory:
    def __init__(self, waypoints, vel_constraints, accel_constraints):
        """
            waypoints = [[q11, q12, q13, ..., q1m], 
                         [q21, q22, q23, ..., q2m], 
                         ....,
                         [qn1, qn2, qn3, ..., qnm]]
            n = number of waypoints
            m = number of joints
        """
        self.waypoints = waypoints
        self.vel_constraints = vel_constraints
        self.accel_constraints = accel_constraints
        self.joint_paths = self.generate_joint_paths()
        
        self.vel_limit_curve = VelocityLimitCurve(self.joint_paths)
        self.accel_limit_curve = AccelerationLimitCurve(self.joint_paths)

        # self.limit_curve = LimitCurve(self.joint_paths)
        # self.inflection_pts = self.limit_curve.find_inflection_points()

        # prev timestep - used for collision-checking
        self.prev_s = 0
        self.prev_s_dot = 0
        self.prev_time = 0

        # curr timestep
        self.curr_s = 0
        self.curr_s_dot = 0
        self.curr_time = 0

        # points where you switch from max accel to max decel
        self.switching_pts = []

        # store (s, s_dot, time) pairs for viz purposes
        self.forward_path = []
        self.backward_path = []
        self.intersection_pts = []

        # matplotlib stuff
        self.fig, self.axs = plt.subplots(3, 2, figsize=(20, 8))

    """
        Definitions:
        - inflection point: s-value where the limit curve bound changes

        Meat of the algorithm (Arbys)
        1. Start at (s=0, s-dot=0)
        2. Integrate forward w/ max accel until you intersect the limit curve
            a. keep note of the s, s-dot value upon intersection
        3. Find the next inflection point (i) and integrate backward from there
            a. if you hit the limit curve at an s-value greater than the 
               s-value where the forward integration hit the limit curve, then join
               the forward and backward paths with the straight line of the limit curve
            b. if you intersect with the forward integration path, then store the intersection
               point as a switching point (max accel to max decel)
        4. Set the inflection point (i) as a switching point as well (max decel to zero)
        5. Continue along velocity limit curve until next inflection point (i+1) and set that
           as switching point
        6. Repeat steps 2-5 until done condition  
    """
    def generate_trajectory(self):
        # step 2
        while not self.done() and \
              not self.find_limit_curve_collisions(accel=False, forward=True) and \
              not self.find_limit_curve_collisions(accel=True, forward=True):
            self.integrate_forward()
        print(self.forward_path[-2])
        self.curr_s, self.curr_s_dot, _ = self.forward_path[-3]
        while not self.done() and \
              not self.find_limit_curve_collisions(accel=False, forward=False) and \
              not self.find_limit_curve_collisions(accel=True, forward=False) and \
              not self.curr_s == 0:
            print(self.curr_s, self.curr_s_dot)
            self.integrate_backward()

        # extract (position, velocity, time) for all s_interval + switching points
        self.final_path = self.forward_path
        # self.output_trajectory()

    def output_trajectory(self):
        s_interval = 1/len(self.waypoints)/constants.discretization
        all_waypoint_s = np.arange(s_interval, 1.000000001, s_interval) # produces list of [s_interval, 2*s_interval, ... , 1]
        print(all_waypoint_s)
        num_waypoints = len(all_waypoint_s)
        num_joints = len(self.joint_paths)
        positions = np.empty((num_waypoints, num_joints))
        velocities = np.empty((num_waypoints, num_joints))
        times = np.empty((num_waypoints,))
        for s_i, s in enumerate(all_waypoint_s):
            # populate positions
            for joint_i, joint_path in enumerate(self.joint_paths):
                positions[s_i, joint_i] = np.deg2rad(joint_path.f(s))

            s_dot, time = self.extract_s_dot_time(s)

            # populate velocities - theta-dot = q_prime(s) * s_dot
            for joint_i, joint_path in enumerate(self.joint_paths):            
                velocities[s_i, joint_i] = np.deg2rad(joint_path.f_prime(s) * s_dot)
            
            # populate time for the waypoint
            times[s_i] = time
        
        self.plot_output_positions(positions, times)
        self.plot_velocity(velocities, times)
        np.set_printoptions(precision=3, suppress=False, floatmode='fixed')
        print(positions)
        print(velocities)
        print(times)
    
    """
        self.final_path is a discrete path of (s, s_dot, time)
        given an arbitrary s (such as one of the s_intervals), we need to interpolate what the corresponding s_dot and time should be
    """
    def extract_s_dot_time(self, s):
        for i, path_elem in enumerate(self.final_path):
            if path_elem[0] == s:
                return path_elem[1], path_elem[2]
            
            """
                keep going until find the first element in the path that is greater than the s value, call this path_elem
                find the prev_path elem by calling self.final_path[i-1]
                calculate the total difference in s_dot, time between prev_path_elem and path_elem
                calculate the percent along the segment formed by prev_path_elem and path_elem that s is
                return the interpolated values of s_dot, time along that line segment
            """
            if path_elem[0] > s:
                prev_path_elem = self.final_path[i-1]
                diff_s_dot = path_elem[1] - prev_path_elem[1]
                diff_time = path_elem[2] - prev_path_elem[2]
                percent_along_seg = (s - prev_path_elem[0])/(path_elem[0] - prev_path_elem[0])
                to_ret_s_dot = prev_path_elem[1] + diff_s_dot * percent_along_seg
                to_ret_time = prev_path_elem[2] + diff_time * percent_along_seg
                return to_ret_s_dot, to_ret_time
        
        raise Exception(f"s value is greater than reached by path. s: {s}")
    
    def fix_back_path_timestamps(self):
        """
            the way the backward path is set up, time starts at 0 at the inflection pt and then continues to be negative until it intersects with forward path
            thus, self.backward_path[0][2] (i.e. the collision point) is the most negative timestamp
            therefore, offset should be the timestamp where the forward path ended (self.forward_path[-1][2]) + the most negative time of the backward path (-1 * self.backward_path[0][2])
        """
        offset = self.forward_path[-1][2] + (-1*self.backward_path[0][2]) 
        for i, path_elem in enumerate(self.backward_path):
            s, s_dot, time = path_elem
            self.backward_path[i] = (s, s_dot, time + offset)

    def find_next_inflection_pt(self):
        self.curr_s = min(1.0, self.curr_s)
        self.curr_time = 0
        for inflection_pt in self.inflection_pts:
            if self.curr_s <= inflection_pt[0]:
                self.curr_s = inflection_pt[0]
                self.curr_s_dot = inflection_pt[1]
                return

        raise Exception(f"shouldn't be here: {self.curr_s} > 1")
        
    """
        given:
            - (self.curr_s, self.curr_s_dot)
            - (self.prev_s, self.prev_s_dot)
            - path: [(s1, s_dot1), (s2, s_dot2), ...]
        goal:
            - find if the line segment formed by (self.curr_s, self.curr_s_dot) and (self.prev_s, self.prev_s_dot)
              intersects with any line segment on the path
            - if it does, then calculate the collision point and return that
        procedure:
            - calculate slope + bias for line segment formed by (self.curr_s, self.curr_s_dot) and (self.prev_s, self.prev_s_dot)
            - for each pair of consecutive points in the path:
                - calculate slop + bias for the line segment
                - find intersection point with curr-prev line segment calculated earlier
                - if intersection point is within rectangle formed by both line segments, then its a collision point and return that
    """
    def find_path_collision(self, path):
        A = (self.curr_s, self.curr_s_dot)
        B = (self.prev_s, self.prev_s_dot)
        for i in range(1, len(path)):
            C = path[i-1]
            D = path[i]
            if self.intersect(A, B, C, D):
                self.forward_path = self.forward_path[0:i]
                return True
        return None
    
    def ccw(self, A,B,C):
        return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

    # Return true if line segments AB and CD intersect
    def intersect(self, A,B,C,D):
        return self.ccw(A,C,D) != self.ccw(B,C,D) and self.ccw(A,B,C) != self.ccw(A,B,D)

    # pt1 = (s1, s_dot1), pt2 = (s1, s_dot2)
    def calc_a_b_c(self, pt1, pt2):
        s1, s_dot1 = pt1
        s2, s_dot2 = pt2
        m = (s_dot2 - s_dot1) / (s2 - s1) # (y2-y1)/(x2-x1)
        """
            point-slope formula: y - y1 = m(x - x1)
            -> y = mx + (-mx1 - y1)
            -> -mx + y = (-mx1 - y1)
            -> mx - y = (mx1 + y1)
            -> a = m, b = -1, c = mx1 + y1
        """
        a = m
        b = -1
        c = m * s1 + s_dot1
        return a, b, c

    """
        two cases for collisions
            a. abs(s_dot - limit_curve(s)) <= epsilon
                i. do nothing - s_dot is within bounds for a collision
            b. s_dot > limit_curve(s) && abs(s_dot - limit_curve(s)) > epsilon
                i. need to go backward to find a precise collision point
    """
    def find_limit_curve_collisions(self, accel, forward):
        if accel:
            lim_curve_s_dot = self.accel_limit_curve.evaluate(self.curr_s)
        else:
            lim_curve_s_dot = self.vel_limit_curve.evaluate(self.curr_s)

        if self.curr_s_dot < lim_curve_s_dot and (lim_curve_s_dot - self.curr_s_dot) > constants.epsilon:
            if forward:
                self.forward_path.append((self.curr_s, self.curr_s_dot, self.curr_time))
            else:
                self.backward_path.append((self.curr_s, self.curr_s_dot, self.curr_time))
            return False
        
        # case a
        if abs(lim_curve_s_dot - self.curr_s_dot) <= constants.epsilon:
            if forward:
                self.forward_path.append((self.curr_s, self.curr_s_dot, self.curr_time))
            else:
                self.backward_path.append((self.curr_s, self.curr_s_dot, self.curr_time))
            return True

        # case b - find (s, s_dot) for precise collision point
        lower_timestep = 0
        upper_timestep = constants.timestep
        while self.curr_s_dot > lim_curve_s_dot and abs(lim_curve_s_dot - self.curr_s_dot) > constants.epsilon: 
            timestep = (lower_timestep + upper_timestep) / 2
            if forward:
                max_s_dot2 = self.calc_max_s_dot2(self.prev_s, self.prev_s_dot)
                self.curr_s_dot = self.prev_s_dot + max_s_dot2 * timestep # v_final = v_initial + a*delta_t
                self.curr_s = self.prev_s + (self.prev_s_dot * timestep) + (0.5 * max_s_dot2 * (timestep ** 2))# s_final = s_initial + v_initial*delta_t + 1/2*a*delta_t^2 
                self.curr_time = self.prev_time + timestep
            else:
                min_s_dot2 = self.calc_min_s_dot2(self.prev_s, self.prev_s_dot)
                self.curr_s_dot = self.prev_s_dot - min_s_dot2 * timestep # v_final = v_initial + a*delta_t
                self.curr_s = self.prev_s - (self.prev_s_dot * timestep) + (0.5 * min_s_dot2 * (timestep ** 2)) # s_final = s_initial + v_initial*delta_t + 1/2*a*delta_t^2 
                self.curr_time = self.prev_time - timestep

            # if resulting s_dot is above lim_curve_s_dot, then set upper_timestep = timestep
            if accel:
                lim_curve_s_dot = self.accel_limit_curve.evaluate(self.curr_s)
            else:
                lim_curve_s_dot = self.vel_limit_curve.evaluate(self.curr_s)

            if self.curr_s_dot > lim_curve_s_dot:
                upper_timestep = timestep
            else:
                lower_timestep = timestep
        if forward:
            self.forward_path.append((self.curr_s, self.curr_s_dot, self.curr_time))
        else:
            self.backward_path.append((self.curr_s, self.curr_s_dot, self.curr_time))
        return True

    def integrate_backward(self):
        # curr_s, curr_s_dot will be overwritten. need to store them in prev_s, prev_s_dot
        self.prev_s = self.curr_s
        self.prev_s_dot = self.curr_s_dot
        self.prev_time = self.curr_time

        min_s_dot2 = self.calc_min_s_dot2(self.prev_s, self.prev_s_dot)
        self.curr_s_dot = self.prev_s_dot - min_s_dot2 * constants.timestep # v_final = v_initial + a*delta_t
        self.curr_s = self.prev_s - (self.prev_s_dot * constants.timestep) + (0.5 * min_s_dot2 * (constants.timestep ** 2)) # s_final = s_initial + v_initial*delta_t + 1/2*a*delta_t^2 
        self.curr_time = self.prev_time - constants.timestep

    def integrate_forward(self):
        # curr_s, curr_s_dot, curr_time will be overwritten. need to store them in prev_s, prev_s_dot, prev_time
        self.prev_s = self.curr_s
        self.prev_s_dot = self.curr_s_dot
        self.prev_time = self.curr_time

        max_s_dot2 = self.calc_max_s_dot2(self.prev_s, self.prev_s_dot)
        print(max_s_dot2)
        self.curr_s_dot = self.prev_s_dot + max_s_dot2 * constants.timestep # v_final = v_initial + a*delta_t
        self.curr_s = self.prev_s + (self.prev_s_dot * constants.timestep) + (0.5 * max_s_dot2 * (constants.timestep ** 2)) # s_final = s_initial + v_initial*delta_t + 1/2*a*delta_t^2 
        self.curr_time = self.prev_time + constants.timestep

    # check if trajectory has reached the end 
    def done(self):
        return self.curr_s < 0 or self.curr_s > 1

    # take waypoints for each joint (i.e. joint angles) and convert them into a continuous path for the joint to follow
    # class JointPath generates piece-wise function that dictates the path of each joint through s-space
    def generate_joint_paths(self):
        all_joints = list(zip(*self.waypoints)) # groups all joint positions together
        assert(len(all_joints) == len(self.vel_constraints) == len(self.accel_constraints))
        joint_paths = []
        for i in range(len(all_joints)):
            joint_i_positions = all_joints[i]
            joint_i_vel_constraint = self.vel_constraints[i]
            joint_i_accel_constraint = self.accel_constraints[i]
            joint_paths.append(JointPath(joint_i_positions, joint_i_vel_constraint, joint_i_accel_constraint))
        return joint_paths
    
    # implementing eqn 22
    def calc_min_s_dot2(self, s, s_dot):
        min_s_dot2 = float('-inf')
        for joint_i in self.joint_paths:
            f_prime_s = joint_i.f_prime(s)
            f_prime2_s = joint_i.f_prime2(s)
            if f_prime_s != 0:
                min_s_dot2 = max(min_s_dot2, (-joint_i.q_dot2_max/abs(f_prime_s)) - ((f_prime2_s*(s_dot ** 2))/f_prime_s))
        return min_s_dot2

    # implementing eqn 23
    def calc_max_s_dot2(self, s, s_dot):
        max_s_dot2 = float('inf')
        for joint_i in self.joint_paths:
            f_prime_s = joint_i.f_prime(s)
            f_prime2_s = joint_i.f_prime2(s)
            if f_prime_s != 0:
                max_s_dot2 = min(max_s_dot2, (joint_i.q_dot2_max/abs(f_prime_s)) - ((f_prime2_s*(s_dot ** 2))/f_prime_s))
        return max_s_dot2

    def plot_segments(self):
        for i, joint_path in enumerate(self.joint_paths):
            s = np.linspace(0, 1, 50000)
            theta = [joint_path.f(i) for i in s]
            theta_prime = [joint_path.f_prime(i) for i in s]
            self.axs[1, 0].scatter(s, theta, s=2)
            self.axs[1, 1].scatter(s, theta_prime, s=2)
            # self.axs[1, 0].set_xlim(0, 1)
            # self.axs[1, 0].set_ylim(-180, 180)
            # self.axs[1, 1].set_xlim(0, 1)
            # self.axs[1, 1].set_ylim(-180, 180)
            
    def plot_output_positions(self, positions, times):
        for i in range(len(self.joint_paths)):
            joint_pos_waypoints = positions[:, i]
            self.axs[0, 0].scatter(times, joint_pos_waypoints, s=4)
    
    def plot_velocity(self, velocities, times):
        for i in range(len(self.joint_paths)):
            joint_velocity_waypoints = velocities[:, i] # velocities: [[
            self.axs[0, 1].scatter(times, joint_velocity_waypoints, s=4)
            
    def plot_limit_curve(self):
        # plot limit curve
        s = np.linspace(0, 1, 500)[:-1]
        s_dot_max_vel = [self.vel_limit_curve.evaluate(i) for i in s]
        s_dot_max_accel = [self.accel_limit_curve.evaluate(i) for i in s]
        self.axs[2, 1].scatter(s, s_dot_max_vel, s=2)
        self.axs[2, 1].scatter(s, s_dot_max_accel, s=2)
        self.axs[2, 1].set_xlim(0, 1)
        self.axs[2, 1].set_ylim(0, 1)

    def plot_forward_path(self):
        # plot forward path
        s, s_dot_max, _ = list(zip(*self.forward_path))
        self.axs[2, 1].scatter(s, s_dot_max, s=2)   
    
    def plot_backward_path(self):
        # plot backward path
        s, s_dot_max, _ = list(zip(*self.backward_path))
        self.axs[2, 1].scatter(s, s_dot_max, s=2)   

    def plot_path(self):
        s, s_dot_max, _ = list(zip(*self.final_path))
        self.axs[2, 1].scatter(s, s_dot_max, s=2)
    
    def plot_intersection_points(self):
        # plot intersection pts
        for i in range(1, 20):
            self.fig, self.axs = plt.subplots(2, 2, figsize=(10, 8))
            self.plot_limit_curve()
            self.plot_path()
            self.plot_inflection_pts()
            lower_bound = int((i-1) * len(self.intersection_pts)/20)
            upper_bound = int(i * len(self.intersection_pts)/20)
            s, s_dot_max = list(zip(*self.intersection_pts[lower_bound:upper_bound]))
            self.axs[2, 1].scatter(s, s_dot_max, s=2)
            plt.show()

    def plot_inflection_pts(self):
        # plot inflection pts
        s, s_dot_max = list(zip(*self.inflection_pts))
        self.axs[2, 1].scatter(s, s_dot_max, color='red', s=8)   
 

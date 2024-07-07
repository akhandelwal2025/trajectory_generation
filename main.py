from Path import *
from Limit import *
import constants
import matplotlib.pyplot as plt
import numpy as np
import time

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
        
        self.limit_curve = LimitCurve(self.joint_paths)

        # prev timestep - used for collision-checking
        self.prev_s = 0
        self.prev_s_dot = 0

        # curr timestep
        self.curr_s = 0
        self.curr_s_dot = 0

        # store (s, s_dot) pairs for viz purposes
        self.path = []

        # matplotlib stuff
        self.fig, self.axs = plt.subplots(2, 2, figsize=(10, 8))

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
        while not self.done():
            # step 2
            while not self.find_collisions():
                self.integrate_forward()

            break
    """
        two cases for collisions
            a. abs(s_dot - limit_curve(s)) <= epsilon
                i. do nothing - s_dot is within bounds for a collision
            b. s_dot > limit_curve(s) && abs(s_dot - limit_curve(s)) > epsilon
                i. need to go backward to find a precise collision point
    """
    def find_collisions(self):
        lim_curve_s_dot = self.limit_curve.evaluate(self.curr_s)
        if self.curr_s_dot < lim_curve_s_dot and (lim_curve_s_dot - self.curr_s_dot) > constants.epsilon:
            print(f"     lim_curve: {lim_curve_s_dot} | self.curr_s_dot: {self.curr_s_dot} | diff: {lim_curve_s_dot - self.curr_s_dot}")
            self.path.append((self.curr_s, self.curr_s_dot))
            return False
        
        # case a
        if abs(lim_curve_s_dot - self.curr_s_dot) <= constants.epsilon:
            self.path.append((self.curr_s, self.curr_s_dot))
            return True

        # case b - find (s, s_dot) for precise collision point
        lower_timestep = 0
        upper_timestep = constants.timestep
        while self.curr_s_dot > lim_curve_s_dot and abs(lim_curve_s_dot - self.curr_s_dot) > constants.epsilon: 
            print("stuck in this loops")
            print(f"     lim_curve: {lim_curve_s_dot} | self.curr_s_dot: {self.curr_s_dot} | diff: {lim_curve_s_dot - self.curr_s_dot}")
            print(f"     lower: {lower_timestep} | upper: {upper_timestep}")
            timestep = (lower_timestep + upper_timestep) / 2
            max_s_dot2 = self.calc_max_s_dot2(self.prev_s)
            self.curr_s_dot = self.prev_s_dot + max_s_dot2 * timestep # v_final = v_initial + a*delta_t
            self.curr_s = self.prev_s + (self.prev_s_dot * timestep) + (0.5 * max_s_dot2 * (timestep ** 2))# s_final = s_initial + v_initial*delta_t + 1/2*a*delta_t^2 
    
            # if resulting s_dot is above lim_curve_s_dot, then set upper_timestep = timestep
            lim_curve_s_dot = self.limit_curve.evaluate(self.curr_s)
            if self.curr_s_dot > lim_curve_s_dot:
                upper_timestep = timestep
            else:
                lower_timestep = timestep
            time.sleep(0.05)
        self.path.append((self.curr_s, self.curr_s_dot))
        return True

    def integrate_forward(self):
        # curr_s, curr_s_dot will be overwritten. need to store them in prev_s, prev_s_dot
        self.prev_s = self.curr_s
        self.prev_s_dot = self.curr_s_dot

        max_s_dot2 = self.calc_max_s_dot2(self.prev_s)
        self.curr_s_dot = self.prev_s_dot + max_s_dot2 * constants.timestep # v_final = v_initial + a*delta_t
        self.curr_s = self.prev_s + (self.prev_s_dot * constants.timestep) + (0.5 * max_s_dot2 * (constants.timestep ** 2)) # s_final = s_initial + v_initial*delta_t + 1/2*a*delta_t^2 
    
    # check if trajectory has reached the end 
    def done(self):
        return self.curr_s >= 1

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
    # technically should be a function of s-dot however since path is not second-order differentiable
    # f-prime2(s) is always 0 so that term goes to zero. thus this is just a function of s
    def calc_min_s_dot2(self, s):
        min_s_dot2 = float('-inf')
        for joint_i in self.joint_paths:
            if joint_i.f_prime(s) != 0:
                min_s_dot2 = max(min_s_dot2, (-1*joint_i.q_dot2_max)/abs(joint_i.f_prime(s)))
        return min_s_dot2

    # implementing eqn 23
    # same note as previous function - technically should be function of s-dot
    def calc_max_s_dot2(self, s):
        max_s_dot2 = float('inf')
        for joint_i in self.joint_paths:
            if joint_i.f_prime(s) != 0:
                max_s_dot2 = min(max_s_dot2, joint_i.q_dot2_max/abs(joint_i.f_prime(s)))
        return max_s_dot2

    def plot_segments(self):
        for i, joint_path in enumerate(self.joint_paths):
            s = np.linspace(0, 1, 500)[:-1]
            theta = [joint_path.f(i) for i in s]

            if i == 0:
                self.axs[0, 0].scatter(s, theta, s=2)
            elif i == 1:
                self.axs[0, 1].scatter(s, theta, s=2)
            else:
                self.axs[1, 0].scatter(s, theta, s=2)
    
    def plot_limit_curve(self):
        # plot limit curve
        s = np.linspace(0, 1, 500)[:-1]
        s_dot_max = [self.limit_curve.evaluate(i) for i in s]
        self.axs[1, 1].scatter(s, s_dot_max, s=2)

        # plot path
        s, s_dot_max = list(zip(*self.path))
        self.axs[1, 1].scatter(s, s_dot_max, s=2)    

def main():
    traj = Trajectory(constants.waypoints, constants.velocity_constraints, constants.accel_constraints)
    traj.plot_segments()
    traj.generate_trajectory()
    traj.plot_limit_curve()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()


        
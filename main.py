from Path import *
from Limit import *
import constants
import matplotlib.pyplot as plt
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
        
        self.limit_curve = LimitCurve(self.joint_paths)

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
        pass

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
                min_s_dot2 = max(min_s_dot2, (-1*joint_i.accel_constraint)/abs(joint_i.f_prime(s)))
        return min_s_dot2

    # implementing eqn 23
    # same note as previous function - technically should be function of s-dot
    def calc_max_s_dot2(self, s):
        max_s_dot2 = float('inf')
        for joint_i in self.joint_paths:
            if joint_i.f_prime(s) != 0:
                max_s_dot2 = min(max_s_dot2, joint_i.accel_constraint/abs(joint_i.f_prime(s)))
        return max_s_dot2

    def plot_segments(self):
        for i, joint_path in enumerate(self.joint_paths):
            s = np.linspace(0, 1, 500)[:-1]
            theta = [joint_path.f(i) for i in s]

            if i == 0:
                self.axs[0, 0].scatter(s, theta)
            elif i == 1:
                self.axs[0, 1].scatter(s, theta)
            else:
                self.axs[1, 0].scatter(s, theta)
    
    def plot_limit_curve(self):
        s = np.linspace(0, 1, 500)[:-1]
        s_dot_max = [self.limit_curve.evaluate(i) for i in s]
        self.axs[1, 1].scatter(s, s_dot_max)

def main():
    traj = Trajectory(constants.waypoints, constants.velocity_constraints, constants.accel_constraints)
    traj.plot_segments()
    traj.plot_limit_curve()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()


        
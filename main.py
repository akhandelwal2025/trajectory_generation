from Path import *
from Limit import *
import matplotlib.pyplot as plt
import numpy as np

class Trajectory:
    def __init__(self, waypoints, velocity_constraints):
        """
            waypoints = [[q11, q12, q13, ..., q1m], 
                         [q21, q22, q23, ..., q2m], 
                         ....,
                         [qn1, qn2, qn3, ..., qnm]]
            n = number of waypoints
            m = number of joints
        """
        self.joint_paths = []
        all_joints = list(zip(*waypoints)) # groups all joint positions together
        assert(len(all_joints) == len(velocity_constraints))
        for i in range(len(all_joints)):
            joint_i_positions = all_joints[i]
            joint_i_vel_constraint = velocity_constraints[i]
            self.joint_paths.append(JointPath(joint_i_positions, joint_i_vel_constraint))
        
        self.limit_curve = LimitCurve(self.joint_paths)

        # matplotlib stuff
        self.fig, self.axs = plt.subplots(2, 2, figsize=(10, 8))

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

waypoints = [
    [30, 20, 10],
    [90, 40, 20],
    [60, 60, 50],
    [60, 80, 10]
] # deg
velocity_constraints = [180, 180, 180] # deg/s

traj = Trajectory(waypoints, velocity_constraints)
traj.plot_segments()
traj.plot_limit_curve()
plt.tight_layout()
plt.show()


        
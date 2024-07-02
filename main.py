from Path import *
import matplotlib.pyplot as plt
import numpy as np

class Trajectory:
    def __init__(self, waypoints):
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
        for joint_positions in all_joints:
            self.joint_paths.append(Path(joint_positions))

    def plot_segments(self):
        for joint_path in self.joint_paths:
            x = np.linspace(0, 1, 500)[:-1]
            y = [joint_path.evaluate(i) for i in x]
            plt.scatter(x, y)
            plt.show()

waypoints = [
    [30, 20, 10],
    [90, 40, 20],
    [60, 60, 50],
    [60, 80, 10]
]

traj = Trajectory(waypoints)
traj.plot_segments()

        
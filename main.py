from Trajectory import *
from Path import *
from Limit import *
import constants
import matplotlib.pyplot as plt

def main():
    traj = Trajectory(constants.waypoints, constants.velocity_constraints, constants.accel_constraints)
    traj.plot_segments()
    traj.generate_trajectory()
    traj.plot_limit_curve()
    traj.plot_path()
    traj.plot_inflection_pts()
    # traj.plot_intersection_points()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()


        
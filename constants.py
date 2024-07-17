import math

timestep = 0.01 #sec
epsilon = 0.005 #s/sec

inflection_epsilon = 0.001 # s
parabolic_epsilon = 0.00001 # s
s_theta_sampling_frequency = 5000 # num_pts
discretization = 1000
blend_radius = 10 # deg
# sandbox examples
# waypoints = [
#     [30, 20, 180],
#     [90, 40, 20],
#     [60, 60, 50],
#     [60, 80, 10]
# ] # deg

# waypoints = [
#     [60, 80, 10],
#     [60, 60, 50],
#     [90, 40, 20],
#     [30, 20, 180],
# ] # deg

# actual examples
waypoints = [
    [10, 30, 90],
    [20, 90, 180],
    [30, 60, 90],
    [40, 30, 0],
] # deg

# waypoints = [
#     [0, 0, 90],
#     [0, 0, 180],
#     [0, 0, 90],
#     [0, 0, 0],
# ] # deg

velocity_constraints = [180, 180, 180] # deg/s
accel_constraints = [50, 50, 50] # deg^2/s
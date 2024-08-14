import math
import utils

timestep = 0.01 #sec
epsilon = 0.005 #s/sec

accel_curve_inflection_epsilon = 5.0 # unitless (technically slope)
vel_curve_inflection_epsilon = 0.001 # s
parabolic_epsilon = 0.00001 # s
s_theta_sampling_frequency = 500 # num_pts
accel_curve_sampling_frequency = 1000 # num_pts
discretization = 4
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
# waypoints = [
#     [180.63, -18.26, 42.01, -30.56, 0.21, -0.08],
#     [195.83, -38.21, 41.92, 18.95, -5.26, -0.08],
#     [195.04, -38.12, 86.34, 19.00, 68.34, -0.03]
# ]
# waypoints = [
#     [10, 30, 90],
#     [20, 90, 180],
#     [30, 60, 90],
#     [40, 30, 0],
# ] # deg

# waypoints = [
#     [0, 0, 90],
#     [0, 0, 180],
#     [0, 0, 90],
#     [0, 0, 0],
# ] # deg
# waypoints = [
#     [3.83501197, -0.89203778, 1.06465, -0.33911747, 2.4389231, -0.08831366],
#     [3.33270621, 0.04799655, 0.8990191, -1.3585643, 3.10598794, -0.08883726],
#     [2.10451801, -1.153139, 0.90861841, -0.8290314, 5.79588938, 0.003665191],
#     [3.1222195, -0.38763763, 0.8857546, -3.51596578, 6.2081361, -0.04049164],
#     [3.60445397, -0.23282692, 0.88784899, -0.86515971, 6.1845742, -0.0191986],
#     [3.83501197, -0.89203778, 1.06465, -0.33911747, 2.4389231, -0.08831366],
# ]
# waypoints = utils.waypoints_rad_2_deg(waypoints)
waypoints = [
    [350.0, -37.64, 44.20, -55.36, 112.92, -5.20],
    [310.52, -121.87, 87.25, -36.36, -9.41, -5.36],
    [240.65, -31.84, 87.21, -81.80, -1.27, -5.26],
    [116.6, -29.48, 76.17, -183.45, -1.69, -23.34],
    [186.67, 5.72, -5.74, -90.82, -1.62, -23.33],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
]
velocity_constraints = [180, 180, 180, 180, 180, 180] # deg/s
accel_constraints = [30, 30, 30, 30, 30, 30] # deg^2/s
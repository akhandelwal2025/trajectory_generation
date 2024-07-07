timestep = 0.05 #sec
epsilon = 0.005 #s/sec

# waypoints = [
#     [30, 20, 180],
#     [90, 40, 20],
#     [60, 60, 50],
#     [60, 80, 10]
# ] # deg

waypoints = [
    [60, 80, 10],
    [60, 60, 50],
    [90, 40, 20],
    [30, 20, 180],
] # deg


velocity_constraints = [180, 180, 180] # deg/s
accel_constraints = [100, 100, 100] # deg^2/s
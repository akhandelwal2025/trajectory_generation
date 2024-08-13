import math

def waypoints_rad_2_deg(waypoints):
    for i in range(len(waypoints)):
        for j in range(len(waypoints[i])):
            waypoints[i][j] = math.degrees(waypoints[i][j])
    return waypoints
import math

def l2_norm(p1, p2):
    segment = (p2[0]-p1[0], p2[1]-p1[1]) # p2 - p1
    return math.sqrt(segment[0] ** 2 + segment[1] ** 2) # calculate magnitude of the segment
 
def unit_vector(p1, p2):
    segment = (p2[0]-p1[0], p2[1]-p1[1]) # p2 - p1
    mag = l2_norm(p1, p2)
    return (segment[0]/mag, segment[1]/mag) # normalize the segment by dividing by its magnitude

def dot(unit1, unit2):
    return (unit1[0] * unit2[0]) + (unit1[1] * unit2[1]) 
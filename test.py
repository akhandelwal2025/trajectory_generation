import numpy as np
import matplotlib.pyplot as plt

# Define the path segments
segments = [
    ((0, 0), (1, 2)),
    ((1, 2), (3, 3)),
    ((3, 3), (4, 0)),
    # Add more segments as needed
]

# Length of the blending region
blend_length = 0.75  # Adjust this length as needed

# Compute blending regions
blending_regions = []
for i in range(len(segments) - 1):
    seg1 = segments[i]
    seg2 = segments[i + 1]
    
    # Define the end of the first segment and start of the second segment
    blending_regions.append((seg1[1], seg2[0]))

# Parabolic blend function
def parabolic_blend(p1, p2, blend_length, t):
    """
    Compute the parabolic blend between points p1 and p2 over blend_length at parameter t.
    t ranges from 0 to 1 within the blending region.
    """
    x1, y1 = p1
    x2, y2 = p2
    
    # Parabolic blend equations
    a = (y2 - y1) / (blend_length ** 2)
    b = -2 * a * blend_length
    c = y1
    
    xt = x1 + (x2 - x1) * t
    yt = a * (xt - x1) ** 2 + b * (xt - x1) + c
    
    return xt, yt

# Generate the path function
def generate_path(segments, blending_regions, blend_length, steps_per_segment=100):
    path = []
    
    for i in range(len(segments) - 1):
        seg1 = segments[i]
        seg2 = segments[i + 1]
        blend_region = blending_regions[i]
        
        # Linear segment before blending region
        for j in range(steps_per_segment - int(blend_length * steps_per_segment)):
            t = j / (steps_per_segment - int(blend_length * steps_per_segment))
            x = seg1[0][0] + t * (seg1[1][0] - seg1[0][0])
            y = seg1[0][1] + t * (seg1[1][1] - seg1[0][1])
            path.append((x, y))
        
        # Parabolic blend region
        for j in range(int(blend_length * steps_per_segment)):
            t = j / (blend_length * steps_per_segment)
            x, y = parabolic_blend(blend_region[0], blend_region[1], blend_length, t)
            path.append((x, y))
    
    # Add the last segment
    last_segment = segments[-1]
    for j in range(steps_per_segment):
        t = j / steps_per_segment
        x = last_segment[0][0] + t * (last_segment[1][0] - last_segment[0][0])
        y = last_segment[0][1] + t * (last_segment[1][1] - last_segment[0][1])
        path.append((x, y))
    
    return path

# Visualize the path
path = generate_path(segments, blending_regions, blend_length)

x, y = zip(*path)
plt.plot(x, y)
plt.scatter(*zip(*[seg[0] for seg in segments] + [segments[-1][1]]), color='red')  # Mark segment points
plt.title("Multi-Segment Path with Parabolic Blends")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

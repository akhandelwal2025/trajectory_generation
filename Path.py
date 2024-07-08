import math

class LinearJointSegment:
    def __init__(self, m, offset, b):
        self.m = m
        self.offset = offset
        self.b = b

    def f(self, s):
        return self.m * (s - self.offset) + self.b

    def f_prime(self, s):
        return self.m

    def f_prime2(self, s):
        return 0

class JointPath:
    def __init__(self, joint_positions, vel_constraint, accel_constraint):
        # joint_positions = list of thetas for a single joint
        # ex. joint1 = [30, 90, 60, 60]
        # vel_constraint = q-dot-max for the joint
        self.q_dot_max = vel_constraint
        self.q_dot2_max = accel_constraint
        self.segments = []
        self.s_interval = 1/len(joint_positions)
        start_pos = 0
        for i, end_pos in enumerate(joint_positions):
            m = (end_pos - start_pos)/self.s_interval
            offset = i * self.s_interval
            b = start_pos
            self.segments.append(LinearJointSegment(m, offset, b))
            start_pos = end_pos
        print(len(self.segments))
    
    def f(self, s):
        if s == 1:
            return self.segments[-1].f(s)
            
        segment_idx = math.floor(s/self.s_interval)
        return self.segments[segment_idx].f(s)

    def f_prime(self, s):
        if s == 0:
            return self.segments[0].f_prime(s)
        if s == 1:
            return self.segments[-1].f_prime(s)
        
        # since joint paths are implemented as linear segments, the path is fully continuous, but not differentiable at the s_intervals
        # need to handle these cases separately
        if s % self.s_interval == 0:
            inflection_idx = int(s/self.s_interval)
            seg1_m = self.segments[inflection_idx-1].f_prime(s)
            seg2_m = self.segments[inflection_idx].f_prime(s)
            # if the sign of the slope of consecutive segments changes, then s-dot should be 0
            if math.copysign(1, seg1_m) != math.copysign(1, seg2_m):     
                return 0
            else:
                return (seg1_m + seg2_m)/2 # TODO: average of two slopes. should this be changed? 
        else:
            segment_idx = math.floor(s/self.s_interval)
            return self.segments[segment_idx].f_prime(s)

    def f_prime2(self, s):
        segment_idx = math.floor(s/self.s_interval)
        return self.segments[segment_idx].f_prime2(s)
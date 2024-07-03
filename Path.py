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

    def f_double_prime(self, s):
        return 0

class JointPath:
    def __init__(self, joint_positions, velocity_constraint):
        # joint_positions = list of thetas for a single joint
        # ex. joint1 = [30, 90, 60, 60]
        # velocity_constraint = q-dot-max for the joint
        self.q_dot_max = velocity_constraint
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
        segment_idx = math.floor(s/self.s_interval)
        return self.segments[segment_idx].f(s)

    def f_prime(self, s):
        segment_idx = math.floor(s/self.s_interval)
        return self.segments[segment_idx].f_prime(s)

    def f_double_prime(self, s):
        segment_idx = math.floor(s/self.s_interval)
        return self.segments[segment_idx].f_double_prime(s)
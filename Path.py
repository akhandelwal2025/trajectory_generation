import math

class LinearSegment:
    def __init__(self, m, offset, b):
        self.m = m
        self.offset = offset
        self.b = b

    def evaluate(self, x):
        return self.m * (x - self.offset) + self.b

    def evaluate_first_deriv(self, x):
        return self.m

    def evaluate_second_deriv(self, x):
        return 0

class Path:
    def __init__(self, joint_positions):
        # joint_positions = list of thetas for a single joint
        # ex. joint1 = [30, 90, 60, 60]
        self.segments = []
        self.s_interval = 1/len(joint_positions)
        start_pos = 0
        for i, end_pos in enumerate(joint_positions):
            m = (end_pos - start_pos)/self.s_interval
            offset = i * self.s_interval
            b = start_pos
            self.segments.append(LinearSegment(m, offset, b))
            start_pos = end_pos
        print(len(self.segments))
    def evaluate(self, x):
        segment_idx = math.floor(x/self.s_interval)
        return self.segments[segment_idx].evaluate(x)
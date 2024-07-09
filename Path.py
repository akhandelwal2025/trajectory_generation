import math

class LinearJointSegment:
    def __init__(self, start, end):
        # start = (s1, s_dot1), end = (s2, s_dot2)
        self.start = start
        self.end = end
        self.m, self.offset, self.b = self.calc_m_offset_b()

    def calc_m_offset_b(self):
        m = (self.end[1] - self.start[1]) / (self.end[0] - self.start[0])
        offset = self.start[0]
        b = self.start[1]
        return m, offset, b

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
        self.joint_positions = joint_positions
        self.s_interval = 1/len(self.joint_positions)
        self.q_dot_max = vel_constraint
        self.q_dot2_max = accel_constraint

        self.segments = self.construct_linear_segments()
        # self.segments, self.segment_start_s = self.add_parabolic_blends()
    
    """
        - for each pair of consecutive segments, first identify if a blend is needed by checking if slopes are different
        - if blend needed, then follow algorithm laid out in paper or moveit library: https://github.com/moveit/moveit/blob/master/moveit_core/trajectory_processing/src/time_optimal_trajectory_generation.cpp#L97
    """
    # def construct_parabolic_segments(self):
    #     for i in range(0, len(self.linear_segments)):
    #         seg1, seg2 = self.linear_segments[i], self.linear_segments[i+1]
    #         intersection_s = self.linear_starts[i+1]
    #         if seg1.m != seg2.m:
                
    # def add_parabolic_blends(self, seg1: LinearJointSegment, seg2: LinearJointSegment):
    #     start_dir = 

    def construct_linear_segments(self):
        start = (0, 0) # (s, theta)
        segments = []
        for i, end_theta in enumerate(self.joint_positions):
            end_s = (i+1) * self.s_interval
            end = (end_s, end_theta)
            segments.append(LinearJointSegment(start, end))
            start = end
        return segments
        # start_pos = 0
        # linear_segments = [] # list of LinearJointSegment
        # linear_starts = [] # starting s of each LinearJointSegment
        # for i, end_pos in enumerate(self.joint_positions):
        #     m = (end_pos - start_pos)/self.s_interval
        #     offset = i * self.s_interval
        #     b = start_pos
        #     linear_segments.append(LinearJointSegment(m, offset, b))
        #     linear_starts.append(i * self.s_interval)
        #     start_pos = end_pos
        # return linear_segments, linear_starts
    
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
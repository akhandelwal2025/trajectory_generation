import constants
import math

class AccelerationLimitCurve:
    def __init__(self, joint_paths):
        self.joint_paths = joint_paths
    
    # implementing eqn 31
    def evaluate(self, s):
        s_dot_max = float('inf')
        for i, joint_i in enumerate(self.joint_paths):
            f_i_prime = joint_i.f_prime(s)
            f_i_prime2 = joint_i.f_prime2(s)

            # case 2, f_i_prime = 0 and f_i_prime2 != 0
            if f_i_prime == 0 and f_i_prime2 != 0:
                s_dot_max = min(s_dot_max, math.sqrt(joint_i.q_dot2_max / abs(joint_i.f_prime2(s))))
            
            # case 1
            for joint_j in self.joint_paths[i+1:]:
                f_j_prime = joint_j.f_prime(s)
                f_j_prime2 = joint_j.f_prime2(s)
                if f_i_prime != 0 and f_j_prime != 0:
                    denom = abs((f_i_prime2/f_i_prime) - (f_j_prime2/f_j_prime))
                    if denom != 0:
                        numer = (joint_i.q_dot2_max / abs(f_i_prime)) + (joint_j.q_dot2_max / abs(f_j_prime))
                        s_dot_max = min(s_dot_max, math.sqrt(numer/denom))
        return s_dot_max

class VelocityLimitCurve:
    def __init__(self, joint_paths):
        self.joint_paths = joint_paths

    # inflection pts are s-values where the limit curve changes
    # we want the lower value of each change
    def find_inflection_points(self):
        num_inflection_pts = len(self.joint_paths[0].segments)
        s_interval = 1/num_inflection_pts
        self.inflection_pts = []
        # suppose there are n-inflection pts. we want to find the s-dot value for the first n-1
        # the last inflection pt is always (s=1.0, s_dot=0), because we want the path to end at zero velocity
        for i in range(1, num_inflection_pts):
            s = s_interval * i
            s_dot = min(self.evaluate(s-constants.inflection_epsilon), self.evaluate(s+constants.inflection_epsilon))
            self.inflection_pts.append((s, s_dot))
        self.inflection_pts.append((1.0, 0.0))
        return self.inflection_pts
    
    # this is just the velocity limit imposed by constraints on velocity
    # implementing eqn 36
    def evaluate(self, s):
        s_dot_max = float('inf')
        for joint_i in self.joint_paths:
            # if f_prime(s) is 0, then that means joint isn't moving
            # therefore, it imposes no constraint on the path velocity
            if joint_i.f_prime(s) != 0: 
                s_dot_max = min(s_dot_max, joint_i.q_dot_max/abs(joint_i.f_prime(s)))
        return s_dot_max
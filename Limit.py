class LimitCurve:
    def __init__(self, joint_paths):
        self.joint_paths = joint_paths
    
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
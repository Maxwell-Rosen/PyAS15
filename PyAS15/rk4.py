def push(velocity: callable, Zi, step_size):
    # Define a function rk4_step(f, Ri, Zi, velocity: callable, step_size)
    #  that takes in a function f, arguments Ri and Zi, a velocity function velocity: callable,
    #  and a step size and returns the next value of R,Z using an rk4 algorithm.
    vZ1 = velocity(Zi)
    k1_Z = step_size * vZ1

    vZ2 = velocity(Zi + 0.5 * k1_Z)
    k2_Z = step_size * vZ2

    vZ3 = velocity(Zi + 0.5 * k2_Z)
    k3_Z = step_size * vZ3

    vZ4 = velocity(Zi + k3_Z)
    k4_Z = step_size * vZ4

    Z_new = Zi + (1/6) * (k1_Z + 2*k2_Z + 2*k3_Z + k4_Z)
    return Z_new
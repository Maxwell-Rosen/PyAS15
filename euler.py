def push(velocity: callable, Zi, step_size):
    # Define a function rk4_step(f, Ri, Zi, velocity: callable, step_size)
    #  that takes in a function f, arguments Ri and Zi, a velocity function velocity: callable,
    #  and a step size and returns the next value of R,Z using an rk4 algorithm.
    vZi = velocity(Zi)
    Z_new = Zi + vZi * step_size
    return Z_new
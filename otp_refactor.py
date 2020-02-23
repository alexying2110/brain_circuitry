import numpy as np

# constants
a_1 = 1.57e1
a_2 = 8.877e1
a_3 = 2.1e0

b_1 = 8.00e-1
b_2 = 9.789e1
b_3 = 1.200e0

g = 1.750e-1
k = 7.089e3
c = 7.534e-2
p = 1
I = 7.774e1

# define R, O, and N
R = 1
O = 1
N = 21

# define time
t_step = 1e-1
t_duration = 3.0
t_range = np.arange(0, t_duration, t_step)

b = 1
d = 132

def peri_receptor_out(alpha, beta, t_range, R, O, N, stimulus = "parabola", form = "intensity"):
    # init x_1, x_2, and u (odor function)
    t_steps = len(t_range)
    t_step = t_range[1] - t_range[0]
    
    x_1 = np.zeros((t_steps, R, O, N))
    x_2 = np.zeros((t_steps, R, O, N))
    u = np.zeros((t_steps, R, O, N))

    # different input parameters

    if stimulus == 'constant':
        for i in range(N):
            u[:,:,:, i] = 1. + 5. * i

    if stimulus == 'step':
        for i in range(N):
            u[(t_range >= 0.5) * (t_range <= 2.5),:,:, i] += 1. + 5. * i
    
    if stimulus == 'ramp':
        for i in range(N):
            ramp = (1. + 5. * i) / 1.8 * ((t_range - 0.5) * (t_range >= 0.5) * (t_range < 2.3) + (1 - 5 * (t_range - 2.3)) * (t_range >= 2.3) * (t_range <= 2.5))
            ramp = np.reshape(ramp, (len(ramp), 1, 1))
            ramp = np.repeat(ramp, R, axis = 1)
            ramp = np.repeat(ramp, O, axis = 2)
            u[:,:,:, i] = ramp

    if stimulus == 'parabola':
        for i in range(N):
            parabola = (1. + 5. * i) / 1.9 * ((t_range - 0.5) ** 2 * (t_range >= 0.5) * (t_range < 2.4) + (1 - 10 * (t_range - 2.4) ** 2) * (t_range >= 2.4) * (t_range <= 2.5))
            parabola = np.reshape(parabola, (len(parabola), 1, 1))
            parabola = np.repeat(parabola, R, axis = 1)
            parabola = np.repeat(parabola, O, axis = 2)
            u[:,:,:, i] = parabola

    if form == 'gradient':
        u[:-1] = (u[1:] - u[:-1]) / t_step
        u[-1] = 0

    for i in range(t_steps - 1):
        dx_1 = x_2[i,:,:,:]
        dx_1[np.isnan(dx_1)] = 0.

        dx_2 = alpha ** 2 * (u[i,:,:,:] - x_1[i,:,:,:]) - 2 * alpha * beta * x_2[i,:,:,:]
        dx_2[np.isnan(dx_2)] = 0.

        x_1[i + 1,:,:,:] = x_1[i,:,:,:] + t_step * dx_1
        x_2[i + 1,:,:,:] = x_2[i,:,:,:] + t_step * dx_2

    print(u)
    print(x_1)
    return x_1

def otp_model(alpha_1, beta_1, alpha_2, beta_2, alpha_3, beta_3, kappa, g, c, p, I, t_range, R, O, N, b, d, filter_h = peri_receptor_out):
    u_h = filter_h(alpha_1, beta_1, t_range, R, O, N, stimulus = "parabola")
    u_h_gradient = filter_h(alpha_1, beta_1, t_range, R, O, N, stimulus = "parabola", form = "gradient")
    v = np.maximum(0., u_h + u_h_gradient * g)
    # a weird fucking hack to get vectors filled with b what the fuck?
    b = v * 0.0 + b
    d = v * 0.0 + d

    t_steps = len(t_range)
    t_step = t_range[1] - t_range[0]
    x_1 = v * 0.0
    x_2 = v * 0.0
    x_3 = v * 0.0
    I = v * 0.0

    for i in range(t_steps - 1):
        dx_1 = b[i,:,:,:] * v[i,:,:,:] * (1 - x_1[i,:,:,:]) - d[i,:,:,:] * x_1[i,:,:,:]
        dx_2 = alpha_2 * x_1[i,:,:,:] * (1 - x_2[i,:,:,:]) * beta_2 * x_2[i,:,:,:] - kappa * np.power(np.maximum(x_2[i,:,:,:], 0), 2 / 3) * np.power(np.maximum(x_3[i,:,:,:], 0), 2 /3)
        dx_3 = alpha_3 * x_2[i,:,:,:] - beta_3 * x_3[i,:,:,:]

        dx_1[np.isnan(dx_1)] = 0.
        dx_2[np.isnan(dx_2)] = 0.
        dx_3[np.isnan(dx_3)] = 0.

        x_1[i + 1,:,:,:] = x_1[i,:,:,:] + t_step * dx_1
        x_2[i + 1,:,:,:] = x_2[i,:,:,:] + t_step * dx_2
        x_3[i + 1,:,:,:] = x_3[i,:,:,:] + t_step * dx_3



peri_receptor_out(a_1, b_1, t_range, R, O, N, stimulus = "parabola")
otp_model(a_1, b_1, a_2, b_2, a_3, b_3, k, g, c, p, I, t_range, R, O, N, b, d)

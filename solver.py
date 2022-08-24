import numpy as np
from numba import jit


# =============================================================================
# Diffusion coefficient (using Prandtl analogy)
# =============================================================================
@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def D_C(tau_b, H,S) :
    #u_s = (1000*9.81*H*S)**(0.5) # shear velocity
    u_s = (tau_b / 1000)**0.5 # shear velocity
    D = 1/15 * u_s * H # diffusion coefficient
    return D

# =============================================================================
# Diffusion solving forward-difference in time, central-difference in space
# =============================================================================
@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def diffusion(b,D,dx,dy,dt):
    b[2:-2, 2:-2] = b[2:-2, 2:-2] + D[2:-2, 2:-2] * dt * \
          ((b[3:-1, 2:-2] - 2*b[2:-2, 2:-2] + b[1:-3, 2:-2])/dx**2 \
          + (b[2:-2, 3:-1] - 2*b[2:-2, 2:-2] + b[2:-2, 1:-3])/dy**2 )    
    return b

# =============================================================================
# Lax-Wendroff method (vectorization)
# =============================================================================
@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def LW(s1d, u1d, dx, nx, dt):
    s2 = np.zeros(nx)
    courant = u1d*dt/dx
    sigma = courant*courant / 2.0
    s2 = s1d[2:-2]-courant*0.5*(s1d[3:-1]-s1d[1:-3]) \
    + sigma*(s1d[3:-1]-2*s1d[2:-2]+s1d[1:-3])

    return s2

# =============================================================================
# Advection solving with directional splitting 
# =============================================================================
@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def advection_x(a, b, c, u, nx, dx, dt):
    for y in range(2, nx+2) :
        s1d_x = a[b,y,:]
        x1d = u[c,y-2,:]
        a[b+1,y,2:-2] = LW(s1d_x,x1d,dx,nx,dt)   
    return a[b+1,2:-2,2:-2]

@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def advection_y(a, b, c, u, ny, dy, dt):
    for k in range(2, ny+2):
        s1d_y = a[b+1,:,k]
        y1d = u[c,:,k-2]
        a[b+1,2:-2,k] = LW(s1d_y,y1d,dy,ny,dt)   
    return a[b+1,2:-2,2:-2]
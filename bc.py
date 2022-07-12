# Adding two-boundaries on each side of the model domain as berms (zero-gradient)
# At the flume side (left), let water flow downstream
def two_boundary(a,i):
    a[i,0,:] = a[i,2,:] # top
    a[i,1,:] = a[i,2,:] # top
    a[i,:,0] = 0  # left
    a[i,:,1] = 0  # left
        
    a[i,-1,:] = a[i,-3,:] # bottom
    a[i,-2,:] = a[i,-3,:] # bottom
    a[i,:,-1] = a[i,:,-3] # right
    a[i,:,-2] = a[i,:,-3] # right    
    return a

# Adding two-boundaries on each side of the model domain (zero-gradient)
# This can be no-flux boundary condition when solving advection and diffisuion seperately.
def two_boundary_origin(a,i):
    a[i,0,:] = a[i,2,:] # top
    a[i,1,:] = a[i,2,:] # top
    a[i,:,0] = a[i,:,2] # left
    a[i,:,1] = a[i,:,2] # left
        
    a[i,-1,:] = a[i,-3,:] # bottom
    a[i,-2,:] = a[i,-3,:] # bottom
    a[i,:,-1] = a[i,:,-3] # right
    a[i,:,-2] = a[i,:,-3] # right      
    return a
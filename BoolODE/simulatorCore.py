import numpy as np

def noise(x,t):
    # Controls noise proportional to
    # square root of activity
    c = 10.#4.
    return (c*np.sqrt(abs(x)))

def deltaW(N, m, h,seed=0):
    """Generate sequence of Wiener increments for m independent Wiener
    processes W_j(t) j=0..m-1 for each of N time intervals of length h.    
    From the sdeint implementation

    :returns:
        - dW : The [n, j] element has the value W_j((n+1)*h) - W_j(n*h) ( has shape (N, m) )
    """
    np.random.seed(seed)
    return np.random.normal(0.0, h, (N, m))

def eulersde(f,G,y0,tspan,pars,seed=0.,dW=None):
    """
    Adapted from sdeint implementation https://github.com/mattja/sdeint/

    :param f: function defining ODE model. Should take vector of current state, current time, and list of parameter values as arguments.
    :type f: function
    :param pars: List of parameter values
    :type pars: list
    :param y0: list of initial values
    :type y0: list
    :param tspan: Array of timepoints to simulate
    :type tspan: ndarray
    :param seed: Seed to initialize random number generator
    :type seed: float
    :returns:
        - y: Array containing the time course of state variables 
    """
    # From sdeint implementation
    N = len(tspan)
    h = (tspan[N-1] - tspan[0])/(N - 1)
    maxtime = tspan[-1]
    # allocate space for result
    d = len(y0)
    y = np.zeros((N+1, d), dtype=type(y0[0]))

    if dW is None:
        # pre-generate Wiener increments (for d independent Wiener processes):
        dW = deltaW(N, d, h, seed=seed)
    y[0] = y0
    currtime = 0
    n = 0
   
    while currtime < maxtime:
        tn = currtime
        yn = y[n]
        dWn = dW[n,:]
        y[n+1] = yn + f(yn, tn,pars)*h + np.multiply(G(yn, tn),dWn)
        # Ensure positive terms
        for i in range(len(y[n+1])):
            if y[n+1][i] < 0:
                y[n+1][i] = yn[i]
        currtime += h
        n += 1 
    return y

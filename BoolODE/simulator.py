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

def simulateModel(Model, y0, parameters,isStochastic, tspan,seed):
    """Call numerical integration functions, either odeint() from Scipy,
    or simulator.eulersde() defined in simulator.py. By default, stochastic simulations are
    carried out using simulator.eulersde.

    :param Model: Function defining ODE model
    :type Model: function
    :param y0: array of initial values for each state variable
    :type y0: array
    :param parameters: list of parameter values to be used in simulations
    :type parameters: list
    :param isStochastic: User defined parameter. Default = True, can be turned off to perform ODE simulations.
    :type isStochastic: bool
    :param tspan: Time points to simulate
    :type tspan: ndarray
    :param seed: Seed to initialize random number generator
    :type seed: float
    :returns: 
        - P: Time course from numerical integration
    :rtype: ndarray

    """
    if not isStochastic:
        P = odeint(Model,y0,tspan,args=(parameters,))
    else:
        P = eulersde(Model,noise,y0,tspan,parameters,seed=seed)
    return(P)

def getInitialCondition(ss, ModelSpec, rnaIndex,
                        proteinIndex,
                        genelist, proteinlist,
                        varmapper,revvarmapper):
    """
    Calculate the initial values of all state variables. 
    Takes into consideration user defined initial conditions, and computes the steady 
    states of the protein variables based on the estimated values of their corresponding genes.

    :param ss: Steady state array
    :type ss: ndarray
    :param ModelSpec: Dictionary of dictionary specifying the ODE model, containing parameters, initial conditions and equations.
    :type ModelSpec: dict
    :param rnaIndex: list of indices of genes
    :type rnaIndex: list
    :param proteinIndex: List of indices of proteins
    :type proteinIndex: list
    :param genelist: List of names of all genes in the model
    :type genelist: list
    :param proteinlist: List of names of all proteins in the model
    :type proteinlist: list
    :param varmapper: Mapper: {variable name : index}
    :type varmapper: dict
    :param revvarmapper: Mapper: {index : variable name}
    :type revvarmapper: dict
    :returns:
        - newics: List containing new initial conditions
    """

    # Initialize
    new_ics = [0 for _ in range(len(varmapper.keys()))]
    # Set the mRNA ics
    for ind in rnaIndex:
        if ss[ind] < 0:
            ss[ind] = 0.0
        new_ics[ind] =  ss[ind]
        if new_ics[ind] < 0:
            new_ics[ind] = 0
    for p in proteinlist:
        ind = revvarmapper['p_'+p]
        if ss[ind] < 0:
            ss[ind] = 0.0
        new_ics[ind] =  ss[ind]
        if new_ics[ind] < 0:
            new_ics[ind] = 0            
    # Calculate the Protein ics based on mRNA levels
    for genename in genelist:
        pss = ((ModelSpec['pars']['r_' + genename])/\
                                      (ModelSpec['pars']['l_p_' + genename]))\
                                      *new_ics[revvarmapper['x_' + genename]]
        new_ics[revvarmapper['p_' + genename.replace('_','')]] = pss
    return(new_ics)

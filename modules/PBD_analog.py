# Implements the PBD analogy and simulations 

import numpy as np
import scipy 

def simul_prot_etienne(l1, l2, l3, m1, m2, n_g0, n_i0, n_sim = 10000, step = 100, horizon = None):
    issues = ['GI', 'IG', 'II', 'GE', 'IE']
    T, Speciation, Extinction = np.zeros(n_sim), np.zeros(n_sim, dtype=bool), np.zeros(n_sim, dtype = bool)
    T_mean, T_mean_spec,F_speciation,F_extinction = np.zeros(n_sim//step), np.zeros(n_sim//step), np.zeros(n_sim//step),np.zeros(n_sim//step)
    Nb_direct_incipient = np.zeros(n_sim)
    if horizon is None:
        tmax = np.inf
    else:
        tmax = horizon
    for n in range(n_sim):
        n_g, n_i = (n_g0, n_i0)
        m = n_i
        t = 0
        while (n_i + n_g > 0) and not(Speciation[n]) and (t < tmax):
            rates = np.array([l1*n_g, l2 * n_i, l3 * n_i, m1 * n_g, m2*n_i])
            sum_rate = np.sum(rates)
            wait_time = np.random.exponential(1/sum_rate)
            if t + wait_time < tmax:
                event = np.random.choice(issues, p=rates / sum_rate)
                if event == 'GI' or event == 'II':
                    n_i += 1
                    m += int(event == 'GI')
                elif event == 'IG':
                    n_i -= 1
                    n_g += 1
                    Speciation[n] = True
                elif event == 'GE':
                    n_g -= 1
                elif event == 'IE':
                    n_i -= 1
                t += wait_time
            else:
                t = tmax
        Extinction[n] = (n_i + n_g == 0)
        T[n] = t
        Nb_direct_incipient[n] = m
        
        if n % step == step - 1:
            T_mean[n//step] = np.mean(T[0:n+1])
            T_mean_spec[n//step] = np.mean(T[Speciation])
            F_speciation[n//step] = np.sum(Speciation) / (n+1)
            F_extinction[n//step] = np.sum(Extinction) / (n+1)
        
    return dict(T = T, Speciation = Speciation, T_mean = T_mean, T_mean_spec = T_mean_spec, 
                F_speciation = F_speciation, F_extinction = F_extinction, 
                Nb_direct_incipient = Nb_direct_incipient, Extinction = Extinction)

# ---- Constant BD rates and infinite-time probabilities ---- # 
def pi(l1, l2, l3, m1, m2):
    pi_ = (l3 + m2 + l2)/(2*l3) * (1- np.sqrt(1 - 4 * m2 * l3 / ((l3 + m2 + l2)**2))) + np.zeros_like(l1)*np.zeros_like(m1) 
    if type(pi_) == np.ndarray:
        pi_[np.where(l2 == 0)[0]] = 1.0
    elif l2 == 0:
        return 1.0
    return pi_
# the zeros added are for having the same shape as the parameters

def p_speciation(l1, l2, l3, m1, m2):
    pi_ = pi(l1, l2, l3, m1, m2)
    return (1-pi_) * l1 / (m1 + (1-pi_) * l1)

def D(l1, l2, l3, m1, m2):
    return np.sqrt((l2+l3)**2 + 2 * (l2-l3)*m2 + m2**2) + np.zeros_like(l1)*np.zeros_like(m1) 
    # We add zeros for the unused parameters to have an output with the same shape 

def phi(l1, l2, l3, m1, m2):
    return l2 - l3 + m2

def tau(l1, l2, l3, m1, m2):
    D_ = D(l1, l2, l3, m1, m2)
    return 2 / (D_ - l2 + l3 - m2) * np.log(2 / (1+(l2-l3+m2)/D_))

def tau_pdf(t, l1, l2, l3, m1, m2):
    D_ = D(l1, l2, l3, m1, m2)
    phi_ = phi(l1, l2, l3, m1, m2)
    return 2 * D_**2 * np.exp(-D_ * t) * (D_ + phi_) / ((D_ + phi_ + np.exp(-D_ * t) * (D_-phi_))**2) 

def approx_expected_T(l1, l2, l3, m1, m2, branch_at_initiation = False):
    if branch_at_initiation:
        return 1/((1-pi(l1, l2, l3, m1, m2)) * l1 + m1)
    return 1/((1-pi(l1, l2, l3, m1, m2)) * l1 + m1) + tau(l1, l2, l3, m1, m2)

def analog_BD_rates(l1, l2, l3, m1, m2, branch_at_initiation = True):
    p_speciation_ = p_speciation(l1, l2, l3, m1, m2)
    approx_expected_T_ = approx_expected_T(l1, l2, l3, m1, m2, branch_at_initiation)
    spe_rate_ = p_speciation_ / approx_expected_T_
    ext_rate_ = (1-p_speciation_) / approx_expected_T_
    return spe_rate_, ext_rate_

# partial derivative of the function analog_BD_rates
def partial_birth_l1(l1, l2, l3, m1, m2):
    return 1 - pi(l1,l2,l3,m1,m2)

def partial_pi_l2(l1,l2,l3,m1,m2):
    Lambda = l2+l3+m2
    return pi(l1,l2,l3,m1,m2)/Lambda - 2*m2 / (Lambda*np.sqrt(Lambda**2 - 4*l3*m2))

def partial_birth_l2(b,l2,l3,m1,m2):
    return -b*partial_pi_l2(b,l2,l3,m1,m2)

def partial_pi_l3(b,l2,l3,m1,m2):
    Lambda = l2 + l3 + m2
    s = np.sqrt(1 - 4 * l3*m2 / (Lambda**2))
    return -0.5*(l2+m2)/l3**2 * (1-s) + m2 * (l2 + m2 - l3) / (Lambda**2 * l3 * s)

def partial_birth_l3(b, l2,l3,m1,m2):
    return -b * partial_pi_l3(b,l2,l3,m1,m2)

def partial_pi_m2(b,l2,l3,m1,m2):
    Lambda = l2 + l3 + m2
    s = np.sqrt(1 - 4 * l3*m2 / (Lambda**2))
    return (1-s)/(2*l3) - (m2-l3-l2)/(Lambda**2 * s)

def partial_birth_m2(l1,l2,l3,m1,m2):
    return -l1*partial_pi_m2(l1,l2,l3,m1,m2)

def partial_birth_m1(l1, l2, l3, m1, m2):
    return 0.0

def jacobian_analog_bd(i_bd, i_pbd, l1, l2, l3, m1, m2):
    """ Calculates the partial derivatives of the equivalent constant time
    birth rate with respect to the PBD parameter. 

    Args:
        i_bd (int, 0 or 1): 0 for birth rate, 1 for death rate
        i_pbd (int between 0 and 4): variable with respect to which to calculate 
           the derivative (0 for l1, 1 for l2, 2 for l3, 3 for m1, 4 for m2)
        l1, l2, l3, m1, m2 (float): parameters of the PBD model

    Returns:
        float or array: value of the partial derivative of the equivalent birth 
           (i_bd = 0) or death (i_bd = 1) with respect to l1, l2, l3, m1 or m2 
           (depending on i_pbd) evaluated in the parameters l1,l2,l3,m1 and m2
           passed in arguments.
    """
    if i_bd == 0:
        partial_der = [partial_birth_l1, partial_birth_l2, partial_birth_l3, 
                        partial_birth_m1, partial_birth_m2]
        return partial_der[i_pbd](l1, l2, l3, m1, m2)
    elif i_bd == 1:
        if i_pbd == 3:
            return 1.0
        elif i_pbd in (0,1,2,4):
            return 0.0
    raise ValueError("i_bd must be 0 or 1 and i_pbd must be in [0,1,2,3,4].")

def simp_jacobian_analog_bd(i_bd, i_simp_pbd, b, l2, e):
    """Calculates the partial derivatives of the equivalent constant time
    birth rate with respect to the PBD parameter in the simplified framework 
    where l1 = l3 =: b and m1 = m2 =: e. 

    Args:
        i_bd (int, 0 or 1): 0 for birth rate, 1 for death rate
        i_simp_pbd (int 0 or 1 or 2): variable with respect to which to calculate 
           the derivative (0 for b, 1 for l2, 3 for e)
        b, l2, e (float): parameters of the simplified PBD model

    Returns:
        float or array: value of the partial derivative of the equivalent birth 
           (i_bd = 0) or death (i_bd = 1) with respect to b or l2 or e 
           (depending on i_pbd) evaluated in the parameters b, l2 and e passed 
           in arguments.
    """
    corres_simp = [[0,2], # i_simp_pbd = 0 --> l1 and l3
                   [1],   # i_simp_pbd = 1 --> l2
                   [3,4]]  # i_simp_pbd = 2 --> m1 and m2
    if i_simp_pbd in [0,1,2]:
        return np.sum([jacobian_analog_bd(i_bd, i_pbd, b, l2, b, e, e) for i_pbd in corres_simp[i_simp_pbd]])

    raise ValueError("i_bd must be 0 or 1 and i_simp_pbd must be in [0,1,2].")

# ---- Time-dependant BD rates ---- #

# Numerically calculate the probabilities of extinction/completion/speciation
# pIE : probability of extinction of an incipient species within a time t
# pIG : probability of extinction of a good species within a time t
# pIC : probability of completion of an incipient species within a time t
# pGS : probability of speciation from a good species within a time t

# Define the ODEs to be solved (after the log transformation)
def ODE_pGE_log(x, y, l1, l2, l3, m1, m2):
    Theta = l1 + m1
    return pIE_sol(np.log(x) / Theta, l1, l2, l3, m1, m2) * (m1  / Theta * (1-1/x) + l1*y/(Theta*x))

def ODE_pIC_log(x, y, l1, l2, l3, m1, m2):
    Lambda = l2 + l3 + m2
    return (m2 / Lambda * (1-1/x) + 1/x * (1+l3/Lambda * y))**2

def ODE_pGS_log(x, y, T, pIC, l1, l2, l3, m1, m2):
    Theta = l1 + m1
    t = np.log(x) / Theta
    pIC_t = np.interp(t, T, pIC)
    return pIC_t + (1-pIC_t) * l1 * y / (Theta * x)

def pIE_sol(T, l1, l2, l3, m1, m2):
    """ Calculate the exact probability of extinction of an incipient species 
        before a time T. 

    Args:
        T (float or array-like): time horizon (must be non negative)
        l1, l2, l3, m2, m2 (float): parameters of the PBD model

    Returns:
        float or array-like: pIE(T)
    """
    if type(T) == list:
        T = np.array(T)
    
    l = l2 + l3 + m2
    a = m2
    b = l3
    k = np.sqrt(l**2 - 4*a*b)
    c1 = 2.0/(k-l) - 1/k
    etk = np.exp(T * k)
    num = c1 * etk * (l-k) + l/k
    den = np.power(c1 * etk + 1/k, 2)
    pIE =  1/b * np.sqrt(num/den - 0.5*(l*(k-l) + 2*a*b))
    if type(T) == np.ndarray: # if T = 0, proba is zero
        pIE[T == 0.0] = 0.0
    elif T == 0.0:
        pIE = 0.0
    return pIE

def PBD_to_probs(T, l1, l2, l3, m1, m2, solver_method = "BDF", 
                 solver_kwargs = dict(atol = 1e-6, rtol = 1e-9)):
    """ Calculate the probabilities of extinction/speciation/completion before a 
        time T in the PBD model.

    Args:
        T (array): time horizons, must be sorted and T[0] must be 0.
        l1, l2, l3, m2, m2 (float): parameters of the PBD model
        solver_method (str or Solver, optional): method used by solver. Must be 
           accepted by scipy.integrate.solve_ivp. Defaults to "BDF".
        solver_kwargs (dict, optional): other kwargs passed to the ODE solver, 
           in particular relative and absolute tolerances. Consider reducing these 
           tolerances to avoid numerical instabilities. 
           Defaults to dict(atol = 1e-6, rtol = 1e-9).

    Returns:
        dict with following fields:
            T (array):  recall of T
            pIE (array): probability of extinction of an incipient species within a time t
            pIG (array): probability of extinction of a good species within a time t
            pIC (array): probability of completion of an incipient species within a time t
            pGS (array): probability of speciation from a good species within a time t.     
    """
    if T[0] != 0.0:
        raise ValueError("T[0] must be 0.")
    
    par = (l1, l2, l3, m1, m2)
    Theta = l1 + m1
    Lambda = l2 + l3 + m2
    x_l = np.exp(Lambda * T)
    x_t = np.exp(Theta * T)
    
    # Step 1 : solving for pGE
    sol = scipy.integrate.solve_ivp(ODE_pGE_log, t_span = (x_t[0], x_t[-1]), y0 = [0.0], 
                                    method = solver_method, args = par, **solver_kwargs)
    x, g = sol["t"], sol["y"][0,:]
    dgdx = np.diff(g) / np.diff(x)
    xGE = 0.5*(x[1:] + x[:-1])
    tGE = np.log(xGE) / Theta
    pGE_num =  dgdx / pIE_sol(tGE, *par)
    pGE = np.interp(T, tGE, pGE_num, left = 0.0)

    # Step 2 : solving for pIC
    sol = scipy.integrate.solve_ivp(ODE_pIC_log, t_span = (x_l[0], x_l[-1]), y0 = [0.0], 
                                    method = solver_method, args = par, **solver_kwargs)
    x, h = sol["t"], sol["y"][0,:]
    dhdx = np.diff(h) / np.diff(x)
    xIC = 0.5*(x[1:] + x[:-1])
    tIC = np.log(xIC) / Lambda 
    pIC_num = 1 - np.sqrt(dhdx)
    pIC = np.interp(T, tIC, pIC_num, left = 0.0)

    # Step 3 : solving for pGS
    sol = scipy.integrate.solve_ivp(ODE_pGS_log, t_span = (x_t[0], x_t[-1]), y0 = [0.0], 
                                    method = solver_method, args = (T, pIC, *par), 
                                    **solver_kwargs)
    x, m = sol["t"], sol["y"][0,:]
    dmdx = np.diff(m) / np.diff(x)
    xGS = 0.5*(x[1:] + x[:-1])
    tGS = np.log(xGS) / Theta
    pIC_interp = np.interp(tGS, tIC, pIC_num)
    pGS_num = (dmdx - pIC_interp)/(1-pIC_interp)
    pGS = np.interp(T, tGS, pGS_num, left = 0.0)

    return dict(T = T, pIE = pIE_sol(T, *par), pIC = pIC,  pGE = pGE, pGS = pGS)



# Fitting BD model to a PBD model

def time_dep_bd_horiz(pS, pE, t):
    """ Calculate BD rates from probability of speciation/extinction and time.

    Args:
        pS (float or array): probability of speciation before a time t.
        pE (float or array): probability of extinction before a time t.
        t (float or array): time horizon. pS, pE and t must have equal size.

    Returns:
        l: time-dependant birth rate(s), same size as t.
        m: time-dependant death rates(s), same size as t. 
    
    Note:
        if pS+pE == 0.0, nan is returned. 
    """
    if np.size(t) == 1 and pS + pE == 0.0:
        return np.nan 
    a = -np.log(1 - pS - pE) / t
    b = pS / pE
    if type(t) == np.ndarray:
        a[pS + pE == 0.0] = np.nan
        b[pS + pE == 0.0] = np.nan
    return a*b/(1.0+b), a/(1.0+b)

def PBD_to_time_dep_BD(T, l1, l2, l3, m1, m2, solver_method = "BDF", 
                       solver_kwargs = dict(atol = 1e-6, rtol = 1e-9)):
    """ Calculates time dependant BD rates from PBD parameters.

    Args:
        T (array): time values, must be sorted.
        l1 (float > 0): initiation rate from a good species.
        l2 (float > 0): completion rate of an incipient species.
        l3 (float >= 0): initiation rate from an incipient species.
        m1 (float >= 0): extinction rate of a good species. 
        m2 (float >= 0): extinction rate of an incipient species.
        solver_method (str or Solver, optional): method used by solver. Must be 
           accepted by scipy.integrate.solve_ivp. Defaults to "BDF".
        solver_kwargs (dict, optional): other kwargs passed to the ODE solver, 
           in particular relative and absolute tolerances. Consider reducing these 
           tolerances to avoid numerical instabilities. 
           Defaults to dict(atol = 1e-6, rtol = 1e-9).

    Returns:
        l: time-dependant birth rate(s), same size as t.
        m: time-dependant death rates(s), same size as t. 
    """
    probs = PBD_to_probs(T, l1, l2, l3, m1, m2, solver_method=solver_method,
                         solver_kwargs=solver_kwargs)
    l_equiv, m_equiv = time_dep_bd_horiz(probs["pGS"], probs["pGE"], T)
    return np.flip(l_equiv), np.flip(m_equiv)

def find_convergence_point(x, tol = 0.0001, w = 5, silent = False):
    """Find the minimum index i such as local convergence is satisfied. The 
    convergence is satisfied when 
    max(x[i : i+w]) - min(x[i : i+w]) < tol
    The test is made on all dimensions of the first axis of x. The convergence must
    be satisfied for all axes of x. 

    Args:
        x (np.ndarray, 2 dimensional): array to test for convergence. The 
           convergence is tested on the second axis of the object. 
        tol (float > 0, optional): tolerance for the convergence test. 
           Defaults to 0.0001.
        w (int > 1, optional): window size for the convergence test. 
            Defaults to 5.
        silent (bool, optional): if False, prompts a message if the convergence 
            criterion is not satisfied. Defaults to False.




    Returns:
        int or None: minimal index such that the local convergence criterion 
            is met, or None if this criterion is not met. 
    """
    i = 0
    n = np.shape(x)[1]
    dim = np.shape(x)[0]
    cv = False
    while i < n - w and not(cv):
        try_cv = True 
        d = 0
        while d < dim and try_cv:
            window = x[d, i:(i+w)]
            try_cv = max(window) - min(window) < tol # test convergence 
            d+=1
        cv = try_cv
        i += 1 
    if cv:
        return i
    else:
        if not(silent):
            print("Convergence test failed")
        return None



def constant_after_conv(rates, tol = 0.0001, w = 5, silent = False):
    """ Sets an array constant to a local minimum after a certain convergence 
    criteria is met.

    Args:
        rates (np.ndarray, 2 dimensional): array to test for convergence. The 
           convergence is tested on the second axis of the object. 
        tol (float > 0, optional): tolerance for the convergence test. 
           Defaults to 0.0001.
        w (int > 1, optional): window size for the convergence test. 
            Defaults to 5.
        silent (bool, optional): if False, prompts a message if the convergence 
            criterion is not satisfied. Defaults to False.

    Returns:
        np.ndarray, same size as rate: rates with values after the local convergence
            criterion is satisfied set to a local mean, or a copy of rates if 
            no convergence is detected. 
    """
    i_cv = find_convergence_point(rates, tol, w, silent)
    dim = np.shape(rates)[0]
    if not(i_cv is None):
        mean = np.mean(rates[:, i_cv:(i_cv+w)], axis=1)
        rates_cv = np.copy(rates)
        for d in range(dim):
            rates_cv[d,i_cv+w-1:] = mean[d]
        return rates_cv
    return rates 

def forward_PBD_to_time_dep_BD(T, l1, l2, l3, m1, m2, solver_method = "BDF", 
                       solver_kwargs = dict(atol = 1e-6, rtol = 1e-9),
                       smooth_convergence = True,
                       conv_crit = dict(tol = 1e-4, w = 5, silent = False)):
    """ Calculate equivalent BD rates from PBD model based on probabilities 
    of speciation and extinction using a forward approach

    Args:
        T (array): time values, must be sorted.
        l1 (float > 0): initiation rate from a good species.
        l2 (float > 0): completion rate of an incipient species.
        l3 (float >= 0): initiation rate from an incipient species.
        m1 (float >= 0): extinction rate of a good species. 
        m2 (float >= 0): extinction rate of an incipient species.
        solver_method (str or Solver, optional): method used by solver. Must be 
           accepted by scipy.integrate.solve_ivp. Defaults to "BDF".
        solver_kwargs (dict, optional): other kwargs passed to the ODE solver, 
           in particular relative and absolute tolerances. Consider reducing these 
           tolerances to avoid numerical instabilities. 
           Defaults to dict(atol = 1e-6, rtol = 1e-9).
        smooth_convergence (bool, optional): 
        conv_crit (dict, optional): kwargs passed to the convergence function that 
           sets the value of the rates to constant after a certain convergence 
           criterion is met. Useless if smooth_convergence is False. For details 
           see the functions constant_after_conv and find_convergence_point.
           Defaults to dict(tol = 1e-4, w = 5, silent = False).

    Returns:
        l: time-dependant birth rate(s), same size as t.
        m: time-dependant death rates(s), same size as t. 
    """
    probs = PBD_to_probs(T, l1, l2, l3, m1, m2, solver_method=solver_method,
                         solver_kwargs=solver_kwargs)
    pGS, pGE = probs["pGS"], probs["pGE"]
    dpS = np.diff(pGS) / np.diff(T)
    dpE = np.diff(pGE) / np.diff(T)
    rates = np.zeros((2, len(T))) # rates[0,:] = birth, rates[1,:] = death
    rates[:, :-1] = np.array([dpS, dpE]) / (1 - pGS[1:] - pGE[1:])
    
    # the last two points are unstable because of derivation twice
    rates[0, -2:] = rates[0, -3]
    rates[1, -2:] = rates[1, -3]

    # smoothing if convergence detected 
    if smooth_convergence:
        rates = constant_after_conv(rates, **conv_crit)
    return np.flip(rates[0,:]), np.flip(rates[1,:]) # flip to have the present = last value

def fixedZero_PBD_to_varBD(T, l1, l2, l3, m1, m2, solver_method = "BDF", 
                       solver_kwargs = dict(atol = 1e-6, rtol = 1e-9),
                       smooth_convergence = True,
                       conv_crit = dict(tol = 1e-4, w = 5, silent = False)):
    probs = PBD_to_probs(T, l1, l2, l3, m1, m2, solver_method=solver_method,
                         solver_kwargs=solver_kwargs)
    p, q = probs["pGS"], probs["pGE"]
    f = p + q
    dT = np.diff(T)
    dfdt = np.diff(f)/dT
    dqdt = np.diff(q) / dT
    rates = np.zeros((2, len(T))) # rates[0,:] = birth, rates[1,:] = death
    sumrates = dfdt / (1 - f[:-1])

    mu = dqdt + sumrates * q[:-1]
    l = sumrates - mu
    rates[0, :-1] = l
    rates[1, :-1] = mu

    # the last two points are unstable because of derivation twice
    rates[0, -2:] = rates[0, -3]
    rates[1, -2:] = rates[1, -3]
    # stabilize if convergence
    if smooth_convergence:
        rates = constant_after_conv(rates, **conv_crit)
    return np.flip(rates[0,:]), np.flip(rates[1,:]) # flip to have the present = last value

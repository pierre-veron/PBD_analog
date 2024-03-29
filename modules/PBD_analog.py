# Implements the PBD analogy and simulations 

import numpy as np

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

def approx_expected_T(l1, l2, l3, m1, m2):
    return 1/((1-pi(l1, l2, l3, m1, m2)) * l1 + m1) + tau(l1, l2, l3, m1, m2)

def analog_BD_rates(l1, l2, l3, m1, m2):
    p_speciation_ = p_speciation(l1, l2, l3, m1, m2)
    approx_expected_T_ = approx_expected_T(l1, l2, l3, m1, m2)
    spe_rate_ = p_speciation_ / approx_expected_T_
    ext_rate_ = (1-p_speciation_) / approx_expected_T_
    return spe_rate_, ext_rate_
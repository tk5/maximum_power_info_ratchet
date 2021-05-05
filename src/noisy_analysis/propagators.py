#!/usr/bin/env python3
#@author: jlucero
#date created: Fri Mar 19 22:08:03 PDT 2021
# purpose: define propagators needed for analysis of noisy system

from numpy import (
    pi, exp, sqrt, sinh, sign, abs, logical_and, where, finfo,
    linspace, logspace, log, cosh, tanh
    )

from scipy.integrate import trapz
from scipy.special import erf, erfc, xlogy
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import eigs
from scipy.optimize import minimize_scalar

from sys import stderr

def N_distr(x, mu, s2):
    """Define normal distribution of x with mean mu and variance s2"""
    return exp(-0.5*((x-mu)**2)/s2)/sqrt(2.0*pi*s2)

def compute_T(xkpr, u, dg=0.8, sg=0.1, t=1e-3, alpha=2):
    return (
            (exp(t - (-dg + (-1 + alpha)*u[None,:] + exp(t)*(dg + xkpr[:,None]))**2/(2.*(-1 + exp(2*t) + (alpha**2)*(sg**2))))
                *(1 + erf(((-1 + exp(2*t) + alpha*(sg**2))*u[None,:] + alpha*(sg**2)*(dg - exp(t)*(dg + xkpr[:,None])))/(sqrt(2)*sg*sqrt((-1 + exp(2*t))*(-1 + exp(2*t) + (alpha**2)*(sg**2)))))))
            /(2.*sqrt(2*pi)*sqrt(-1 + exp(2*t) + (alpha**2)*(sg**2)))
            + (exp(t/2. - ((dg + u[None,:] - exp(t)*(dg + xkpr[:,None]))**2/sinh(t))/(4.*exp(t)))*(1 - erf(u[None,:]/(sqrt(2)*sg))))/(4.*sqrt(pi)*sqrt(sinh(t)))
        )*abs(xkpr[1]-xkpr[0])

def compute_Ttilde(xkr, v, dg=0.8, sig=0.1, t=1e-3, alpha=2):
    return (
        -(exp(t/2. - ((dg + v[None,:] - exp(t)*(dg + xkr[:,None]))**2/sinh(t))/(4.*exp(t)))*(-1 + erf(xkr[:,None]/(sqrt(2)*sig))))/(4.*sqrt(pi)*sqrt(sinh(t)))
        + (exp((4*t - (2*xkr[:,None]**2)/((alpha**2)*(sig**2)) - ((dg - dg*exp(t) + v[None,:])**2/sinh(t))/exp(t) + ((1.0/sinh(t))*((alpha**2)*(sig**2)*(dg - dg*exp(t) + v[None,:]) - 2*(-1 + alpha)*xkr[:,None]*sinh(t))**2)/((alpha**2)*(sig**2)*((alpha**2)*(sig**2)*cosh(t) + (2 + alpha*(-4 + alpha*(2 + (sig**2))))*sinh(t))))/4.)*(-1 + (1.0/tanh(t)))*
        ((1.0/sinh(t))*(-2 + exp(2*t)*(1 + erf((-((alpha**2)*exp(t)*(sig**2)*(dg*(-1 + exp(t)) - v[None,:])) + (-1 + alpha)*(1 - exp(2*t))*xkr[:,None])/(sqrt(2)*alpha*sig*sqrt((-1 + exp(2*t))*(-(-1 + alpha)**2 + exp(2*t)*(1 + alpha*(-2 + alpha + alpha*(sig**2)))))))) +
        erfc((exp(t)*((alpha**2)*(sig**2)*(dg - dg*exp(t) + v[None,:]) - 2*(-1 + alpha)*xkr[:,None]*sinh(t)))/(sqrt(2)*alpha*sig*sqrt((-1 + exp(2*t))*(-(-1 + alpha)**2 + exp(2*t)*(1 + alpha*(-2 + alpha + alpha*(sig**2)))))))) +
        2*exp(t)*erf((abs((alpha**2)*exp(t)*(sig**2)*(dg*(-1 + exp(t)) - v[None,:]) + (-1 + alpha)*(-1 + exp(2*t))*xkr[:,None])*sqrt(2 + alpha*(-4 + alpha*(2 + (sig**2))) + (alpha**2)*(sig**2)*(1.0/tanh(t))))/(2.*alpha*sig*(-(-1 + alpha)**2 + exp(2*t)*(1 + alpha*(-2 + alpha + alpha*(sig**2))))))*
        sign((alpha**2)*exp(t)*(sig**2)*(dg*(-1 + exp(t)) - v[None,:]) + (-1 + alpha)*(-1 + exp(2*t))*xkr[:,None]) - 2*exp(t)*erf((abs(alpha*exp(t)*(sig**2)*(dg*(-1 + exp(t)) - v[None,:]) + xkr[:,None] - alpha*xkr[:,None] + exp(2*t)*(-1 + alpha + alpha*(sig**2))*xkr[:,None])*sqrt(2 + alpha*(-4 + alpha*(2 + (sig**2))) + (alpha**2)*(sig**2)*(1.0/tanh(t))))/
        (2.*exp(t)*((alpha**2)*exp(t)*(sig**3) + 2*(-1 + alpha)**2*sig*sinh(t))))*sign(alpha*exp(t)*(sig**2)*(dg*(-1 + exp(t)) - v[None,:]) + xkr[:,None] - alpha*xkr[:,None] + exp(2*t)*(-1 + alpha + alpha*(sig**2))*xkr[:,None]))*sinh(t))/(4.*sqrt(-2*(-1 + alpha)**2*pi + 2*exp(2*t)*pi*(1 + alpha*(-2 + alpha + alpha*(sig**2)))))
    )*abs(xkr[1]-xkr[0])

def find_steady_states(out_grid, in_grid, dg=0.8, sg=0.1, dt=1e-3, alpha=2):

    # compute transition matrices
    T = compute_T(out_grid, in_grid, dg, sg, dt, alpha)
    Ttilde = compute_Ttilde(out_grid, in_grid, dg, sg, dt, alpha)

    # turn matrices into sparse format
    sp_T = csc_matrix(T)
    sp_Ttilde = csc_matrix(Ttilde)

    # find the 3 largest eigenvalues and associated eigenvectors
    wT_ss, vT_ss = eigs(sp_T, k=3)
    wTtilde_ss, vTtilde_ss = eigs(sp_Ttilde, k=3)

    # find the eigenvector with eigenvalue 1
    p_xkpr = vT_ss[:, where((wT_ss - 1.0).__abs__() < finfo("float32").eps)[0][0]]
    p_xkr = vTtilde_ss[:, where((wTtilde_ss - 1.0).__abs__() < finfo("float32").eps)[0][0]]

    # re-normalize the eigenvectors to make them into distributions
    p_xkpr /= trapz(p_xkpr, out_grid)
    p_xkr /= trapz(p_xkr, out_grid)

    return p_xkpr.real, p_xkr.real


def compute_thermo_quants(ngrid=4000, dg=0.8, sn=0.1, ts=1e-3, alpha=2):

    LOW_LIMIT, UP_LIMIT = -60.0, 60.0
    out_grid = in_grid = linspace(LOW_LIMIT, UP_LIMIT, int(ngrid))
    p_xkpr, p_xkr = find_steady_states(out_grid, in_grid, dg, sn, ts, alpha)

    # regularization
    p_xkr[logical_and(p_xkr > -finfo("float32").eps, p_xkr < 0.0)] = 0.0
    p_xkpr[logical_and(p_xkpr > -finfo("float32").eps, p_xkpr < 0.0)] = 0.0

    # checks on distribution
    assert (p_xkr >= 0.0).all(), "p_xkr has non-positive entries!"
    assert (p_xkpr >= 0.0).all(), "p_xkpr has non-positive entries!"

    p_xkr_norm = trapz(p_xkr, out_grid)
    p_xkpr_norm = trapz(p_xkpr, out_grid)

    if (abs(p_xkr_norm - 1.0) > (finfo("float32").eps)):
        print(f"p_xkr not normalized! Normalization value {p_xkr_norm:.8f}", file=stderr)
    if (abs(p_xkpr_norm - 1.0) > (finfo("float32").eps)):
        print(f"p_xkpr not normalized! Normalization value {p_xkpr_norm:.8f}", file=stderr)

    # compute quantities
    mu_xkpr = trapz(out_grid*p_xkpr, out_grid)
    mu_xkr = trapz(out_grid*p_xkr, out_grid)
    ms_xkpr = trapz((out_grid**2)*p_xkpr, out_grid)
    ms_xkr = trapz((out_grid**2)*p_xkr, out_grid)

    Q_in = 0.5*(ms_xkpr-ms_xkr)
    W_in = 0.5*(ms_xkr-ms_xkpr)
    W_out = dg*(mu_xkpr-mu_xkr)

    return Q_in, W_in, W_out

def compute_info_flow(ngrid=4000, dg=0.8, sn=0.1, ts=1e-3, alpha=2):

    LOW_LIMIT, UP_LIMIT = -60.0, 60.0

    temp_grid = linspace(LOW_LIMIT, UP_LIMIT, int(ngrid))

    # ========== finding distributions ==========

    # find steady-state of relative coordinates
    p_xkpr, p_xkr = find_steady_states(temp_grid, temp_grid, dg, sn, ts, alpha)

    # regularization: zero out entries that are too small
    p_xkr[logical_and(p_xkr > -finfo("float32").eps, p_xkr < 0.0)] = 0.0
    p_xkpr[logical_and(p_xkpr > -finfo("float32").eps, p_xkpr < 0.0)] = 0.0

    # before proceeding to computations check that the distributions are behaving properly
    p_xkpr_norm = trapz(p_xkpr, temp_grid)
    p_xkr_norm = trapz(p_xkr, temp_grid)

    assert (p_xkpr >= 0.0).all(), "p_xkpr has non-positive entries!"
    assert (p_xkr >= 0.0).all(), "p_xkr has non-positive entries!"

    if (abs(p_xkpr_norm - 1.0) > (finfo("float32").eps)):
        print(f"p_xkpr not normalized! Normalization value {p_xkpr_norm:.8f}", file=stderr)
    if (abs(p_xkr_norm - 1.0) > (finfo("float32").eps)):
        print(f"p_xkr not normalized! Normalization value {p_xkr_norm:.8f}", file=stderr)

    # ========== computing entropies ==========

    # computing the conditional entropies
    H_xkr = -trapz(xlogy(p_xkr, p_xkr), temp_grid)
    H_xkpr = -trapz(xlogy(p_xkpr, p_xkpr), temp_grid)

    return -(H_xkpr-H_xkr)

def compute_transfer_entropy(sn=0.1, ts=1e-3):
    return 0.5*log((-1 + exp(2*ts) + (1 + exp(2*ts))*(sn**2)
        + exp(ts)*sqrt(((-1 + exp(2*ts))*(-1 + exp(2*ts) + 2*(1 + exp(2*ts))*(sn**2)
        + (-1 + exp(2*ts))*(sn**4)))/exp(2*ts)))/(2.0*exp(2*ts)*(sn**2))
        )

def optim(sn_in, ts_in, ngrid_in=4000):

    def objective_func(alpha, sn, ts, ngrid):
        _, win, _ = compute_thermo_quants(ngrid=ngrid, sn=sn, ts=ts, alpha=alpha)
        return abs(win)

    res = minimize_scalar(
        objective_func, bounds=(5e-3, 2.5),
        args=(sn_in, ts_in, ngrid_in), method="bounded"
        )

    return res.x, res.fun
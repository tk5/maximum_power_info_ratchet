#!/usr/bin/env python3

from numpy import (
    heaviside, exp as nexp, sqrt as nsqrt, pi as npi, finfo,
    linspace, logspace, where, abs as nabs
    )
from scipy.integrate import trapz

from scipy.sparse import csc_matrix
from scipy.sparse.linalg import eigs

from os import environ, mkdir
from os.path import isdir

from pltconfig import *


def N_distr(x, mu, sigma2):
    """Return the normal distribution with mean mu and variance sigma2."""
    return nexp(-0.5*((x-mu)**2)/sigma2)/nsqrt(2.0*npi*sigma2)


def compute_p_r_rr(r, rr, delta_g, dt):
    return (
        heaviside(-rr[None, :], 0.0)*N_distr(
            r[:, None], (rr[None, :]+delta_g)*nexp(-dt) - delta_g, 1.0-nexp(-2.0*dt)
            )
        + heaviside(rr[None, :], 0.0)*N_distr(
            r[:, None], -delta_g - (rr[None, :]-delta_g)*nexp(-dt), 1.0-nexp(-2.0*dt)
            )
        )*nabs(r[1]-r[0])


def compute_p_s_ss(s, ss, delta_g, dt):
    return heaviside(-s[:, None], 0.0)*(
        N_distr(s[:, None], delta_g - (ss[None, :]+delta_g)*nexp(-dt), 1.0-nexp(-2.0*dt))
        + N_distr(s[:, None], (ss[None, :]+delta_g)*nexp(-dt) - delta_g, 1.0-nexp(-2.0*dt))
    )*nabs(s[1]-s[0])


def find_steady_states(r, rr, delta_g, dt):

    # compute transition matrices
    pr_rr = compute_p_r_rr(r, rr, delta_g, dt)
    ps_ss = compute_p_s_ss(r, rr, delta_g, dt)

    # turn matrices into sparse format
    sp_pr_rr = csc_matrix(pr_rr)
    sp_ps_ss = csc_matrix(ps_ss)

    # find the 3 largest eigenvalues and associated eigenvectors
    wr_rr, vr_rr = eigs(sp_pr_rr, k=3)
    ws_ss, vs_ss = eigs(sp_ps_ss, k=3)

    # find the eigenvector with eigenvalue 1
    pr = vr_rr[:, where((wr_rr - 1.0).__abs__() < finfo("float32").eps)[0][0]]
    ps = vs_ss[:, where((ws_ss - 1.0).__abs__() < finfo("float32").eps)[0][0]]

    # re-normalize the eigenvectors to make them into distributions
    pr /= trapz(pr, r)
    ps /= trapz(ps, r)

    return pr, ps


def compute_means(target_dir):

    delta_g = 0.8
    low_limit = -20.0
    up_limit = 20.0
    ngrid = 2e3
    grid_from = linspace(low_limit, up_limit, int(ngrid))
    grid_to = linspace(low_limit, up_limit, int(ngrid))

    nscan = 1
    times = logspace(1.0, -3.0, nscan)

    target_file = f"/ref_deltag_{delta_g:.2f}_outfile.txt"

    with open(target_dir + target_file, "w") as ofile:
        for time in times:

            pr, ps = find_steady_states(grid_to, grid_from, delta_g, time)

            # compute the mean work
            mean_power_out = delta_g*(
                trapz(grid_to*pr.real, grid_to)
                - trapz(grid_to*ps.real, grid_to)
                )/time

            ofile.write(f"{time:.15e}\t{mean_power_out:.15e}\n")
            ofile.flush()


if __name__ == "__main__":

    target_repo = "/Users/jlucero/data_dir/"
    today_dir = environ["TODAY"]

    target_dir = target_repo + "/" + today_dir + "/"
    if not isdir(target_dir):
        print("Target directory doesn't exist. Making it now.")
        mkdir(target_dir)

    compute_means(target_dir)

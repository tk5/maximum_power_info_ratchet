#!/usr/bin/env python3
#@author: jlucero
#date created: Fri Mar 19 22:19:48 PDT 2021
# purpose: launcher for noisy system theoretical analysis

from numpy import linspace

from os import environ, makedirs
from os.path import isdir

from propagators import (
    compute_thermo_quants, compute_info_flow, optim, compute_transfer_entropy
    )

def get_params():

    ts = 1e-3
    sn_min = 1e-1
    sn_max = 1.0
    sn_scan = 1

    sn_array = linspace(sn_min, sn_max, int(sn_scan))

    return (ts, sn_array)

def main():

    ts, snvals = get_params()

    for i, sn in enumerate(snvals):

        # choose the value of astar = \alpha^{*} for input work = 0
        astar, _ = optim(sn, ts, ngrid_in=15000)

        # given value of astar compute the input heat, input work and output work,
        qin, win, wout = compute_thermo_quants(
            ngrid=15000, sn=sn, ts=ts, alpha=astar
            )
        # compute the information flow
        inf_flow = compute_info_flow(ngrid=15000, sn=sn, ts=ts, alpha=astar)
        # compute the transfer entropy
        trans_entr = compute_transfer_entropy(sn=sn, ts=ts)

        # compute the efficiencies using the two information-theoretic quantities
        eta_flow = wout/(win-inf_flow)
        eta_trans = wout/(win+trans_entr)

        # print results to file
        print(
            f"{sn:.3f}\t"                 # noise magnitude
            + f"{astar:.15e}\t"            # astar => input work = 0
            + f"{qin/ts:.15e}\t"           # input heat
            + f"{win/ts:.15e}\t"           # input (trap) work [minus input heat in steady-state]
            + f"{wout/ts:.15e}\t"          # output (gravitational) work
            + f"{inf_flow/ts:.15e}\t"      # information flow between x to \lambda
            + f"{trans_entr/ts:.15e}\t"    # transfer entropy between x to y
            + f"{eta_flow:.15e}\t"         # efficiency using information flow
            + f"{eta_trans:.15e}"          # efficiency using transfer entropy
        )

if __name__ == "__main__":
    main()

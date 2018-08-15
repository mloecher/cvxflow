import numpy as np
cimport numpy as np
import time

cdef extern from "gropt_proc.cpp":
    void _opt3 "opt3"(int N, double dt, double *G_in, double *G, double *resid, double *moments, int n_it, double cushion, double gmax, double smax, double tmax)
    void _get_stim "get_stim"(int N, double dt, double dt2, double c, double *G_in, double *stim_out)

def opt3(N, d_M0, d_M1, G_in = None, dt = 10.0e-3, n_it = 2000, cushion = 0.95, gmax = 80.0, smax = 200.0, tmax = 1.0):

    mm = []
    for i in range(3):
        mm.append(d_M0[i])
        mm.append(d_M1[i])

    m = np.array(mm)
    m = np.ascontiguousarray(np.ravel(m), np.float64)
    cdef np.ndarray[np.float64_t, ndim=1, mode="c"] m_c = m

    G = np.ascontiguousarray(np.ravel(np.zeros(N*3, dtype=np.float64)))
    cdef np.ndarray[np.float64_t, ndim=1, mode="c"] G_c = G

    if G_in is None:
        G_in = np.ascontiguousarray(np.ravel(np.zeros(N*3, dtype=np.float64)))
    cdef np.ndarray[np.float64_t, ndim=1, mode="c"] G_in_c = np.ascontiguousarray(np.ravel(G_in), np.float64)

    cdef double resid

    _opt3(N, dt, &G_in_c[0], &G_c[0], &resid, &m_c[0], n_it, cushion, gmax, smax, tmax)

    return G, resid

def get_stim(G_in, dt = 10.0e-3, n_it = 2000, cushion = 1.0, gmax = 80.0, smax = 200.0, tmax = 1.0):

    N = G_in.size//3

    cdef np.ndarray[np.float64_t, ndim=1, mode="c"] G_in_c = np.ascontiguousarray(np.ravel(G_in), np.float64)

    stim_out = np.ascontiguousarray(np.ravel(np.zeros(N*3, dtype=np.float64)))
    cdef np.ndarray[np.float64_t, ndim=1, mode="c"] stim_out_c = np.ascontiguousarray(np.ravel(stim_out), np.float64)

    c = 1.0

    _get_stim(N, dt, dt,c, &G_in_c[0], &stim_out_c[0])

    return stim_out
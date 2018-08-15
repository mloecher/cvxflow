import numpy as np
cimport numpy as np
import time

cdef extern from "gropt_proc_spoiler.cpp":
    void _opt3_spoil "opt3"(int N, double dt, double *G_in, double *G, double *resid, double *moments, int n_it, double cushion, double gmax, double smax, double tmax, double resx)


def opt3_spoil(N, d_M0, G_in = None, dt = 10.0e-3, n_it = 2000, cushion = 0.95, gmax = 80.0, smax = 200.0, tmax = 1.0):

    mm = []
    for i in range(3):
        mm.append(d_M0[i])

    m = np.array(mm)
    m = np.ascontiguousarray(np.ravel(m), np.float64)
    cdef np.ndarray[np.float64_t, ndim=1, mode="c"] m_c = m

    G = np.ascontiguousarray(np.ravel(np.zeros(N*3, dtype=np.float64)))
    cdef np.ndarray[np.float64_t, ndim=1, mode="c"] G_c = G

    if G_in is None:
        G_in = np.ascontiguousarray(np.ravel(np.zeros(N*3, dtype=np.float64)))
    cdef np.ndarray[np.float64_t, ndim=1, mode="c"] G_in_c = np.ascontiguousarray(np.ravel(G_in), np.float64)

    cdef double resid

    resx = d_M0[0]/11.830749999999998

    _opt3_spoil(N, dt, &G_in_c[0], &G_c[0], &resid, &m_c[0], n_it, cushion, gmax, smax, tmax, resx)

    return G, resid

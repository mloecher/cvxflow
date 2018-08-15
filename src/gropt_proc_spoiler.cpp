#include <cmath>
#include <iostream>
#include <fstream>
#include "Eigen/Core"
#include "Eigen/Sparse"

// #include <sys/time.h>
// typedef unsigned long long timestamp_t;

// static timestamp_t
// get_timestamp ()
// {
//   struct timeval now;
//   gettimeofday (&now, NULL);
//   return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
// }

using namespace Eigen;

class StimOpt {
    public:    
        StimOpt(int init_N, double init_dt);
        
        VectorXd forward(const VectorXd &G);
        VectorXd transpose(const VectorXd &G);
        VectorXd clip_stim(const VectorXd &x, double tmax);
        
        int N;
        float dt;

        float t1x, t2x, t3x;
        float t1y, t2y, t3y;
        float t1z, t2z, t3z;

        float a1x, a2x, a3x;
        float a1y, a2y, a3y;
        float a1z, a2z, a3z;
        

        MatrixXd S;

        MatrixXd M1x;
        MatrixXd M1y;
        MatrixXd M1z;

        MatrixXd M2x;
        MatrixXd M2y;
        MatrixXd M2z;

        MatrixXd M3x;
        MatrixXd M3y;
        MatrixXd M3z;

        MatrixXd tS3;
        MatrixXd tM13;
        MatrixXd tM23;
        MatrixXd tM33;

        SparseMatrix<double> S3;
        SparseMatrix<double> M13;
        SparseMatrix<double> M23;
        SparseMatrix<double> M33;

        SparseMatrix<double> T;
};

StimOpt::StimOpt(int init_N, double init_dt)
{
    N = init_N;
    dt = init_dt;

    t1x = 1.88401e-01;
    t2x = 1.19991e+01;
    t3x = 8.48963e-01;
    t1y = 2.32399e-01;
    t2y = 1.20000e+01;
    t3y = 9.40258e-01;
    t1z = 1.60200e-01;
    t2z = 1.20014e+01;
    t3z = 7.43985e-01;
    a1x = 2.19094e-02;
    a2x = 1.53346e-04;
    a3x = 3.14602e-03;
    a1y = 2.12378e-02;
    a2y = 3.04470e-04;
    a3y = 5.68286e-03;
    a1z = 2.68527e-02;
    a2z = 1.46336e-04;
    a3z = 2.91029e-03;

    S.resize(N,N);
    S.setZero();

    for (int i = 0; i < N-1; i++) {
        S(i, i) = -1.0;
        S(i, i+1) = 1.0;
    }

    // -------
    M1x.resize(N,N);
    M1x.setZero();
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < j; i++) {
            M1x(j, i) = a1x * exp(-(j-i-1) * dt / t1x);
        }
    }

    M1y.resize(N,N);
    M1y.setZero();
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < j; i++) {
            M1y(j, i) = a1y * exp(-(j-i-1) * dt / t1y);
        }
    }

    M1z.resize(N,N);
    M1z.setZero();
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < j; i++) {
            M1z(j, i) = a1z * exp(-(j-i-1) * dt / t1z);
        }
    }

    // -------
    M2x.resize(N,N);
    M2x.setZero();
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < j; i++) {
            M2x(j, i) = a2x * exp(-(j-i-1) * dt / t2x);
        }
    }

    M2y.resize(N,N);
    M2y.setZero();
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < j; i++) {
            M2y(j, i) = a2y * exp(-(j-i-1) * dt / t2y);
        }
    }

    M2z.resize(N,N);
    M2z.setZero();
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < j; i++) {
            M2z(j, i) = a2z * exp(-(j-i-1) * dt / t2z);
        }
    }

    // -------
    M3x.resize(N,N);
    M3x.setZero();
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < j; i++) {
            M3x(j, i) = a3x * exp(-(j-i-1) * dt / t3x);
        }
    }

    M3y.resize(N,N);
    M3y.setZero();
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < j; i++) {
            M3y(j, i) = a3y * exp(-(j-i-1) * dt / t3y);
        }
    }

    M3z.resize(N,N);
    M3z.setZero();
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < j; i++) {
            M3z(j, i) = a3z * exp(-(j-i-1) * dt / t3z);
        }
    }

    // -------
    tS3.resize(N*3, N*3);
    tS3.setZero();
    tS3.block(0 * N, 0 * N, N, N) = S;
    tS3.block(1 * N, 1 * N, N, N) = S;
    tS3.block(2 * N, 2 * N, N, N) = S;

    tM13.resize(N*3, N*3);
    tM13.setZero();
    tM13.block(0 * N, 0 * N, N, N) = M1x;
    tM13.block(1 * N, 1 * N, N, N) = M1y;
    tM13.block(2 * N, 2 * N, N, N) = M1z;

    tM23.resize(N*3, N*3);
    tM23.setZero();
    tM23.block(0 * N, 0 * N, N, N) = M2x;
    tM23.block(1 * N, 1 * N, N, N) = M2y;
    tM23.block(2 * N, 2 * N, N, N) = M2z;

    tM33.resize(N*3, N*3);
    tM33.setZero();
    tM33.block(0 * N, 0 * N, N, N) = M3x;
    tM33.block(1 * N, 1 * N, N, N) = M3y;
    tM33.block(2 * N, 2 * N, N, N) = M3z;

    S3 = tS3.sparseView(1.0, 1.0e-6);
    M13 = tM13.sparseView(1.0, 1.0e-6);
    M23 = tM23.sparseView(1.0, 1.0e-6);
    M33 = tM33.sparseView(1.0, 1.0e-6);

    // -------
    T = M13*S3 + M23*S3 + M33*S3;
}

VectorXd StimOpt::forward(const VectorXd &G)
{

    VectorXd out;

    VectorXd S3G = S3*G;

    VectorXd S3G_sign(S3G.size());

    for (int i = 0; i < S3G.size(); i++) {
        if (S3G(i)<0) {S3G_sign(i) = -1.0;}
        else {S3G_sign(i) = 1.0;}
    }
    
    out = ( (M13*(S3G)) ).array() + 
          ( M23*((S3G).cwiseAbs()) ).array() * S3G_sign.array() + 
          ( (M33*(S3G)) ).array();
    

    return out;
}

VectorXd StimOpt::transpose(const VectorXd &G)
{

    VectorXd out;

    VectorXd S3G = S3.transpose()*G;
    
    VectorXd S3G_sign(S3G.size());

    for (int i = 0; i < S3G.size(); i++) {
        if (S3G(i)<0) {S3G_sign(i) = -1.0;}
        else {S3G_sign(i) = 1.0;}
    }

    out = ( (M13.transpose()*(S3G)) ).array() + 
          ( M23.transpose()*((S3G).cwiseAbs()) ).array() * S3G_sign.array() + 
          ( (M33.transpose()*(S3G)) ).array();
    

    return out;
}

VectorXd StimOpt::clip_stim(const VectorXd &x, double tmax)
{
    VectorXd out(x.size());

    double s_norm;
    for (int i = 0; i < N; i++) {
        s_norm = sqrt(x(i)*x(i) + x(i + N)*x(i + N) + x(i + 2*N)*x(i + 2*N));
        if (s_norm > tmax) {
            out(i) = tmax * x(i) / s_norm;
            out(i + N) = tmax * x(i + N) / s_norm;
            out(i + 2*N) = tmax * x(i + 2*N) / s_norm;
        } else {
            out(i) = x(i);
            out(i + N) = x(i + N);
            out(i + 2*N) = x(i + 2*N);
        }
    }

    return out;
}

MatrixXd make_moment_mat(int N, double dt, int N_moments)
{

    VectorXd tvec(N);

    for (int i = 0; i < N; i++) {
        tvec(i) = i * dt;
    }

    MatrixXd M(N_moments, N);

    for (int i = 0; i < N_moments; i++) {
        for (int j = 0; j < N; j++) {
            M(i, j) = dt * pow(tvec(j), i);
        }
    }
    MatrixXd M3 = MatrixXd::Zero(M.rows() * 3, M.cols() * 3);
    for (int i = 0; i < 3; i++)
    {
        M3.block(i * M.rows(), i * M.cols(), M.rows(), M.cols()) = M;
    }
    return M3;
}

VectorXd calc_M_norms(MatrixXd *M)
{
    VectorXd out((*M).rows());
    out.setZero();
    
    for(int i = 0; i < (*M).rows(); i++) {
        for (int j = 0; j < (*M).cols(); j++) {
            out(i) += (*M)(i,j) * (*M)(i,j);
        }
    }

    for(int i = 0; i < (*M).rows(); i++) {
        out(i) = sqrt(out(i));
    }

    for(int i = 0; i < (*M).rows(); i++) {
        for (int j = 0; j < (*M).cols(); j++) {
            (*M)(i,j) /= out(i);
        }
    }

    return out;
}

SparseMatrix<double> make_slew_mat(int N)
{
    SparseMatrix<double> S(N*3, N*3); 

    for (int i = 0; i < N-1; i++) {
        S.insert(i, i) = -1.0;
        S.insert(i, i+1) = 1.0;

        S.insert(i + N, i + N) = -1.0;
        S.insert(i + N, i + N + 1) = 1.0;

        S.insert(i + 2*N, i + 2*N) = -1.0;
        S.insert(i + 2*N, i + 2*N + 1) = 1.0;
    }

    return S;

}

VectorXd clip_gmax(const VectorXd &x, int N, double gmax, double dt, double resx)
{
    VectorXd out(x.size());
    for (int i = 0; i < x.size(); i++) {
        if (x(i) > gmax) {out(i) = gmax;}
        else if (x(i) < -gmax) {out(i) = -gmax;}
        else {out(i) = x(i);}
    }

    out(0) = 13.1069/resx;
    out(N-1) = 0.0;

    out(N) = 0.0;
    out(2*N-1) = 0.0;

    out(2*N) = 0.0;
    out(3*N-1) = 0.0;


    return out;
}

VectorXd clip_smax(const VectorXd &x, double smax)
{
    VectorXd out(x.size());

    for (int i = 0; i < x.size(); i++) {
        if (x(i) > smax) {out(i) = smax;}
        else if (x(i) < -smax) {out(i) = -smax;}
        else {out(i) = x(i);}
    }

    return out;
}

double calc_resid(const VectorXd &G, const MatrixXd &M, const VectorXd &moments, const VectorXd &normM, const MatrixXd &S, StimOpt &stimopt, int N, double dt, double gmax, double smax, double tmax)
{
    double resid = 0.0;

    for (int i = 0; i < G.size(); i++) {
        if ( fabs(G(i)) > gmax ) {
            resid = 10.0e16;
            // std::cout << "gmax violation" << std::endl;  
        }
    }

    VectorXd SG = S*G;
    for (int i = 0; i < SG.size(); i++) {
        if ( fabs(SG(i)) > dt * smax ) {
            resid = 10.0e16;
            // std::cout << "smax violation" << std::endl;
        }
    }

    VectorXd TG = stimopt.forward(G);
    double s_norm;
    for (int i = 0; i < N; i++) {
        s_norm = sqrt(TG(i)*TG(i) + TG(i + N)*TG(i + N) + TG(i + 2*N)*TG(i + 2*N));
        if (s_norm > tmax) {
            resid = 10.0e16;
            // std::cout << "tmax violation" << std::endl;
        }
    }

    VectorXd dM = M*G - moments;
    dM.array() *= normM.array();
    resid += dM.norm();

    return resid;
}

void opt3(int N, double dt,  double *G_in, double *G_out, double *resid, double *moments_in, int n_iter, double cushion, double gmax, double smax, double tmax, double resx)
{
    VectorXd moments = VectorXd::Zero(3);
    for (int i = 0; i < 3; i++) {
        moments(i) = moments_in[i];
    }

    MatrixXd M = make_moment_mat(N, dt, 1);
    VectorXd normM = calc_M_norms(&M);
    moments.array() /= normM.array();
    VectorXd M_colsum = M.cwiseAbs().colwise().sum();
    SparseMatrix<double> S = make_slew_mat(N);    
    VectorXd S_colsum = RowVectorXd::Ones(S.rows()) * S.cwiseAbs();

    StimOpt stimopt(N, dt);
    MatrixXd T = stimopt.T;
    VectorXd T_colsum = T.cwiseAbs().colwise().sum();

    VectorXd tau = (M_colsum + S_colsum + T_colsum).cwiseInverse();

    VectorXd sigM = M.cwiseAbs().rowwise().sum();
    sigM = 10.0 * sigM.cwiseInverse();

    VectorXd sigS = S.cwiseAbs() * VectorXd::Ones(S.cols());
    sigS = sigS.cwiseInverse();

    VectorXd sigT = T.cwiseAbs().rowwise().sum();
    sigT = 1.0 * sigT.cwiseInverse();
    for (int i = 0; i < sigT.size(); i++) {
        if ( std::isinf(sigT(i)) ) {sigT(i) = 1.0;}
    }

    VectorXd G = VectorXd::Zero(N*3);

    for (int i = 0; i < N*3; i++) {
        G(i) = G_in[i];
    }

    VectorXd xbar = G;

    VectorXd zM = VectorXd::Zero(M.rows());
    VectorXd zS = VectorXd::Zero(S.rows());
    VectorXd zT = VectorXd::Zero(T.rows());
    VectorXd zK = VectorXd::Zero(G.rows());

    VectorXd txmx;
    VectorXd zMbuff;
    VectorXd zSbuff;
    VectorXd zTbuff;
    VectorXd zKbuff;

    VectorXd zMbar;
    VectorXd zSbar;
    VectorXd zTbar;
    VectorXd zKbar;

    int count = 0;

    double rr = 0;

    while (count < n_iter) {
        xbar = G.array() - tau.array() * (M.transpose()*zM + S.transpose()*zS + stimopt.transpose(zT)).array();


        xbar = clip_gmax(xbar, N, 0.99*cushion*gmax, dt, resx);

        txmx    = 2.0*xbar-G;
        zMbuff  = zM.array() + sigM.array()*(M*txmx).array();
        zSbuff  = zS.array() + sigS.array()*(S*txmx).array();
        zTbuff  = zT.array() + sigT.array()*(stimopt.forward(txmx)).array();

        zMbar = zMbuff.array() - sigM.array()*moments.array();
        zSbar = zSbuff.array() - sigS.array()*clip_smax(zSbuff.array()/sigS.array(), 0.99*cushion*dt*smax).array();
        zTbar = zTbuff.array() - sigT.array()*(stimopt.clip_stim(zTbuff.array()/sigT.array(), 0.99*cushion*tmax)).array();

        G = xbar;
        zM = zMbar;
        zS = zSbar; 
        zT = zTbar;

        count++;

        if (count % 20 == 0) {
            rr = calc_resid(G, M, moments, normM, S, stimopt, N, dt, cushion*gmax, cushion*smax, cushion*tmax);
            if (rr < 1.0e-3) {break;}
        }
    }

    for (int i = 0; i < N*3; i++) {
        G_out[i] = G(i);
    }

    *resid = calc_resid(G, M, moments, normM, S, stimopt, N, dt, cushion*gmax, cushion*smax, cushion*tmax);

}


int main(int argc, char *argv[])
{
    std::cout << " --- Standalone Execution --- " << argc << std::endl;

   

    float resx;
    float resy;
    float resz;
    float venc;
    char prot_type[32];
    int Ny = 160;
    int Nz = 32;

    int mode; 
    if (argc > 6) {
        mode = atoi(argv[6]);
    } else {
        mode = 1;
    }
    std::cout << "mode = " << mode << std::endl;


    if (mode == 1) {
        resx = 2.0;
        resy = 2.0;
        resz = 3.0;
        venc = 160.0;
        Ny = 160;
        Nz = 32;
        sprintf(prot_type, "cardiac");
    } else if (mode == 2) {
        resx = 1.0;
        resy = 1.0;
        resz = 2.0;
        venc = 80.0;
        Ny = 168;
        Nz = 32;
        sprintf(prot_type, "neuro");
    } else if (mode == 3) {
        resx = 2.2;
        resy = 2.2;
        resz = 3.0;
        venc = 160.0;
        Ny = 160;
        Nz = 32;
        sprintf(prot_type, "cardiac22");
    } else if (mode == 4) {
        resx = 3.0;
        resy = 3.0;
        resz = 3.0;
        venc = 160.0;
        Ny = 160;
        Nz = 32;
        sprintf(prot_type, "ph160");
    } else if (mode == 5) {
        resx = 3.0;
        resy = 3.0;
        resz = 3.0;
        venc = 80.0;
        Ny = 160;
        Nz = 32;
        sprintf(prot_type, "ph80");
    } else if (mode == 6) {
        resx = 0.75;
        resy = 0.75;
        resz = 0.75;
        venc = 80.0;
        Ny = 160;
        Nz = 32;
        sprintf(prot_type, "timing75");
    } else if (mode == 7) {
        resx = 1.0;
        resy = 1.0;
        resz = 1.0;
        venc = 80.0;
        Ny = 160;
        Nz = 32;
        sprintf(prot_type, "timing100");
    } else if (mode == 8) {
        resx = 1.5;
        resy = 1.5;
        resz = 1.5;
        venc = 80.0;
        Ny = 160;
        Nz = 32;
        sprintf(prot_type, "timing150");
    } else if (mode == 9) {
        resx = 2.0;
        resy = 2.0;
        resz = 2.0;
        venc = 80.0;
        Ny = 160;
        Nz = 32;
        sprintf(prot_type, "timing200");
    }
    
    double ss_m0;
    if (argc > 1) {
        ss_m0 = atof(argv[1]);
    } else {
        ss_m0 = 0.0;
    }
    std::cout << "ss_m0 = " << ss_m0 << std::endl;

    float cushion;
    if (argc > 2) {
        cushion = atof(argv[2]);
    } else {
        cushion = 0.9f;
    }
    std::cout << "cushion = " << cushion << std::endl;

    float tmax;
    if (argc > 3) {
        tmax = atof(argv[3]);
    } else {
        tmax = 0.7f;
    }
    std::cout << "tmax = " << tmax << std::endl;

    double gmax; 
    if (argc > 4) {
        gmax = atof(argv[4]);
    } else {
        gmax = 45.0f;
    }
    std::cout << "gmax = " << gmax << std::endl;

    int N_in; 
    if (argc > 5) {
        N_in = atoi(argv[5]);
    } else {
        N_in = 8;
    }
    if (N_in < 0) {
        N_in = 8;
    }
    std::cout << "N_in = " << N_in << std::endl;

    char filename [100];
    int n;
    n = sprintf (filename, "%s_spoil_%.1f_%d_%d_%d.raw", prot_type, ss_m0, (int)(cushion*100), (int)(tmax*100), (int)gmax);
    std::ofstream myFile (filename, std::ios::out | std::ios::binary);

    std::cout << filename << std::endl;

    float moments[3];
    moments[0] = 11.74 * 1.0 / resx;
    moments[1] = -11.74 * 1.0 / resy;
    moments[2] = -11.74 * 1.0 / resz;



    std::cout << "moments = " << moments[0] << " " << moments[1] << " " << moments[2] << std::endl;
    
    int n_iter = 3000;
    double smax = 200.0; 
    double resid = 0.0;
    double dt = 40.0e-3;

    std::cout << "gmax = " << gmax << std::endl;


    double line_c = 0.0;
    double par_c = 0.0;


    int N = (int)N_in + 1;

    double* G_in;
    G_in = new double [1000];

    double* G_out;
    G_out = new double [1000];

    double* G_0;
    G_0 = new double [1000];
    for (int i = 0; i < N*3; i++) {
        G_0[i] = 0.0;
    }

    double* G_1;
    G_1 = new double [1000];
    for (int i = 0; i < N*3; i++) {
        G_1[i] = 0.0;
    }

    double moments_in[3];
    Matrix<double,3,1>  lin_par;
    Matrix<double,3,1>  d_M0;
    Matrix<double,3,1>  init_m0; 
    init_m0 << moments[0], moments[1], moments[2];

    bool error_flag = false;



    float Nf = (float) (4.0 * (N-1)) + 1.0;

    float* header_vals = new float [32];

    header_vals[0] = Nf;
    header_vals[1] = Ny;
    header_vals[2] = Nz;

    
    
    float Ny2 = (double) Ny / 2.0;
    float Nz2 = (double) Nz / 2.0;

    bool still_working = true;

    while (still_working) {
        std::cout << "Testing N = " << N << std::endl;

        Nf = (float) (4.0 * (N-1)) + 1.0;

        header_vals[0] = Nf;
        header_vals[1] = Ny;
        header_vals[2] = Nz;

        myFile.seekp (0);

        for (int i = 0; i < 32; i++) {
            float write_temp = header_vals[i];
            myFile.write ((char*) &write_temp, sizeof(float));
        }

        for (int i = 0; i < N*3; i++) {
            G_0[i] = 0.0;
        }

        for (int i = 0; i < N*3; i++) {
            G_1[i] = 0.0;
        }

        error_flag = false;

        for (int part = 0; part < Nz; part++) {
            if (error_flag) break;
            for (int line = 0; line < Ny; line++) {
                if (error_flag) break;

                line_c = -(Ny2-line) / Ny2;
                par_c = -(Nz2-part) / Nz2;
                
                lin_par << 1.0, line_c, par_c;
                d_M0 = init_m0.array() * lin_par.array();


                if (line == 0) {
                    for (int i = 0; i < N*3; i++) {
                        G_in[i] = G_0[i];
                    }
                } else {
                    for (int i = 0; i < N*3; i++) {
                        G_in[i] = G_1[i];
                    }
                }
                
                moments_in[0] = d_M0(0); 
                moments_in[1] = d_M0(1); 
                moments_in[2] = d_M0(2); 

                opt3(N, dt, G_in, G_out, &resid, moments_in, n_iter, cushion, gmax, smax, tmax, resx);
                
                if (resid > .01) {
                    std::cout << std::endl << "ERROR: Line, Part: " << line << ", " << part << " resid too high = " << resid << std::endl << std::endl;

                    std::cout << dt << " " << n_iter << " " << cushion << " " << gmax << " " << smax << " " << tmax << " " << std::endl;

                    error_flag = true;

                    N += 1;
                    break;
                }

                if (line == 0) {
                    for (int i = 0; i < N*3; i++) {
                        G_0[i] = G_out[i];
                    }
                } 

                for (int i = 0; i < N*3; i++) {
                    G_1[i] = G_out[i];
                }

                for (int j = 0; j < 3; j++) {
                    for (int i = 0; i < N-1; i++) {
                        float gint0 = G_out[i + j*N];
                        float gint1 = G_out[i + j*N + 1];
                        for (int k = 0; k < 4; k++) {
                            float write_temp = (float)  (gint0 * (4.0-k)/4.0 + gint1 * ((float)k)/4.0);
                            myFile.write ((char*) & write_temp, sizeof(float));
                        }
                    }
                    float write_temp = (float) G_out[N-1 + j*N];
                    myFile.write ((char*) & write_temp, sizeof(float));
                }
            }
            if ((part % 10) == 0) {
                std::cout << d_M0.transpose() << std::endl;
                std::cout << N << "  Part = " << part << "   resid = " << resid << std::endl << std::endl;
            }
        }
        if (error_flag == false) {still_working = false;}
    }

    myFile.close();

    return 0;
}
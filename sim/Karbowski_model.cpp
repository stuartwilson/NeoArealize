#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <math.h>

using namespace std;

double randFloat (void)
{
    return ((double) rand()) / (double)RAND_MAX;
}

vector<double> space (vector<double> x, double dx)
{
    int n = x.size();
    vector<double> X(1,x[0]);

    for (int xi=0; xi<n; xi++) {
        X.push_back (x[xi]);
    }
    X.push_back(x[n-1]);

    for (int xi=0; xi<n; xi++) {
        x[xi] = (X[xi+2]-X[xi]) * dx;
    }

    return x;
};

// Does: f = (alpha * ci) + betaterm. c.f. Karb2004, Eq 1
vector<double> compute_dci_dt (vector<double> ci, double alpha, vector<double> betaterm)
{
    vector<double> dci_dt;
    for (unsigned int xi=0; xi<ci.size(); xi++) {
        dci_dt.push_back (alpha * ci[xi] + betaterm[xi]);
    }
    return dci_dt;
};

// Computes the "flux of axonal branches" term, and adds "addon",
// which is alpha * ci - betaterm. Computes dJdX, then multiplies by
// 2.dx as the sum is carried out over two distance increments.
vector<double> compute_axonalbranchflux (vector<double> a,
                                         vector<double> g,
                                         double D,
                                         double dx,
                                         vector<double> addon)
{
    int n = a.size();
    double dJdX;
    vector<double> output(n,0.);
    for (int j=0; j<n; j++) {
        if(j==0){
            dJdX =
                D*(a[j+2]-a[j])     -a[j+1]*g[j+1];
        }
        else if(j==1) {
            dJdX =
                D*(a[j+2]-a[j])     -a[j+1]*g[j+1]
                -D*(a[j]-a[j-1])    +a[j-1]*g[j-1];
        }
        else if(j==n-1){
            dJdX =
                -D*(a[j]-a[j-1])    +a[j]*g[j];
        }
        else if (j==n-2) {
            dJdX =
                D*(a[j+1]-a[j])     -a[j+1]*g[j+1]
                -D*(a[j]-a[j-2])    +a[j-1]*g[j-1];
        }
        else {
            dJdX =
                D*(a[j+2]-a[j])     -a[j+1]*g[j+1]
                -D*(a[j]-a[j-2])    +a[j-1]*g[j-1];
        }
        output[j] = dJdX*2.*dx + addon[j];
    }
    return output;
}


int main (int argc, char **argv)
{
    if (argc<2) {
        cout << "supply steps and filename" << endl;
        return 0;
    }

    // PARAMETERS
    double dt = 0.1;
    double halfdt = 0.5*dt, sixthdt = dt/6.0;

    // Number of TC axon types
    int axontypes = 5;

    double Aemx = 1.34, Apax = 1.4, Afgf = 0.9;
    double Cemx = 25.6, Cpax = 27.3, Cfgf = 26.4;
    double w1 = 2.4, w2 = 2.1, v1 = 2.6, v2 = 2.7;

    vector<double> gammaA(axontypes,0.), gammaB(axontypes,0.), gammaC(axontypes,0.);
    gammaA[0]=1.6, gammaA[1]=-0.4, gammaA[2]=-2.21, gammaA[3]=-2.1, gammaA[4]=-2.45;
    gammaB[0]=-0.6, gammaB[1]=-0.5, gammaB[2]=0.4, gammaB[3]=-0.5, gammaB[4]=-1.0;
    gammaC[0]=-2.9, gammaC[1]=-2.5, gammaC[2]=-2.23, gammaC[3]=-0.6, gammaC[4]=1.7;

    double theta1 = 0.77, theta2 = 0.5, theta3 = 0.39, theta4 = 0.08;
    double kA = 0.58, kB = 0.9, kC = 0.55;
    double sigmaA = 0.2, sigmaB = 0.2, sigmaC = 0.2;
    double k = 3., D = 0.1, L = 40., dx = 0.25, alpha = 3.0, beta = 3.0;


    // INIT
    int nx = int((L+dx)/dx);

    vector<double> eta_emx(nx,0.), eta_pax(nx,0.), eta_fgf(nx,0.);
    for (int i=0; i<nx; i++) {
        double x = (double)i*((double)L/(double)(nx-1));
        eta_emx[i] = Aemx * exp (-((x-L) * (x-L)) / (Cemx * Cemx));
        eta_pax[i] = Apax * exp (-(x*x) / (Cpax * Cpax));
        eta_fgf[i] = Afgf * exp (-(x*x) / (Cfgf * Cfgf));
    }

    // Precompute the density of the guidance molecules.
    vector<double> s(nx,0.), r(nx,0.), f(nx,0.);
    double taus = 0.0001, taur = 0.0001, tauf = 0.0001;
    for(int t=0;t<300000;t++){
        for(int i=0;i<nx;i++){
            s[i] += taus * (-s[i] + eta_emx[i] / (1. + w2 * f[i] + v2 * r[i]));
            r[i] += taur * (-r[i] + eta_pax[i] / (1. + v1 * s[i]));
            f[i] += tauf * (-f[i] + eta_fgf[i] / (1. + w1 * s[i]));
        }
    }

    // chemo-attraction gradient
    vector<double> pA(nx,0.), pB(nx,0.), pC(nx,0.);
    for (int i=0;i<nx;i++) {
        pA[i] = (kA/2.)*(1.+tanh((f[i]-theta1)/sigmaA));
        pB[i] = (kB/2.)*(1.+tanh((theta2-f[i])/sigmaB))*(kB/2.)*(1.+tanh((f[i]-theta3)/sigmaB));
        pC[i] = (kC/2.)*(1.+tanh((theta4-f[i])/sigmaC));
    }

    // Set up the chemo-interaction terms used to compute the axon branching flux, J
    vector<vector<double> > G(axontypes);
    for (int i=0; i<axontypes; ++i) {
        G[i].resize(nx,0.);
    }
    for (int i=0; i<axontypes; ++i) {
        vector<double> pAs = space (pA, gammaA[i]);
        vector<double> pBs = space (pB, gammaB[i]);
        vector<double> pCs = space (pC, gammaC[i]);
        for (int xi=0; xi<nx; xi++) {
            G[i][xi] = pAs[xi] + pBs[xi] + pCs[xi];
        }
    }
    // Initialise C to zero
    vector<vector<double> > C(axontypes);
    for (int i=0; i<axontypes; i++) {
        C[i].resize (nx,0.);
        for (int xi=0; xi<nx; xi++) {
            C[i][xi] = 0.;
        }
    }
    // Initialise A with random noise
    vector<vector<double> > A(axontypes);
    for (int i=0; i<axontypes; i++) {
        A[i].resize (nx,0.);
        for (int xi=0; xi<nx; xi++) {
            A[i][xi] = (randFloat()*0.1) + 0.8;
        }
    }

    // storage for results
    vector< vector <double*> > Q;
    {
        for (int i=0; i<axontypes; i++) {
            vector <double*> q;
            for (int xi=0; xi<nx; xi++) { q.push_back(&A[i][xi]); }
            Q.push_back (q);
        }
        for (int i=0; i<axontypes; i++){
            vector <double*> q;
            for (int xi=0; xi<nx; xi++) { q.push_back(&C[i][xi]); }
            Q.push_back (q);
        }
        {
            vector <double*> q;
            for (int xi=0; xi<nx; xi++) { q.push_back(&pA[xi]); }
            Q.push_back (q);
        }
        {
            vector <double*> q;
            for (int xi=0; xi<nx; xi++) { q.push_back(&pB[xi]); }
            Q.push_back (q);
        }
        {
            vector <double*> q;
            for (int xi=0; xi<nx;xi++) { q.push_back(&pC[xi]); }
            Q.push_back (q);
        }
        {
            vector <double*> q;
            for (int xi=0; xi<nx; xi++) { q.push_back(&f[xi]); }
            Q.push_back (q);
        }
    }


    // MAIN SIMULATION LOOP
    for (int t=0; t<stoi(argv[1]); t++) {

        // Compute Karb2004 Eq 3:
        vector<double> n(nx,0.);
        for (int xi=0; xi<nx; xi++) {
            for (int i=0; i<axontypes; i++) {
                n[xi] += C[i][xi];
            }
            n[xi] = 1. - n[xi];
        }

        // For each TC type (axontypes is 5, I think)
        for (int i=0; i<axontypes; i++) {
            vector<double> alpha_C_beta_nA(nx,0.);
            for (int xi=0; xi<nx; xi++) {
                alpha_C_beta_nA[xi] = alpha * C[i][xi] - beta * n[xi] * pow (A[i][xi], k);
            }
            // Runge-Kutta integration for A
            vector<double> q(nx,0.);
            vector<double> k1 = compute_axonalbranchflux (A[i], G[i], D, dx, alpha_C_beta_nA);
            for (int xi=0; xi<nx; xi++) {
                q[xi] = A[i][xi] + k1[xi] * halfdt;
            }

            vector<double> k2 = compute_axonalbranchflux (q, G[i], D, dx, alpha_C_beta_nA);
            for (int xi=0; xi<nx; xi++) {
                q[xi] = A[i][xi] + k2[xi] * halfdt;
            }

            vector<double> k3 = compute_axonalbranchflux (q, G[i], D, dx, alpha_C_beta_nA);
            for (int xi=0; xi<nx; xi++) {
                q[xi] = A[i][xi] + k3[xi] * dt;
            }

            vector<double> k4=compute_axonalbranchflux (q, G[i], D, dx, alpha_C_beta_nA);
            for (int xi=0; xi<nx; xi++) {
                A[i][xi] += (k1[xi] + 2. * (k2[xi] + k3[xi]) + k4[xi]) * sixthdt;
            }
        }

        for (int i=0; i<axontypes; i++) {
            vector<double> beta_n_A(nx,0.);
            for (int xi=0; xi<nx; xi++) {
                beta_n_A[xi] = beta * n[xi] * pow(A[i][xi],k);
            }
            // Runge-Kutta integration for C (or ci)
            vector<double> q(nx,0.);
            vector<double> k1 = compute_dci_dt (C[i], -alpha, beta_n_A);
            for (int xi=0; xi<nx; xi++) {
                q[xi] = C[i][xi] + k1[xi] * halfdt;
            }

            vector<double> k2 = compute_dci_dt (q, -alpha, beta_n_A);
            for (int xi=0; xi<nx; xi++) {
                q[xi] = C[i][xi] + k2[xi] * halfdt;
            }

            vector<double> k3 = compute_dci_dt (q, -alpha, beta_n_A);
            for (int xi=0; xi<nx; xi++) {
                q[xi] = C[i][xi] + k3[xi] * dt;
            }

            vector<double> k4 = compute_dci_dt (q, -alpha, beta_n_A);
            for (int xi=0; xi<nx; xi++) {
                C[i][xi] += (k1[xi]+2. * (k2[xi] + k3[xi]) + k4[xi]) * sixthdt;
            }
        }

    }

    // STORE DATA TO BINARY FILE
    {
        ofstream outFile;
        outFile.open(argv[2], ios::out|ios::trunc|ios::binary);
        for (unsigned int i=0; i<Q.size(); i++) {
            for (unsigned int j=0; j<Q[i].size(); j++) {
                outFile.write ((char*)Q[i][j], sizeof(double));
            }
        }
        outFile.close();
    }

    return 0;
}

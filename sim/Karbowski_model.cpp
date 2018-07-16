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

    for (int i=0; i<n; i++) {
        X.push_back (x[i]);
    }
    X.push_back(x[n-1]);

    for (int i=0; i<n; i++){
        x[i] = (X[i+2]-X[i]) * dx;
    }

    return x;
};

vector<double> funct1 (vector<double> x, double a, vector<double> b)
{
    vector<double> y;
    for (int i=0; i<x.size(); i++) {
        y.push_back (a*x[i]+b[i]);
    }
    return y;
};

vector<double> funct2 (vector<double> a, vector<double> g, double D, double dx, vector<double> extra)
{
    int n = a.size();
    double dJdX;
    vector<double> output(n,0.);
    for(int j=0;j<n;j++){
        if(j==0){ dJdX =
            D*(a[j+2]-a[j])     -a[j+1]*g[j+1];
        }
        else if(j==1) { dJdX =
            D*(a[j+2]-a[j])     -a[j+1]*g[j+1]
            -D*(a[j]-a[j-1])    +a[j-1]*g[j-1];
        }
        else if(j==n-1){ dJdX =
            -D*(a[j]-a[j-1])    +a[j]*g[j];
        }
        else if (j==n-2){ dJdX =
            D*(a[j+1]-a[j])     -a[j+1]*g[j+1]
            -D*(a[j]-a[j-2])    +a[j-1]*g[j-1];
        }
        else { dJdX =
            D*(a[j+2]-a[j])     -a[j+1]*g[j+1]
            -D*(a[j]-a[j-2])    +a[j-1]*g[j-1];
        }
        output[j] = dJdX*2.*dx+extra[j];
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

    int ng = 5;

    double Aemx = 1.34, Apax = 1.4, Afgf = 0.9;
    double Cemx = 25.6, Cpax = 27.3, Cfgf = 26.4;
    double w1 = 2.4, w2 = 2.1, v1 = 2.6, v2 = 2.7;

    vector<double> gammaA(ng,0.), gammaB(ng,0.), gammaC(ng,0.);
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
        eta_emx[i] = Aemx*exp(-((x-L)*(x-L)) / (Cemx*Cemx));
        eta_pax[i] = Apax*exp(-(x*x) / (Cpax*Cpax));
        eta_fgf[i] = Afgf*exp(-(x*x) / (Cfgf*Cfgf));
    }

    vector<double> s(nx,0.), r(nx,0.), f(nx,0.);
    double taus = 0.0001, taur = 0.0001, tauf = 0.0001;
    for(int t=0;t<300000;t++){
        for(int i=0;i<nx;i++){
            s[i]+=taus*(-s[i]+eta_emx[i]/(1.+w2*f[i]+v2*r[i]));
            r[i]+=taur*(-r[i]+eta_pax[i]/(1.+v1*s[i]));
            f[i]+=tauf*(-f[i]+eta_fgf[i]/(1.+w1*s[i]));
        }
    }

    // chemo-attraction gradient
    vector<double> pA(nx,0.), pB(nx,0.), pC(nx,0.);
    for (int i=0;i<nx;i++) {
        pA[i] = (kA/2.)*(1.+tanh((f[i]-theta1)/sigmaA));
        pB[i] = (kB/2.)*(1.+tanh((theta2-f[i])/sigmaB))*(kB/2.)*(1.+tanh((f[i]-theta3)/sigmaB));
        pC[i] = (kC/2.)*(1.+tanh((theta4-f[i])/sigmaC));
    }

    vector<vector<double> > G(ng);
    for (int i=0;i<ng;i++){
        G[i].resize(nx,0.);
    }
    for (int i=0;i<ng;i++) {
        vector<double> pAs = space(pA,gammaA[i]);
        vector<double> pBs = space(pB,gammaB[i]);
        vector<double> pCs = space(pC,gammaC[i]);
        for(int j=0;j<nx;j++){
            G[i][j] = pAs[j]+pBs[j]+pCs[j];
        }
    }
    vector<vector<double> > C(ng);
    for (int i=0; i<ng; i++) {
        C[i].resize (nx,0.);
        for (int j=0; j<nx; j++) {
            C[i][j] = 0.;
        }
    }
    vector<vector<double> > A(ng);
    for (int i=0; i<ng; i++) {
        A[i].resize (nx,0.);
        for (int j=0; j<nx; j++) {
            A[i][j] = (randFloat()*0.1) + 0.8;
        }
    }

    // storage for results
    vector< vector <double*> > Q;
    {
        for(int i=0;i<ng;i++){
            vector <double*> q;
            for (int j=0; j<nx; j++) { q.push_back(&A[i][j]); }
            Q.push_back (q);
        }
        for(int i=0;i<ng;i++){
            vector <double*> q;
            for (int j=0; j<nx; j++) { q.push_back(&C[i][j]); }
            Q.push_back (q);
        }
        {
            vector <double*> q;
            for (int j=0; j<nx; j++) { q.push_back(&pA[j]); }
            Q.push_back (q);
        }
        {
            vector <double*> q;
            for (int j=0; j<nx; j++) { q.push_back(&pB[j]); }
            Q.push_back (q);
        }
        {
            vector <double*> q;
            for (int j=0; j<nx;j++) { q.push_back(&pC[j]); }
            Q.push_back (q);
        }
        {
            vector <double*> q;
            for (int j=0; j<nx; j++) { q.push_back(&f[j]); }
            Q.push_back (q);
        }
    }


    // MAIN SIMULATION LOOP
    for (int t=0; t<stoi(argv[1]); t++) {

        vector<double> n(nx,0.);
        for (int j=0; j<nx; j++) {
            for (int i=0; i<ng; i++) {
                n[j] += C[i][j];
            }
            n[j] = 1. - n[j];
        }

        for (int i=0; i<ng; i++) {
            vector<double> extra(nx,0.);
            for (int j=0; j<nx; j++) {
                extra[j] = alpha * C[i][j] - beta * n[j] * pow (A[i][j], k);
            }
            // Runge-Kutta integration for A
            vector<double> q(nx,0.);
            vector<double> k1 = funct2 (A[i], G[i], D, dx, extra);
            for (int j=0; j<nx; j++) {
                q[j] = A[i][j] + k1[j] * halfdt;
            }

            vector<double> k2 = funct2 (q, G[i], D, dx, extra);
            for (int j=0; j<nx; j++) {
                q[j] = A[i][j] + k2[j] * halfdt;
            }

            vector<double> k3 = funct2 (q, G[i], D, dx, extra);
            for (int j=0; j<nx; j++) {
                q[j] = A[i][j] + k3[j] * dt;
            }

            vector<double> k4=funct2 (q, G[i], D, dx, extra);
            for (int j=0; j<nx; j++) {
                A[i][j] += (k1[j] + 2. * (k2[j] + k3[j]) + k4[j]) * sixthdt;
            }
        }

        for (int i=0; i<ng; i++) {
            vector<double> bd(nx,0.);
            for (int j=0; j<nx; j++) {
                bd[j] = beta * n[j] * pow(A[i][j],k);
            }
            // Runge-Kutta integration for C
            vector<double> q(nx,0.);
            vector<double> k1 = funct1 (C[i], -alpha, bd);
            for (int j=0; j<nx; j++) {
                q[j] = C[i][j] + k1[j] * halfdt;
            }

            vector<double> k2 = funct1 (q, -alpha, bd);
            for (int j=0; j<nx; j++) {
                q[j] = C[i][j] + k2[j] * halfdt;
            }

            vector<double> k3 = funct1 (q, -alpha, bd);
            for (int j=0; j<nx; j++) {
                q[j] = C[i][j] + k3[j] * dt;
            }

            vector<double> k4 = funct1 (q, -alpha, bd);
            for (int j=0; j<nx; j++) {
                C[i][j] += (k1[j]+2. * (k2[j] + k3[j]) + k4[j]) * sixthdt;
            }
        }

    }


    { // STORE DATA TO BINARY FILE
        ofstream outFile;
        outFile.open(argv[2], ios::out|ios::trunc|ios::binary);
        for (int i=0; i<Q.size(); i++) {
            for (int j=0; j<Q[i].size(); j++) {
                outFile.write ((char*)Q[i][j], sizeof(double));
            }
        }
        outFile.close();
    }

    return 0;
}

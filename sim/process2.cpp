#include "morph/world.h"
#include "morph/sockserve.h"
#include "morph/display.h"
#include "morph/tools.h"
#include "morph/ReadCurves.h"
#include "HexGrid.h"
#include <iostream>
//#include <ofstream>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>

using namespace std;

using morph::HexGrid;
using morph::ReadCurves;

/*!
 * Reaction diffusion system; 2-D Karbowski 2004.
 */
class RD_2D_Karb
{
public:
    /*!
     * Holds the number of hexes in the populated HexGrid
     */
    unsigned int nhex = 0;

    /*!
     * how many thalamo-cortical axon types are there? Denoted by N in
     * the paper, and so we use N here too.
     */
    unsigned int N = 5;

    /*!
     * These are the c_i(x,t) variables from the Karb2004 paper. x is
     * a vector in two-space.
     */
    vector<vector<double> > c;

    /*!
     * These are the a_i(x,t) variables from the Karb2004 paper. x is
     * a vector in two-space. The first vector is over the different
     * TC axon types, enumerated by i, the second vector are the a_i
     * values, indexed by the vi in the Hexes in HexGrid.
     */
    vector<vector<double> > a;

    /*!
     * For each TC axon type, this holds the two components of the
     * gradient field of the scalar value a(x,t) (where this x is a
     * vector in two-space)
     */
    vector<array<vector<double>, 2> > grad_a;

    /*!
     * Contains the chemo-attractant modifiers which are applied to
     * a_i(x,t) in Eq 4.
     */
    vector<array<vector<double>, 2> > chemo;

    /*!
     * n(x,t) variable from the Karb2004 paper.
     */
    vector<double> n;

    /*!
     * J_i(x,t) variables - the "flux current of axonal branches of
     * type i". This is a vector field.
     */
    vector<array<vector<double>, 2> > J;

    /*!
     * Holds the divergence of the J_i(x)s
     */
    vector<vector<double> > divJ;

    /*!
     * Our choice of dt.
     */
    double dt = 0.1;

    /*!
     * Compute half and sixth dt in constructor.
     */
    //@{
    double halfdt = 0.0;
    double sixthdt = 0.0;
    //@}

    /*!
     * The power to which a_i(x,t) is raised in Eqs 1 and 2 in the
     * paper.
     */
    double k = 3.0;

    /*!
     * The diffusion parameter.
     */
    double D = 0.1;

    /*!
     * alpha_i parameters
     */
    vector<double> alpha;

    /*!
     * beta_i parameters
     */
    vector<double> beta;

    /*!
     * gamma_A/B/C_i parameters from Eq 4
     */
    //@{
    vector<double> gammaA;
    vector<double> gammaB;
    vector<double> gammaC;
    //@}

    /*!
     * Variables for factor expression dynamics (Eqs 5-7)
     */
    //@{
    vector<double> eta_emx;
    vector<double> eta_pax;
    vector<double> eta_fgf;

    /*!
     * These are s(x), r(x) and f(x) in Karb2004.
     */
    //@}
    vector<double> emx;
    vector<double> pax;
    vector<double> fgf;
    //@}
    //@}

    /*!
     * Parameters for factor expression dynamics (Eqs 5-7)
     */
    //@{
    double Aemx = 1.34;
    double Apax = 1.4;
    double Afgf = 0.9;

    double Chiemx = 2.0;//0.094; //25.6;
    double Chipax = 2.2;//0.1;  //27.3;
    double Chifgf = 2.1;//0.098; //26.4

    double v1 = 2.6;
    double v2 = 2.7;
    double w1 = 2.4;
    double w2 = 2.1;

    /*!
     * Note: Using tau_emx, tau_pax, tau_fgf in place of tau_s, tau_r, tau_f
     */
    //@{
    double tau_emx = 0.0001;
    double tau_pax = 0.0001;
    double tau_fgf = 0.0001;
    //@}

    /*!
     * The directions of the change (in radians) in uncoupled factor
     * concentrations
     */
    //@{
    float diremx = 3.141593;
    float dirpax = 0.785; // norm 0
    float dirfgf = -0.785; // norm 0
    //@}

    double sigmaA = 0.2;
    double sigmaB = 0.2;
    double sigmaC = 0.2;

    double kA = 0.58;
    double kB = 0.9;
    double kC = 0.55;

    double theta1 = 0.77;
    double theta2 = 0.5;
    double theta3 = 0.39;
    double theta4 = 0.08;

    //@}

    /*!
     * Rho_A/B/C variables in Eq 4 - the concentrations of axon
     * guidance molecules A, B and C. In Karbowski 2004, these are
     * time independent and we will treat time as such, populating
     * them at initialisation.
     */
    vector<double> rhoA;
    vector<double> rhoB;
    vector<double> rhoC;

    /*!
     * Into grad_rhoA/B/C put the two components of the gradient of
     * rhoA/B/C computed across the HexGrid surface.
     */
    array<vector<double>, 2> grad_rhoA;
    array<vector<double>, 2> grad_rhoB;
    array<vector<double>, 2> grad_rhoC;

    /*!
     * The HexGrid "background" for the Reaction Diffusion system.
     */
    HexGrid* hg;

    /*!
     * Hex to hex distance. Populate this from hg.d after hg has been
     * initialised.
     */
    double dx = 1.0;

    /*!
     * Memory to hold an intermediate result
     */
    vector<vector<double> > betaterm;

    /*!
     * Holds an intermediate value for the computation of Eqs 1 and 2.
     */
    vector<vector<double> > alpha_c_beta_na;

    unsigned int stepCount = 0;

    /*!
     * Simple constructor; no arguments.
     */
    RD_2D_Karb (void) {
        this->halfdt = this->dt/2.0;
        this->sixthdt = this->dt/6.0;
    }

    /*!
     * Destructor required to free up HexGrid memory
     */
    ~RD_2D_Karb (void) {
        delete (this->hg);
    }

    /*!
     * A utility function to resize the vector-vectors that hold a
     * variable for the N different thalamo-cortical axon types.
     */
    void resize_vector_vector (vector<vector<double> >& vv) {
        vv.resize (this->N);
        for (unsigned int i =0; i<this->N; ++i) {
            vv[i].resize (this->nhex, 0.0);
        }
    }

    /*!
     * Resize a variable that'll be nhex elements long
     */
    void resize_vector_variable (vector<double>& v) {
        v.resize (this->nhex, 0.0);
    }

    /*!
     * Resize a parameter that'll be N elements long
     */
    void resize_vector_param (vector<double>& p) {
        p.resize (this->N, 0.0);
    }

    /*!
     * Resize a gradient field
     */
    void resize_gradient_field (array<vector<double>, 2>& g) {
        g[0].resize (this->nhex, 0.0);
        g[1].resize (this->nhex, 0.0);
    }

    /*!
     * Resize a vector (over TC types i) of an array of two
     * vector<double>s which are the x and y components of a
     * (mathematical) vector field.
     */
    void resize_vector_array_vector (vector<array<vector<double>, 2> >& vav) {
        vav.resize (this->N);
        for (unsigned int i = 0; i<this->N; ++i) {
            this->resize_gradient_field (vav[i]);
        }
    }

    /*!
     * Initialise this vector of vectors with noise. This is a
     * model-specific function.
     */
    void noiseify_vector_vector (vector<vector<double> >& vv) {
        for (unsigned int i = 0; i<this->N; ++i) {
            for (unsigned int h = 0; h < this->nhex; ++h) {
                // Note the model-specific choice of multiplier and offset here:
                vv[i][h] = morph::Tools::randDouble() * 0.1 + 0.8;
            }
        }
    }

    /*!
     * Initialise HexGrid, variables and parameters. Carry out
     * one-time computations of the model.
     */
    void init (vector<morph::Gdisplay>& displays) {

//        cout << "init() called" << endl;

        // Create a HexGrid
        this->hg = new HexGrid (0.01, 3);
        // Read the curves which make a boundary
        ReadCurves r("./trial.svg");
        // Set the boundary in the HexGrid
        this->hg->setBoundary (r.getCorticalPath());
        // Vector size comes from number of Hexes in the HexGrid
        this->nhex = this->hg->num();
        // Spatial dx comes from the HexGrid, too.
        this->dx = this->hg->getd();

        // Resize and zero-initialise the various containers
        this->resize_vector_vector (this->c);
        this->resize_vector_vector (this->a);
        this->resize_vector_vector (this->betaterm);
        this->resize_vector_vector (this->alpha_c_beta_na);
        this->resize_vector_vector (this->divJ);

        this->resize_vector_variable (this->n);
        this->resize_vector_variable (this->rhoA);
        this->resize_vector_variable (this->rhoB);
        this->resize_vector_variable (this->rhoC);

        this->resize_vector_variable (this->eta_emx);
        this->resize_vector_variable (this->eta_pax);
        this->resize_vector_variable (this->eta_fgf);

        this->resize_vector_variable (this->emx);
        this->resize_vector_variable (this->pax);
        this->resize_vector_variable (this->fgf);

        this->resize_vector_param (this->alpha);
        this->resize_vector_param (this->beta);
        this->resize_vector_param (this->gammaA);
        this->resize_vector_param (this->gammaB);
        this->resize_vector_param (this->gammaC);

        this->resize_gradient_field (this->grad_rhoA);
        this->resize_gradient_field (this->grad_rhoB);
        this->resize_gradient_field (this->grad_rhoC);

        // Resize grad_a and other vector-array-vectors
        this->resize_vector_array_vector (this->grad_a);
        this->resize_vector_array_vector (this->chemo);
        this->resize_vector_array_vector (this->J);

        // Initialise a with noise
        this->noiseify_vector_vector (this->a);

        // Populate parameters
        this->gammaA[0] =  1.6;
        this->gammaA[1] = -0.4;
        this->gammaA[2] = -2.21;
        this->gammaA[3] = -2.1;
        this->gammaA[4] = -2.45;

        this->gammaB[0] = -0.6;
        this->gammaB[1] = -0.5;
        this->gammaB[2] =  0.4;
        this->gammaB[3] = -0.5;
        this->gammaB[4] = -1.0;

        this->gammaC[0] = -2.9;
        this->gammaC[1] = -2.5;
        this->gammaC[2] = -2.23;
        this->gammaC[3] = -0.6;
        this->gammaC[4] =  1.7;

        this->alpha[0] = 1;
        this->alpha[1] = 1;
        this->alpha[2] = 1;
        this->alpha[3] = 1;
        this->alpha[4] = 1;

        this->beta[0] = 1;
        this->beta[1] = 1;
        this->beta[2] = 1;
        this->beta[3] = 1;
        this->beta[4] = 1;

        // Generate the assumed uncoupled concentrations of growth/transcription factors
        this->createFactorInitialConc (this->diremx, this->Aemx, this->Chiemx, this->eta_emx);
        this->createFactorInitialConc (this->dirpax, this->Apax, this->Chipax, this->eta_pax);
        this->createFactorInitialConc (this->dirfgf, this->Afgf, this->Chifgf, this->eta_fgf);

        // Run the expression dynamics, showing images as we go.
        this->runExpressionDynamics (displays);

        // Can now populate rhoA, rhoB and rhoC according to the paper.
        this->populateChemoAttractants (displays);

        // Compute gradients of guidance molecule concentrations once only
        this->spacegrad2D (this->rhoA, this->grad_rhoA);
        this->spacegrad2D (this->rhoB, this->grad_rhoB);
        this->spacegrad2D (this->rhoC, this->grad_rhoC);

        // Having computed gradients, build this->chemo; has to be done once only.
        for (unsigned int i=0; i<this->N; ++i) {
            for (unsigned int h=0; h<this->nhex; ++h) {
                // Here I'm adding, but does Ji have components?
                // Revisit this question. Yes, it does, but then we
                // compute the divergence of J to get
                // da_i/dt|diffusion
                this->chemo[i][0][h] = this->gammaA[i] * this->grad_rhoA[0][h]
                    + this->gammaB[i] * this->grad_rhoB[0][h]
                    + this->gammaC[i] * this->grad_rhoC[0][h];
                this->chemo[i][1][h] = this->gammaA[i] * this->grad_rhoA[1][h]
                    + this->gammaB[i] * this->grad_rhoB[1][h]
                    + this->gammaC[i] * this->grad_rhoC[1][h];
            }
        }
    }

    /*!
     * Compute the divergence of the vector field vf, placing the
     * result in div_vf.
     */
    void divergence (array<vector<double>, 2>& vf, vector<double>& div_vf) {
        // It's dvf_x/dx + dvf_y/dy. Our incoming vector field, vf is
        // already in plain old x,y coordinates, but we do need to
        // compute change of values between adjacent hexes.
        for (auto h : this->hg->hexen) {
            div_vf[h.vi] = 0.0;
            // Find x gradient of vf[0]
            if (h.has_ne && h.has_nw) {
                div_vf[h.vi] = (vf[0][h.ne->vi] - vf[0][h.nw->vi]) / ((double)h.d * 2.0);
            } else if (h.has_ne) {
                div_vf[h.vi] = (vf[0][h.ne->vi] - vf[0][h.vi]) / (double)h.d;
            } else if (h.has_nw) {
                div_vf[h.vi] = (vf[0][h.vi] - vf[0][h.nw->vi]) / (double)h.d;
            } else {
                // zero gradient in x direction as no neighbours in
                // those directions? Or possibly use the average of
                // the gradient between the nw,ne and sw,se neighbours
            }

            // Add on the y gradient of vf[1]
            if (h.has_nnw && h.has_nne && h.has_nsw && h.has_nse) {
                // Full complement. Compute the mean of the nse->nne and nsw->nnw gradients
                div_vf[h.vi] += ((vf[1][h.nne->vi] - vf[1][h.nse->vi]) + (vf[1][h.nnw->vi] - vf[1][h.nsw->vi])) / (double)h.getDv();

            } else if (h.has_nnw && h.has_nne ) {
                div_vf[h.vi] += ( (vf[1][h.nne->vi] + vf[1][h.nnw->vi]) / 2.0 - vf[1][h.vi]) / (double)h.getDv();

            } else if (h.has_nsw && h.has_nse) {
                div_vf[h.vi] += (vf[1][h.vi] - (vf[1][h.nse->vi] + vf[1][h.nsw->vi]) / 2.0) / (double)h.getDv();

            } else if (h.has_nnw && h.has_nsw) {
                div_vf[h.vi] += (vf[1][h.nnw->vi] - vf[1][h.nsw->vi]) / (double)h.getTwoDv();

            } else if (h.has_nne && h.has_nse) {
                div_vf[h.vi] += (vf[1][h.nne->vi] - vf[1][h.nse->vi]) / (double)h.getTwoDv();
            } else {
                // Leave grady at 0
            }
        }
    }

    /*!
     * 2D spatial integration of the function f. Result placed in gradf.
     *
     * For each Hex, work out the gradient in x and y directions
     * using whatever neighbours can contribute to an estimate.
     */
    void spacegrad2D (vector<double>& f, array<vector<double>, 2>& gradf) {


        // Note - East is positive x; North is positive y.
        for (auto h : this->hg->hexen) {

            gradf[0][h.vi] = 0.0;
            gradf[1][h.vi] = 0.0;

            // Find x gradient
            if (h.has_ne && h.has_nw) {
                gradf[0][h.vi] = (f[h.ne->vi] - f[h.nw->vi]) / ((double)h.d * 2.0);
            } else if (h.has_ne) {
                gradf[0][h.vi] = (f[h.ne->vi] - f[h.vi]) / (double)h.d;
            } else if (h.has_nw) {
                gradf[0][h.vi] = (f[h.vi] - f[h.nw->vi]) / (double)h.d;
            } else {
                // zero gradient in x direction as no neighbours in
                // those directions? Or possibly use the average of
                // the gradient between the nw,ne and sw,se neighbours
            }

            // Find y gradient
            if (h.has_nnw && h.has_nne && h.has_nsw && h.has_nse) {
                // Full complement. Compute the mean of the nse->nne and nsw->nnw gradients
                gradf[1][h.vi] = ((f[h.nne->vi] - f[h.nse->vi]) + (f[h.nnw->vi] - f[h.nsw->vi])) / (double)h.getDv();

            } else if (h.has_nnw && h.has_nne ) {
                gradf[1][h.vi] = ( (f[h.nne->vi] + f[h.nnw->vi]) / 2.0 - f[h.vi]) / (double)h.getDv();

            } else if (h.has_nsw && h.has_nse) {
                gradf[1][h.vi] = (f[h.vi] - (f[h.nse->vi] + f[h.nsw->vi]) / 2.0) / (double)h.getDv();

            } else if (h.has_nnw && h.has_nsw) {
                gradf[1][h.vi] = (f[h.nnw->vi] - f[h.nsw->vi]) / (double)h.getTwoDv();

            } else if (h.has_nne && h.has_nse) {
                gradf[1][h.vi] = (f[h.nne->vi] - f[h.nse->vi]) / (double)h.getTwoDv();
            } else {
                // Leave grady at 0
            }
        }
    }

#if 0
    vector<double> getLaplacian (vector<double> Q, double dx) {
        double overdxSquare = 1.0/(dx*dx);
        vector<double> L(nhex, 0.0);
        for (int i=0; i<nhex; i++) {
            L[i] = (Q[N[i][0]] + Q[N[i][1]] + Q[N[i][2]] + Q[N[i][3]] + Q[N[i][4]] + Q[N[i][5]] - 6.0 * Q[i]) * overdxSquare;
        }
        return L;
    }
#endif

    /*!
     * Do a step through the model.
     *
     * T460s manages about 20 of these per second.
     */
    void step (void) {

        this->stepCount++;

//        if (this->stepCount % 100 == 0) {
//            cout << "100 steps done..." << endl;
//        }

        // 1. Compute Karb2004 Eq 3. (coupling between connections made by each TC type)
        for (auto h : this->hg->hexen) {
            for (unsigned int i=0; i<N; ++i) {
                n[h.vi] += c[i][h.vi];
            }
            n[h.vi] = 1. - n[h.vi];
        }

        // 2. Do integration of a (RK in the 1D model). Involves computing axon branching flux.

        // Pre-compute intermediate val:
        // FIXME Could precompute betaterm first, then use it in this avoiding more pow(a, k)s
        for (unsigned int i=0; i<this->N; ++i) {
            for (unsigned int h=0; h<this->nhex; ++h) {
                this->alpha_c_beta_na[i][h] = alpha[i] * c[i][h] - beta[i] * n[h] * pow (a[i][h], k);
            }
        }

        // Runge-Kutta:
        for (unsigned int i=0; i<this->N; ++i) {

            // Runge-Kutta integration for A
            vector<double> q(this->nhex, 0.0);
            this->compute_axonalbranchflux (a[i], i); // populates divJ[i]
            vector<double> k1(this->nhex, 0.0);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k1[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                q[h] = this->a[i][h] + k1[h] * halfdt;
            }

            vector<double> k2(this->nhex, 0.0);
            this->compute_axonalbranchflux (q, i);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k2[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                q[h] = this->a[i][h] + k2[h] * halfdt;
            }

            vector<double> k3(this->nhex, 0.0);
            this->compute_axonalbranchflux (q, i);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k3[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                q[h] = this->a[i][h] + k3[h] * dt;
            }

            vector<double> k4(this->nhex, 0.0);
            this->compute_axonalbranchflux (q, i);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k4[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                a[i][h] += (k1[h] + 2.0 * (k2[h] + k3[h]) + k4[h]) * sixthdt;
//                if (abs(a[i][h]) > 0) {
//                    cout << "a[i][h]>0 : " << a[i][h] << endl;
//                }
            }
        }

        // 3. Do integration of c
        for (unsigned int i=0; i<this->N; ++i) {

            for (unsigned int h=0; h<nhex; h++) {
                this->betaterm[i][h] = beta[i] * n[h] * pow (a[i][h], k);
            }
            // Runge-Kutta integration for C (or ci)
            vector<double> q(nhex,0.);
            vector<double> k1 = compute_dci_dt (c[i], i);
            for (unsigned int h=0; h<nhex; h++) {
                q[h] = c[i][h] + k1[h] * halfdt;
            }

            vector<double> k2 = compute_dci_dt (q, i);
            for (unsigned int h=0; h<nhex; h++) {
                q[h] = c[i][h] + k2[h] * halfdt;
            }

            vector<double> k3 = compute_dci_dt (q, i);
            for (unsigned int h=0; h<nhex; h++) {
                q[h] = c[i][h] + k3[h] * dt;
            }

            vector<double> k4 = compute_dci_dt (q, i);
            for (unsigned int h=0; h<nhex; h++) {
                c[i][h] += (k1[h]+2. * (k2[h] + k3[h]) + k4[h]) * sixthdt;
            }
        }
    }

    /*!
     * Plot the system on @a disps
     */
    void plot (vector<morph::Gdisplay>& disps) {
        this->plot_f (this->a, disps, 6);
        this->plot_f (this->c, disps, 11);
    }

    /*!
     * Plot a or c
     */
    void plot_f (vector<vector<double> >& f, vector<morph::Gdisplay>& disps, unsigned int display_offset) {

        vector<double> fix(3, 0.0);
        vector<double> eye(3, 0.0);
        eye[2] = -0.4;
        vector<double> rot(3, 0.0);

        // Copies data to plot out of the model
        vector<double> maxa (5, -1e7);
        vector<double> mina (5, +1e7);
        // Determines min and max
        for (auto h : this->hg->hexen) {
            if (h.onBoundary() == false) {
                for (unsigned int i = 0; i<this->N; ++i) {
                    if (f[i][h.vi]>maxa[i]) { maxa[i] = f[i][h.vi]; }
                    if (f[i][h.vi]<mina[i]) { mina[i] = f[i][h.vi]; }
                }
            }
        }
        vector<double> scalea (5, 0);
        for (unsigned int i = 0; i<this->N; ++i) {
            scalea[i] = 1.0 / (maxa[i]-mina[i]);
        }

        // Determine a colour from min, max and current value
        vector<vector<double> > norm_a;
        this->resize_vector_vector (norm_a);
        for (unsigned int i = 0; i<this->N; ++i) {
            for (unsigned int h=0; h<this->nhex; h++) {
                norm_a[i][h] = fmin (fmax (((f[i][h]) - mina[i]) * scalea[i], 0.0), 1.0);
            }
        }

        // Draw
        for (unsigned int i = 0; i<this->N; ++i) {
            disps[i+display_offset].resetDisplay (fix, eye, rot);
            for (auto h : this->hg->hexen) {
                array<float,3> cl_a = morph::Tools::getJetColorF (norm_a[i][h.vi]);
                disps[i+display_offset].drawHex (h.position(), (h.d/2.0f), cl_a);
            }
            disps[i+display_offset].redrawDisplay();
        }
    }

    /*!
     * Plot expression of emx, pax and fgf
     */
    void plotexpression (vector<morph::Gdisplay>& disps) {

        vector<double> fix(3, 0.0);
        vector<double> eye(3, 0.0);
        eye[2] = -0.4;
        vector<double> rot(3, 0.0);

        // Copies data to plot out of the model
        double maxemx = -1e7;
        double minemx = +1e7;
        double maxpax = -1e7;
        double minpax = +1e7;
        double maxfgf = -1e7;
        double minfgf = +1e7;
        // Determines min and max
        for (auto h : this->hg->hexen) {
            if (h.onBoundary() == false) {
                if (this->emx[h.vi]>maxemx) { maxemx = this->emx[h.vi]; }
                if (this->emx[h.vi]<minemx) { minemx = this->emx[h.vi]; }

                if (this->pax[h.vi]>maxpax) { maxpax = this->pax[h.vi]; }
                if (this->pax[h.vi]<minpax) { minpax = this->pax[h.vi]; }

                if (this->fgf[h.vi]>maxfgf) { maxfgf = this->fgf[h.vi]; }
                if (this->fgf[h.vi]<minfgf) { minfgf = this->fgf[h.vi]; }
            }
        }
        double scaleemx = 1.0 / (maxemx-minemx);
        double scalepax = 1.0 / (maxpax-minpax);
        double scalefgf = 1.0 / (maxfgf-minfgf);

        // Determine a colour from min, max and current value
        vector<double> norm_emx(this->nhex, 0.0);
        vector<double> norm_pax(this->nhex, 0.0);
        vector<double> norm_fgf(this->nhex, 0.0);
        for (unsigned int h=0; h<this->nhex; h++) {
            norm_emx[h] = fmin (fmax (((this->emx[h]) - minemx) * scaleemx, 0.0), 1.0);
            norm_pax[h] = fmin (fmax (((this->pax[h]) - minpax) * scalepax, 0.0), 1.0);
            norm_fgf[h] = fmin (fmax (((this->fgf[h]) - minfgf) * scalefgf, 0.0), 1.0);
        }

        // Step through vectors or iterate through list? The latter should be just fine here.
        disps[0].resetDisplay (fix, eye, rot);
        for (auto h : this->hg->hexen) {
            array<float,3> cl_emx = morph::Tools::getJetColorF (norm_emx[h.vi]);
            disps[0].drawHex (h.position(), (h.d/2.0f), cl_emx);
        }
        disps[0].redrawDisplay();

        disps[1].resetDisplay (fix, eye, rot);
        for (auto h : this->hg->hexen) {
            array<float,3> cl_pax = morph::Tools::getJetColorF (norm_pax[h.vi]);
            disps[1].drawHex (h.position(), (h.d/2.0f), cl_pax);
        }
        disps[1].redrawDisplay();

        disps[2].resetDisplay (fix, eye, rot);
        for (auto h : this->hg->hexen) {
            array<float,3> cl_fgf = morph::Tools::getJetColorF (norm_fgf[h.vi]);
            disps[2].drawHex (h.position(), (h.d/2.0f), cl_fgf);
        }
        disps[2].redrawDisplay();
    }

    /*!
     * Plot concentrations of chemo-attractor molecules A, B and C.
     */
    void plotchemo (vector<morph::Gdisplay>& disps) {

        vector<double> fix(3, 0.0);
        vector<double> eye(3, 0.0);
        eye[2] = -0.4;
        vector<double> rot(3, 0.0);

        // Copies data to plot out of the model
        double maxrhoA = -1e7;
        double minrhoA = +1e7;
        double maxrhoB = -1e7;
        double minrhoB = +1e7;
        double maxrhoC = -1e7;
        double minrhoC = +1e7;
        // Determines min and max
        for (auto h : this->hg->hexen) {
            if (h.onBoundary() == false) {
                if (this->rhoA[h.vi]>maxrhoA) { maxrhoA = this->rhoA[h.vi]; }
                if (this->rhoA[h.vi]<minrhoA) { minrhoA = this->rhoA[h.vi]; }

                if (this->rhoB[h.vi]>maxrhoB) { maxrhoB = this->rhoB[h.vi]; }
                if (this->rhoB[h.vi]<minrhoB) { minrhoB = this->rhoB[h.vi]; }

                if (this->rhoC[h.vi]>maxrhoC) { maxrhoC = this->rhoC[h.vi]; }
                if (this->rhoC[h.vi]<minrhoC) { minrhoC = this->rhoC[h.vi]; }
            }
        }
        double scalerhoA = 1.0 / (maxrhoA-minrhoA);
        double scalerhoB = 1.0 / (maxrhoB-minrhoB);
        double scalerhoC = 1.0 / (maxrhoC-minrhoC);

        // Determine a colour from min, max and current value
        vector<double> norm_rhoA(this->nhex, 0.0);
        vector<double> norm_rhoB(this->nhex, 0.0);
        vector<double> norm_rhoC(this->nhex, 0.0);
        for (unsigned int h=0; h<this->nhex; h++) {
            norm_rhoA[h] = fmin (fmax (((this->rhoA[h]) - minrhoA) * scalerhoA, 0.0), 1.0);
            norm_rhoB[h] = fmin (fmax (((this->rhoB[h]) - minrhoB) * scalerhoB, 0.0), 1.0);
            norm_rhoC[h] = fmin (fmax (((this->rhoC[h]) - minrhoC) * scalerhoC, 0.0), 1.0);
        }

        // Step through vectors or iterate through list? The latter should be just fine here.
        disps[3].resetDisplay (fix, eye, rot);
        for (auto h : this->hg->hexen) {
            array<float,3> cl_rhoA = morph::Tools::getJetColorF (norm_rhoA[h.vi]);
            disps[3].drawHex (h.position(), (h.d/2.0f), cl_rhoA);
        }
        disps[3].redrawDisplay();

        disps[4].resetDisplay (fix, eye, rot);
        for (auto h : this->hg->hexen) {
            array<float,3> cl_rhoB = morph::Tools::getJetColorF (norm_rhoB[h.vi]);
            disps[4].drawHex (h.position(), (h.d/2.0f), cl_rhoB);
        }
        disps[4].redrawDisplay();

        disps[5].resetDisplay (fix, eye, rot);
        for (auto h : this->hg->hexen) {
            array<float,3> cl_rhoC = morph::Tools::getJetColorF (norm_rhoC[h.vi]);
            disps[5].drawHex (h.position(), (h.d/2.0f), cl_rhoC);
        }
        disps[5].redrawDisplay();
    }

    /*!
     * Does: f = (alpha * f) + betaterm. c.f. Karb2004, Eq 1. f is
     * c[i] or q from the RK algorithm.
     */
    vector<double> compute_dci_dt (vector<double>& f, unsigned int i) {
        vector<double> dci_dt (this->nhex, 0.0);
        for (unsigned int h=0; h<this->nhex; h++) {
            dci_dt[h] = this->betaterm[i][h] - this->alpha[i] * f[h];
        }
        return dci_dt;
    }

    /*!
     * Computes the "flux of axonal branches" term, J_i(x) (Eq 4)
     *
     * Inputs: this->chemo, f (which is this->a[i] or a q in the RK
     * algorithm), this->D, @a i, the TC type.  Helper functions:
     * spacegrad2D(), divergence().  Output: this->divJ
     */
    void compute_axonalbranchflux (vector<double>& f, unsigned int i) {

        // Compute gradient of a_i(x)
        this->spacegrad2D (f, this->grad_a[i]);
        // Compute J
        for (unsigned int h = 0; h<this->nhex; ++h) {
            this->J[i][0][h] = this->D * this->grad_a[i][0][h] - f[h] * this->chemo[i][0][h];
            this->J[i][1][h] = this->D * this->grad_a[i][1][h] - f[h] * this->chemo[i][1][h];
        }
        // Compute divergence of J
        this->divergence (this->J[i], this->divJ[i]);
    }

    /*!
     * Create a 2-D scalar field which follows a curve along one
     * direction (at angle @a phi radians, anti-clockwise from East),
     * being constant in the orthogonal direction. Place result into
     * @a result.
     *
     * @param Afac 'A' parameter for factor fac. c.f. Aemx, Apax, etc
     * in Karb2004.
     *
     * @param chifac 'chi' parameter for factor fac. c.f. Chi_emx, Chi_pax, etc
     */
    void createFactorInitialConc (float phi, double Afac, double chifac, vector<double>& result) {

        // Work in a co-ordinate system rotated by phi radians, called x_, y_
        double x_ = 0.0;

        double cosphi = (double) cos (phi);
        double sinphi = (double) sin (phi);
//        cout << "cosphi: " << cosphi << endl;
        // Get minimum x and maximum x in the rotated co-ordinate system.
        double x_min_ = this->hg->getXmin (phi);
//        cout << "x_min_: " << x_min_ << endl;

        for (auto h : this->hg->hexen) {
            // Rotate x, then offset by the minimum along that line
            x_ = (h.x * cosphi) + (h.y * sinphi) - x_min_;
            // x here is x from the Hex.
            result[h.vi] = Afac * exp (-(x_ * x_) / (chifac * chifac));
        }
    }

    /*!
     * Execute Eqs 5-7 of the Karbowski paper to find the steady state
     * of the growth/transcription factors after they have interacted
     * for a long time.
     */
    void runExpressionDynamics (vector<morph::Gdisplay>& displays) {
        for (unsigned int t=0; t<30000; ++t) { // 300000 matches Stuart's 1D Karbowski model
            for (auto h : this->hg->hexen) {
                emx[h.vi] += tau_emx * (-emx[h.vi] + eta_emx[h.vi] / (1. + w2 * fgf[h.vi] + v2 * pax[h.vi]));
                pax[h.vi] += tau_pax * (-pax[h.vi] + eta_pax[h.vi] / (1. + v1 * emx[h.vi]));
                fgf[h.vi] += tau_fgf * (-fgf[h.vi] + eta_fgf[h.vi] / (1. + w1 * emx[h.vi]));
            }
//            if (t%1000 == 0) {
//                cout << "Plot for t=" << t << endl;
///                this->plotexpression (displays);
//            }
        }
    }

    /*!
     * Using this->emx, this->pax and this->fgf, populate rhoA/B/C
     */
    void populateChemoAttractants (vector<morph::Gdisplay>& displays) {
        // chemo-attraction gradient. cf Fig 1 of Karb 2004
        for (unsigned int h=0; h<this->nhex; ++h) {
            this->rhoA[h] = (kA/2.)*(1.+tanh((fgf[h]-theta1)/sigmaA));
            this->rhoB[h] = (kB/2.)*(1.+tanh((theta2-fgf[h])/sigmaB))*(kB/2.)*(1.+tanh((fgf[h]-theta3)/sigmaB));
            this->rhoC[h] = (kC/2.)*(1.+tanh((theta4-fgf[h])/sigmaC));

        }
///        this->plotchemo (displays);
    }

    /*!
     * Create a symmetric, 2D Gaussian hill centred at coordinate (x,y) with
     * width sigma and height gain. Place result into @a result.
     */
    void createGaussian (float x, float y, double gain, double sigma, vector<double>& result) {

        // Once-only parts of the calculation of the Gaussian.
        double root_2_pi = 2.506628275;
        double one_over_sigma_root_2_pi = 1 / sigma * root_2_pi;
        double two_sigma_sq = 2 * sigma * sigma;

        // Gaussian dist. result, and a running sum of the results:
        double gauss = 0.0;
        double sum = 0.0;

        // x and y components of the vector from (x,y) to any given Hex.
        float rx = 0.0f, ry = 0.0f;
        // distance from any Hex to (x,y)
        float r = 0.0f;

        // Calculate each element of the kernel:
        for (auto h : this->hg->hexen) {
            rx = x - h.x;
            ry = y - h.y;
            r = sqrt (rx*rx + ry*ry);
            gauss = gain * (one_over_sigma_root_2_pi
                            * exp ( static_cast<double>(-(r*r))
                                    / two_sigma_sq ));
            result[h.vi] = gauss;
            sum += gauss;
            ++k;
        }

        // Normalise the kernel to 1 by dividing by the sum:
        unsigned int j = this->nhex;
        while (j > 0) {
            --j;
            result[j] = result[j] / sum;
        }
    }

}; // RD_2D_Karb


int main (int argc, char **argv)
{
    srand(atoi(argv[3]));

    // DISPLAYS
    vector<morph::Gdisplay> displays;
    vector<double> fix(3, 0.0);
    vector<double> eye(3, 0.0);
    eye[2] = -0.4;
    vector<double> rot(3, 0.0);
#if 0
    displays.push_back (morph::Gdisplay (600, "emx", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "pax", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "fgf", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "rhoA", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "rhoB", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "rhoC", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "a[0]", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "a[1]", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "a[2]", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "a[3]", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "a[4]", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "c[0]", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "c[1]", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "c[2]", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "c[3]", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    displays.push_back (morph::Gdisplay (600, "c[4]", 0.0, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();
#endif
    // Instantiate the model object
    RD_2D_Karb M;
    try {
        M.init (displays);
    } catch (const exception& e) {
        cerr << "Exception initialising RD_2D_Karb object: " << e.what() << endl;
    }

    morph::World W(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), M.dt);

    // Keep track of the frame number
    unsigned int frameN = 0;

    // Keep track of the time
    double TIME=0.;
    vector <double*> f;
    f.push_back(&TIME);

    // Start the loop
    bool doing = true;
    while (doing) {

        std::stringstream TIMEss;
        TIMEss << setw(10) << setfill('0') << TIME;
        const char* TIMEcs = TIMEss.str().c_str();

        std::stringstream out;
        out.clear();
        out.setf (ios::fixed, ios::floatfield);

        // Start making an output message
        out << TIME << ",";

        // Check for a command from the model world
        vector<string> command;
        string messageI = W.master.exchange (out.str().c_str());
        stringstream ss(messageI);
        while (ss.good()) {
            string substr;
            getline (ss,substr, ',');
            command.push_back (substr);
        }
        ss.clear();

        // Interpret commands:
        switch (stoi(command[0])) {

        case 0: // *** QUIT ***
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 0=QUIT" << endl;

            W.logfile.close();
            for (unsigned int i=0; i<displays.size(); i++) {
                displays[i].closeDisplay();
            }
            W.master.closeSocket();
            for (unsigned int i=0; i<W.ports.size(); i++) {
                W.ports[i].closeSocket();
            }
            doing = false;
            break;
        }

        case 1: // *** STEP ***
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 1=STEP" << endl;

            // Exchange comms [should be vector<vector<string> > commands = W.getcommand()]
            vector<vector<string> > commands(W.ports.size());
            for (unsigned int i=0; i<W.ports.size(); i++) {
                string messageI = W.ports[i].exchange (out.str().c_str());
                stringstream ss(messageI);
                while (ss.good()) {
                    string substr;
                    getline (ss, substr, ',');
                    commands[i].push_back (substr);
                }
                ss.clear();
            }

            // Step the model
            try {
                M.step();
            } catch (const exception& e) {
                W.logfile << "Caught exception calling M.step(): " << e.what() << endl;
            }

            TIME++;
            break;
        }

        case 2: // *** PLOT ***
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 2=PLOT" << endl;
            displays[0].resetDisplay (fix, eye, rot);
            try {
                M.plot (displays);
            } catch (const exception& e) {
                W.logfile << "Caught exception calling M.plot(): " << e.what() << endl;
            }
            break;
        }

        case 3:
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 3=RECD" << endl;
            std::stringstream frameFile1;
            frameFile1 << "logs/" << W.processName << "frameA";
            frameFile1 << setw(5) << setfill('0') << frameN;
            frameFile1 << ".png";
            displays[0].saveImage (frameFile1.str());
            frameN++;
            break;
        }

        case 4:
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 4=RECD" << endl;

            // Set up parameters from sim.py like this
            switch(stoi(command[1])){
            case 0:
            {
                //Dn[stoi(command[2])] = stod(command[3]);
                break;
            }
            case 1:
            {
                //Dc[stoi(command[2])] = stod(command[3]);
                break;
            }
            }
            break;
        }

        case 5: // *** SAVE ***
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 5=SAVE" << endl;

            if (command.size() == 2) {

                std::stringstream oFile; oFile << command[1];
                //M.save (command[1].str());
#if 0
                ofstream outFile;
                std::stringstream oFile; oFile << command[1];
                outFile.open (oFile.str().c_str(), ios::out|ios::binary);

                double n = (double)S.size();
                double* N = &n;
                outFile.write ((char*)N, sizeof(double));
                for (unsigned int i=0; i<S.size(); i++) {
                    double m = (double)S[i].size();
                    double* M = &m;
                    outFile.write ((char*)M, sizeof(double));
                    for (unsigned int j=0; j<S[i].size(); j++) {
                        outFile.write ((char*)S[i][j], sizeof(double));
                    }
                }
                outFile.close();
#endif
            } else {
                W.logfile << "No output filename." << endl;
            }
            break;
        }

        case 6: // *** LOAD ***
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 6=LOAD" << endl;
#if 0
            if (command.size() == 2) {
                double dummy;
                ifstream inFile;
                std::stringstream iFile;
                iFile << command[1];
                inFile.open (iFile.str().c_str(), ios::in|ios::binary);
                inFile.read ((char*)&dummy, sizeof(double));
                int I = (int)dummy;
                if (I == static_cast<int>(S.size())) {
                    for (int i=0; i<I; i++) {
                        inFile.read ((char*)&dummy, sizeof(double));
                        int J = (int)dummy;
                        if (J == static_cast<int>(S[i].size())) {
                            for (int j=0; j<J; j++) {
                                inFile.read ((char*)S[i][j], sizeof(double));
                            }
                        } else {
                            W.logfile << "Wrong dims I." << endl;
                        }
                    }
                } else {
                    W.logfile << "Wrong dims J." << endl;
                }
                inFile.close();
            } else {
                W.logfile << "No input filename." << endl;
            }
#endif
            break;
        }

        }
    }

    return 0;
};

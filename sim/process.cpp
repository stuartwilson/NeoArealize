#include "morph/world.h"
#include "morph/sockserve.h"
#include "morph/display.h"
#include "morph/tools.h"
#include "morph/ReadCurves.h"
#include "morph/HexGrid.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>

//#include <omp.h>

#define DEBUG 1
#define DBGSTREAM std::cout
#include <morph/MorphDbg.h>

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
     * Constants
     */
    //@{
    //! Square root of 3 over 2
    const double R3_OVER_2 = 0.866025403784439;
    //! Square root of 3
    const double ROOT3 = 1.73205080756888;
    //@}

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
    vector<array<vector<double>, 2> > g;

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
    double dt = 0.0001;

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
    double Chipax = 2.5;//0.1;  //27.3;
    double Chifgf = 2.25;//0.098; //26.4

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
    float dirpax = 1.2; // norm 0
    float dirfgf = 0; // norm 0
    //@}

    double sigmaA = 0.1;
    double sigmaB = 0.1;
    double sigmaC = 0.2;

    double kA = 0.1;
    double kB = 6;
    double kC = 0.9;

    double theta1 = 0.1; // 0.77 orig.
    double theta2 = 0.4; // 0.5
    double theta3 = 0.6; // 0.39
    double theta4 = 0.1; // 0.08

    //@} end factor expression dynamics parameters

    /*!
     * Rho_A/B/C variables in Eq 4 - the concentrations of axon
     * guidance molecules A, B and C. In Karbowski 2004, these are
     * time independent and we will treat time as such, populating
     * them at initialisation.
     */
    //@{
    vector<double> rhoA;
    vector<double> rhoB;
    vector<double> rhoC;
    //@}

    /*!
     * Into grad_rhoA/B/C put the two components of the gradient of
     * rhoA/B/C computed across the HexGrid surface.
     */
    //@{
    array<vector<double>, 2> grad_rhoA;
    array<vector<double>, 2> grad_rhoB;
    array<vector<double>, 2> grad_rhoC;
    //@}

    /*!
     * The HexGrid "background" for the Reaction Diffusion system.
     */
    HexGrid* hg;

    /*!
     * Hex to hex distance. Populate this from hg.d after hg has been
     * initialised.
     */
    double d = 1.0;

    /*!
     * Memory to hold an intermediate result
     */
    vector<vector<double> > betaterm;

    /*!
     * Holds an intermediate value for the computation of Eqs 1 and 2.
     */
    vector<vector<double> > alpha_c_beta_na;

    /*!
     * Track the number of computational steps that we've carried
     * out. Only to show a message saying "100 steps done...", but
     * that's reason enough.
     */
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
#define ZERO_TOWARDS_EDGE 1
#ifdef ZERO_TOWARDS_EDGE
            for (auto h : this->hg->hexen) {
                // boundarySigmoid. Jumps sharply (100, larger is
                // sharper) over length scale 0.05 to 1. So if
                // distance from boundary > 0.05, noise has normal
                // value. Close to boundary, noise is less.
                vv[i][h.vi] = morph::Tools::randDouble() * 0.1;// + 0.8;
                if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                    double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary-0.01)) ); // 0.05!
                    vv[i][h.vi] = vv[i][h.vi] * bSig;
                }
            }
#else
            for (unsigned int h = 0; h < this->nhex; ++h) {
                // Note the model-specific choice of multiplier and offset here:
                vv[i][h] = morph::Tools::randDouble() * 0.1 + 0.8;
            }
#endif
        }
    }

    /*!
     * Initialise HexGrid, variables and parameters. Carry out
     * one-time computations of the model.
     */
    void init (vector<morph::Gdisplay>& displays) {

        DBG ("called");

        // Create a HexGrid
        this->hg = new HexGrid (0.01, 3);
        // Read the curves which make a boundary
        ReadCurves r("./trial.svg");
        // Set the boundary in the HexGrid
        this->hg->setBoundary (r.getCorticalPath());
        // Compute the distances from the boundary
        this->hg->computeDistanceToBoundary();
        // Vector size comes from number of Hexes in the HexGrid
        this->nhex = this->hg->num();
        // Spatial d comes from the HexGrid, too.
        this->d = this->hg->getd();

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
        this->resize_vector_array_vector (this->g);
        this->resize_vector_array_vector (this->J);

        // Initialise a with noise
        this->noiseify_vector_vector (this->a);

        // Populate parameters
        double gammagain = 20.0;
        this->gammaA[0] =  1.6 * gammagain;
        this->gammaA[1] = -0.4 * gammagain;
        this->gammaA[2] = -2.21 * gammagain;
        this->gammaA[3] = -2.1 * gammagain;
        this->gammaA[4] = -2.45 * gammagain;

        this->gammaB[0] = -0.6 * gammagain;
        this->gammaB[1] = -0.5 * gammagain;
        this->gammaB[2] =  0.4 * gammagain;
        this->gammaB[3] = -0.5 * gammagain;
        this->gammaB[4] = -1.0 * gammagain;

        this->gammaC[0] = -2.9 * gammagain;
        this->gammaC[1] = -2.5 * gammagain;
        this->gammaC[2] = -2.23 * gammagain;
        this->gammaC[3] = -0.6 * gammagain;
        this->gammaC[4] =  1.7 * gammagain;

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

        // Having computed gradients, build this->g; has to be done once only.
        for (unsigned int i=0; i<this->N; ++i) {
#ifdef ZERO_TOWARDS_EDGE
            for (auto h : this->hg->hexen) {
                // Sigmoid/logistic fn params: 100 sharpness, 0.01 dist offset from boudnary
                double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary-0.02)) );
                this->g[i][0][h.vi] = (this->gammaA[i] * this->grad_rhoA[0][h.vi]
                                    + this->gammaB[i] * this->grad_rhoB[0][h.vi]
                                    + this->gammaC[i] * this->grad_rhoC[0][h.vi]) * bSig;
                this->g[i][1][h.vi] = (this->gammaA[i] * this->grad_rhoA[1][h.vi]
                                    + this->gammaB[i] * this->grad_rhoB[1][h.vi]
                                    + this->gammaC[i] * this->grad_rhoC[1][h.vi]) * bSig;
            }
#else
            for (unsigned int h=0; h<this->nhex; ++h) {
                // this->g is what I called g(x) in my methods writeup.
                this->g[i][0][h] = this->gammaA[i] * this->grad_rhoA[0][h]
                    + this->gammaB[i] * this->grad_rhoB[0][h]
                    + this->gammaC[i] * this->grad_rhoC[0][h];
                this->g[i][1][h] = this->gammaA[i] * this->grad_rhoA[1][h]
                    + this->gammaB[i] * this->grad_rhoB[1][h]
                    + this->gammaC[i] * this->grad_rhoC[1][h];
            }
#endif
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
            if (h.d == 0.0) {
                throw runtime_error ("h.d is unexpectedly 0.");
            }
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
                div_vf[h.vi] += ((vf[1][h.nne->vi] - vf[1][h.nse->vi]) + (vf[1][h.nnw->vi] - vf[1][h.nsw->vi])) / (double)h.getV();

            } else if (h.has_nnw && h.has_nne ) {
                div_vf[h.vi] += ( (vf[1][h.nne->vi] + vf[1][h.nnw->vi]) / 2.0 - vf[1][h.vi]) / (double)h.getV();

            } else if (h.has_nsw && h.has_nse) {
                div_vf[h.vi] += (vf[1][h.vi] - (vf[1][h.nse->vi] + vf[1][h.nsw->vi]) / 2.0) / (double)h.getV();

            } else if (h.has_nnw && h.has_nsw) {
                div_vf[h.vi] += (vf[1][h.nnw->vi] - vf[1][h.nsw->vi]) / (double)h.getTwoV();

            } else if (h.has_nne && h.has_nse) {
                div_vf[h.vi] += (vf[1][h.nne->vi] - vf[1][h.nse->vi]) / (double)h.getTwoV();
            } else {
                // Leave grady at 0
            }
        }
    }

    /*!
     * Using the boundary integral, compute the divergence of the
     * vector field vf, placing the result in div_vf.
     */
    void divergence_boundary (array<vector<double>, 2>& vf, vector<double>& div_vf) {
        // It's dvf_x/dx + dvf_y/dy. Our incoming vector field, vf is
        // already in plain old x,y coordinates, but we do need to
        // compute change of values between adjacent hexes.
        for (auto h : this->hg->hexen) {
            div_vf[h.vi] = 0.0;
            if (h.d == 0.0) {
                throw runtime_error ("h.d is unexpectedly 0.");
            }

            // First sum
            if (h.has_ne) {
                div_vf[h.vi] += /*cos (0)*/ (vf[0][h.ne->vi] + vf[0][h.vi]);
            } else {
                div_vf[h.vi] += (2.0 * vf[0][h.vi]);
            }
            if (h.has_nne) {
                div_vf[h.vi] += /*cos (60)*/ 0.5 * (vf[0][h.nne->vi] + vf[0][h.vi]);
            } else {
                div_vf[h.vi] += (/*0.5 * 2.0 * */ vf[0][h.vi]);
            }
            if (h.has_nnw) {
                div_vf[h.vi] -= /*cos (120)*/ 0.5 * (vf[0][h.nnw->vi] + vf[0][h.vi]);
            } else {
                div_vf[h.vi] -= (/*0.5 * 2.0 * */ vf[0][h.vi]);
            }
            if (h.has_nw) {
                div_vf[h.vi] -= /*cos (180)*/ (vf[0][h.nw->vi] + vf[0][h.vi]);
            } else {
                div_vf[h.vi] -= (2.0 * vf[0][h.vi]);
            }
            if (h.has_nsw) {
                div_vf[h.vi] -= /*cos (240)*/ 0.5 * (vf[0][h.nsw->vi] + vf[0][h.vi]);
            } else {
                div_vf[h.vi] -= (/*0.5 * 2.0 * */ vf[0][h.vi]);
            }
            if (h.has_nse) {
                div_vf[h.vi] += /*cos (300)*/ 0.5 * (vf[0][h.nse->vi] + vf[0][h.vi]);
            } else {
                div_vf[h.vi] += (/*0.5 * 2.0 * */ vf[0][h.vi]);
            }

            // 2nd sum
            //div_vf[h.vi] += sin (0) * (vf[1][h.ne->vi] + vf[1][h.vi]);
            if (h.has_nne) {
                div_vf[h.vi] += /*sin (60)*/ R3_OVER_2 * (vf[1][h.nne->vi] + vf[1][h.vi]);
            } else {
                div_vf[h.vi] += ROOT3 * vf[1][h.vi];
            }
            if (h.has_nnw) {
                div_vf[h.vi] += /*sin (120)*/ R3_OVER_2 * (vf[1][h.nnw->vi] + vf[1][h.vi]);
            } else {
                div_vf[h.vi] += ROOT3 * vf[1][h.vi];
            }
            //div_vf[h.vi] += sin (180) * (vf[1][h.nw->vi] + vf[1][h.vi]);
            if (h.has_nsw) {
                div_vf[h.vi] -= /*sin (240)*/ R3_OVER_2 * (vf[1][h.nsw->vi] + vf[1][h.vi]);
            } else {
                div_vf[h.vi] -= ROOT3 * vf[1][h.vi];
            }
            if (h.has_nse) {
                div_vf[h.vi] -= /*sin (300)*/ R3_OVER_2 * (vf[1][h.nse->vi] + vf[1][h.vi]);
            } else {
                div_vf[h.vi] -= ROOT3 * vf[1][h.vi];
            }
            div_vf[h.vi] /= 2.0;

        }
    }

    /*!
     * Examine the value in each Hex of the hexgrid of the scalar
     * field f. If abs(f[h]) exceeds the size of dangerThresh, then
     * output debugging information.
     */
    void debug_values (vector<double>& f, double dangerThresh) {
        for (auto h : this->hg->hexen) {
            if (abs(f[h.vi]) > dangerThresh) {
                DBG ("Blow-up threshold exceeded at Hex.vi=" << h.vi << " ("<< h.ri <<","<< h.gi <<")" <<  ": " << f[h.vi]);
                unsigned int wait = 0;
                while (wait++ < 120) {
                    usleep (1000000);
                }
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

            DBG2 ("(h.ri,h.gi): (" << h.ri << "," << h.gi << ")");
            // Find x gradient
            if (h.has_ne && h.has_nw) {
                DBG2 ("x case 1 f[h.ne]: " << f[h.ne->vi] << " - f[h.nw]" << f[h.nw->vi] << "/ h.d*2: " << (double)h.d * 2.0);
                gradf[0][h.vi] = (f[h.ne->vi] - f[h.nw->vi]) / ((double)h.d * 2.0);
            } else if (h.has_ne) {
                DBG2 ("x case 2 f[h.ne]: " << f[h.ne->vi] << " - f[h]" << f[h.vi] << "/ h.d: " << (double)h.d);
                gradf[0][h.vi] = (f[h.ne->vi] - f[h.vi]) / (double)h.d;
            } else if (h.has_nw) {
                DBG2 ("x case 3 f[h]: " << f[h.vi] << " - f[h.nw]" << f[h.nw->vi] << "/ h.d: " << (double)h.d);
                gradf[0][h.vi] = (f[h.vi] - f[h.nw->vi]) / (double)h.d;
            } else {
                // zero gradient in x direction as no neighbours in
                // those directions? Or possibly use the average of
                // the gradient between the nw,ne and sw,se neighbours
            }

            // Find y gradient
            if (h.has_nnw && h.has_nne && h.has_nsw && h.has_nse) {
                // Full complement. Compute the mean of the nse->nne and nsw->nnw gradients
#ifdef DEBUG2
                if (h.vi == 0) {
                    double _d = (double)h.getV();
                    double _td = (double)h.getTwoV();
                    DBG2 ("y case 1. getV: " << _d << " getTwoV: " << _td);
                }
#endif
                gradf[1][h.vi] = ((f[h.nne->vi] - f[h.nse->vi]) + (f[h.nnw->vi] - f[h.nsw->vi])) / (double)h.getV();

            } else if (h.has_nnw && h.has_nne ) {
                //if (h.vi == 0) { DBG ("y case 2"); }
                gradf[1][h.vi] = ( (f[h.nne->vi] + f[h.nnw->vi]) / 2.0 - f[h.vi]) / (double)h.getV();

            } else if (h.has_nsw && h.has_nse) {
                //if (h.vi == 0) { DBG ("y case 3"); }
                gradf[1][h.vi] = (f[h.vi] - (f[h.nse->vi] + f[h.nsw->vi]) / 2.0) / (double)h.getV();

            } else if (h.has_nnw && h.has_nsw) {
                //if (h.vi == 0) { DBG ("y case 4"); }
                gradf[1][h.vi] = (f[h.nnw->vi] - f[h.nsw->vi]) / (double)h.getTwoV();

            } else if (h.has_nne && h.has_nse) {
                //if (h.vi == 0) { DBG ("y case 5"); }
                gradf[1][h.vi] = (f[h.nne->vi] - f[h.nse->vi]) / (double)h.getTwoV();
            } else {
                // Leave grady at 0
            }

            //if (h.vi == 0) {
            //    DBG ("gradf[0/1][0]: " << gradf[0][0] << "," << gradf[1][0]);
            //}
        }
    }

    /*!
     * Do a step through the model.
     *
     * T460s manages about 20 of these per second.
     */
    void step (void) {

        this->stepCount++;

        if (this->stepCount % 100 == 0) {
            DBG ("System computed " << this->stepCount << " times so far...");
        }

        // 1. Compute Karb2004 Eq 3. (coupling between connections made by each TC type)
        double nsum = 0.0;
        double csum = 0.0;
        for (auto h : this->hg->hexen) {
            n[h.vi] = 0; // whoops forgot this!
            #pragma omp parallel for
            for (unsigned int i=0; i<N; ++i) {
                n[h.vi] += c[i][h.vi];
            }
            csum += c[0][h.vi];
            n[h.vi] = 1. - n[h.vi];
            nsum += n[h.vi];
        }
        if (this->stepCount % 20 == 0) {
            DBG("sum of all n is " << nsum);
            DBG("sum of all c for i=0 is " << csum);
        }

        // 2. Do integration of a (RK in the 1D model). Involves computing axon branching flux.

        // Pre-compute intermediate val:
        // FIXME Could precompute betaterm first, then use it in this avoiding more pow(a, k)s
        for (unsigned int i=0; i<this->N; ++i) {
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                this->alpha_c_beta_na[i][h] = alpha[i] * c[i][h] - beta[i] * n[h] * pow (a[i][h], k);
            }
        }

        // Runge-Kutta:
        #pragma omp parallel for
        for (unsigned int i=0; i<this->N; ++i) {

            DBG2 ("alpha_c_beta_na["<<i<<"][0] = " << this->alpha_c_beta_na[i][0]);

            // Runge-Kutta integration for A
            vector<double> q(this->nhex, 0.0);
            this->compute_divJ (a[i], i); // populates divJ[i]
            DBG2 ("Computing divJ, divJ[0][0]: " << divJ[0][0]);
            vector<double> k1(this->nhex, 0.0);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k1[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                q[h] = this->a[i][h] + k1[h] * halfdt;
            }
            DBG2 ("After RK stage 1, q[0]: " << q[0]);

            vector<double> k2(this->nhex, 0.0);
            this->compute_divJ (q, i);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k2[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                q[h] = this->a[i][h] + k2[h] * halfdt; // Kaboom!
            }
            DBG2 ("After RK stage 2, q[0]:" << q[0] << " from a["<<i<<"][0]:" << a[i][0] << " divj["<<i<<"][0]:" << divJ[i][0] << " k2[0]:" << k2[0]);

            vector<double> k3(this->nhex, 0.0);
            this->compute_divJ (q, i);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k3[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                q[h] = this->a[i][h] + k3[h] * dt;
            }
            DBG2 ("After RK stage 3, q[0]: " << q[0]);

            vector<double> k4(this->nhex, 0.0);
            this->compute_divJ (q, i);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k4[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                a[i][h] += (k1[h] + 2.0 * (k2[h] + k3[h]) + k4[h]) * sixthdt;
            }
            DBG2 ("After RK stage 4, a[" << i << "][0]: " << a[i][0]);

            DBG2("Debug a["<<i<<"]");
            //this->debug_values (a[i], 1e8);
        }

        // 3. Do integration of c
        #pragma omp parallel for
        for (unsigned int i=0; i<this->N; ++i) {

            for (unsigned int h=0; h<nhex; h++) {
                this->betaterm[i][h] = beta[i] * n[h] * pow (a[i][h], k);
            }
            DBG2 ("betaterm[" << i << "][0]: " << betaterm[i][0]);

            // Runge-Kutta integration for C (or ci)
            vector<double> q(nhex,0.);
            vector<double> k1 = compute_dci_dt (c[i], i);
            for (unsigned int h=0; h<nhex; h++) {
                q[h] = c[i][h] + k1[h] * halfdt;
            }
            DBG2 ("After RK stage 1, q[0]: " << q[0]);

            vector<double> k2 = compute_dci_dt (q, i);
            for (unsigned int h=0; h<nhex; h++) {
                q[h] = c[i][h] + k2[h] * halfdt;
            }
            DBG2 ("After RK stage 2, q[0]: " << q[0]);

            vector<double> k3 = compute_dci_dt (q, i);
            for (unsigned int h=0; h<nhex; h++) {
                q[h] = c[i][h] + k3[h] * dt;
            }
            DBG2 ("After RK stage 3, q[0]: " << q[0]);

            vector<double> k4 = compute_dci_dt (q, i);
            for (unsigned int h=0; h<nhex; h++) {
                c[i][h] += (k1[h]+2. * (k2[h] + k3[h]) + k4[h]) * sixthdt;
            }
            DBG2 ("After RK stage 4, c["<<i<<"][0]: " << c[i][0]);

            DBG2("Debug c["<<i<<"]");
            //this->debug_values (c[i], 1e8);
        }
    }

    /*!
     * Plot the system on @a disps
     */
    void plot (vector<morph::Gdisplay>& disps) {
        this->plot_f (this->a, disps[2]);
        this->plot_f (this->c, disps[3]);

        // To enable examination, keep replotting these:
        this->plotchemo (disps);
        this->plotexpression (disps);
    }

    /*!
     * Plot a or c
     */
    void plot_f (vector<vector<double> >& f, morph::Gdisplay& disp) {
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

        // Create an offset which we'll increment by the width of the
        // map, starting from the left-most map (f[0])
        float hgwidth = this->hg->getXmax()-this->hg->getXmin();
        array<float,3> offset = { 2*(-hgwidth-(hgwidth/20)), 0.0f, 0.0f };

        // Draw
        disp.resetDisplay (fix, eye, rot);
        for (unsigned int i = 0; i<this->N; ++i) {
            for (auto h : this->hg->hexen) {
#ifdef PLOT_WITH_Z
                disp.drawHex (h.position(), offset, (h.d/2.0f), norm_a[i][h.vi]);
#else
                array<float,3> cl_a = morph::Tools::getJetColorF (norm_a[i][h.vi]);
                disp.drawHex (h.position(), offset, (h.d/2.0f), cl_a);
#endif
            }
            offset[0] += hgwidth + (hgwidth/20);
        }
        disp.redrawDisplay();
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

        // Set offsets for the three maps that we'll plot
        float hgwidth = this->hg->getXmax()-this->hg->getXmin();
        array<float,3> offset1 = { -hgwidth-(hgwidth/20), 0.0f, 0.0f };
        array<float,3> offset2 = { 0.0f, 0.0f, 0.0f };
        array<float,3> offset3 = { hgwidth+(hgwidth/20), 0.0f, 0.0f };

#ifdef PLOT_EXPRESSION_WITH_Z
        for (auto h : this->hg->hexen) {
            disps[0].drawHex (h.position(), offset1, (h.d/2.0f), norm_emx[h.vi]);
        }
        for (auto h : this->hg->hexen) {
            disps[0].drawHex (h.position(), offset2, (h.d/2.0f), norm_pax[h.vi]);
        }
        for (auto h : this->hg->hexen) {
            disps[0].drawHex (h.position(), offset3, (h.d/2.0f), norm_fgf[h.vi]);
        }
#else
        for (auto h : this->hg->hexen) {
            array<float,3> cl_emx = morph::Tools::getJetColorF (norm_emx[h.vi]);
            disps[0].drawHex (h.position(), offset1, (h.d/2.0f), cl_emx);
        }
        for (auto h : this->hg->hexen) {
            array<float,3> cl_pax = morph::Tools::getJetColorF (norm_pax[h.vi]);
            disps[0].drawHex (h.position(), offset2, (h.d/2.0f), cl_pax);
        }
        for (auto h : this->hg->hexen) {
            array<float,3> cl_fgf = morph::Tools::getJetColorF (norm_fgf[h.vi]);
            disps[0].drawHex (h.position(), offset3, (h.d/2.0f), cl_fgf);
        }
#endif
        disps[0].redrawDisplay();
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

        // Set offsets for the three maps that we'll plot
        float hgwidth = this->hg->getXmax()-this->hg->getXmin();
        array<float,3> offset1 = { -hgwidth-(hgwidth/20), 0.0f, 0.0f };
        array<float,3> offset2 = { 0.0f, 0.0f, 0.0f };
        array<float,3> offset3 = { hgwidth+(hgwidth/20), 0.0f, 0.0f };

        // Step through vectors or iterate through list? The latter should be just fine here.
        disps[1].resetDisplay (fix, eye, rot);
#ifdef PLOT_CHEMO_WITH_Z
        for (auto h : this->hg->hexen) {
            disps[1].drawHex (h.position(), offset1, (h.d/2.0f), norm_rhoA[h.vi]);
        }

        for (auto h : this->hg->hexen) {
            disps[1].drawHex (h.position(), offset2, (h.d/2.0f), norm_rhoB[h.vi]);
        }

        for (auto h : this->hg->hexen) {
            disps[1].drawHex (h.position(), offset3, (h.d/2.0f), norm_rhoC[h.vi]);
        }
#else
        for (auto h : this->hg->hexen) {
            array<float,3> cl_rhoA = morph::Tools::getJetColorF (norm_rhoA[h.vi]);
            disps[1].drawHex (h.position(), offset1, (h.d/2.0f), cl_rhoA);
        }

        for (auto h : this->hg->hexen) {
            array<float,3> cl_rhoB = morph::Tools::getJetColorF (norm_rhoB[h.vi]);
            disps[1].drawHex (h.position(), offset2, (h.d/2.0f), cl_rhoB);
        }

        for (auto h : this->hg->hexen) {
            array<float,3> cl_rhoC = morph::Tools::getJetColorF (norm_rhoC[h.vi]);
            disps[1].drawHex (h.position(), offset3, (h.d/2.0f), cl_rhoC);
        }
#endif
        disps[1].redrawDisplay();
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
     * Inputs: this->g, fa (which is this->a[i] or a q in the RK
     * algorithm), this->D, @a i, the TC type.  Helper functions:
     * spacegrad2D(), divergence().  Output: this->divJ
     */
    void compute_divJ (vector<double>& fa, unsigned int i) {

// Stable with dt = 0.0001;
#define VECTOR_CALCULUS_EXPANSION_METHOD 1
//#define VECTOR_CALCULUS_EXPANSION_METHOD_BOUNDARY_FORCED 1
// Stable(ish) with dt = 0.001;
//#define NAIVE_METHOD_WITH_BOUNDARY_TESTING 1
// Unstable
//#define NAIVE_METHOD 1

#ifdef VECTOR_CALCULUS_EXPANSION_METHOD
        // Three terms to compute; see Eq. 14 in methods_notes.pdf

        // Compute gradient of a_i(x), for use computing the third term, below.
        this->spacegrad2D (fa, this->grad_a[i]);

        for (auto h : this->hg->hexen) {
            // 1. The D Del^2 a_i term
            // Compute the sum around the neighbours
            double thesum = -6 * fa[h.vi];
            if (h.has_ne) {
                thesum += fa[h.ne->vi];
            } else {
                // Apply boundary condition
            }
            if (h.has_nne) {
                thesum += fa[h.nne->vi];
            } else {
                thesum += fa[h.vi]; // A ghost neighbour-east with same value as Hex_0
            }
            if (h.has_nnw) {
                thesum += fa[h.nnw->vi];
            } else {
                thesum += fa[h.vi];
            }
            if (h.has_nw) {
                thesum += fa[h.nw->vi];
            } else {
                thesum += fa[h.vi];
            }
            if (h.has_nsw) {
                thesum += fa[h.nsw->vi];
            } else {
                thesum += fa[h.vi];
            }
            if (h.has_nse) {
                thesum += fa[h.nse->vi];
            } else {
                thesum += fa[h.vi];
            }
            // Multiply bu 2D/3d^2
            double term1 = (this->D * 2) / (3 * this->d * this->d) * thesum;

            // 2. The a div(g) term. Two sums for this.
            double term2 = 0.0;
            // First sum
            if (h.has_ne) {
                term2 += /*cos (0)*/ (this->g[i][0][h.ne->vi] + this->g[i][0][h.vi]);
            } else {
                // Boundary condition _should_ be satisfied by
                // sigmoidal roll-off of g towards the boundary, so
                // add only g[i][0][h.vi]
                term2 += /*cos (0)*/ (this->g[i][0][h.vi]);
            }
            if (h.has_nne) {
                term2 += /*cos (60)*/ 0.5 * (this->g[i][0][h.nne->vi] + this->g[i][0][h.vi]);
            } else {
                term2 += /*cos (60)*/ 0.5 * (this->g[i][0][h.vi]);
            }
            if (h.has_nnw) {
                term2 -= /*cos (120)*/ 0.5 * (this->g[i][0][h.nnw->vi] + this->g[i][0][h.vi]);
            } else {
                term2 -= /*cos (120)*/ 0.5 * (this->g[i][0][h.vi]);
            }
            if (h.has_nw) {
                term2 -= /*cos (180)*/ (this->g[i][0][h.nw->vi] + this->g[i][0][h.vi]);
            } else {
                term2 -= /*cos (180)*/ (this->g[i][0][h.vi]);
            }
            if (h.has_nsw) {
                term2 -= /*cos (240)*/ 0.5 * (this->g[i][0][h.nsw->vi] + this->g[i][0][h.vi]);
            } else {
                term2 -= /*cos (240)*/ 0.5 * (this->g[i][0][h.vi]);
            }
            if (h.has_nse) {
                term2 += /*cos (300)*/ 0.5 * (this->g[i][0][h.nse->vi] + this->g[i][0][h.vi]);
            } else {
                term2 += /*cos (300)*/ 0.5 * (this->g[i][0][h.vi]);
            }
            // 2nd sum
            //term2 += sin (0) * (this->g[i][1][h.ne->vi] + this->g[i][1][h.vi]);
            if (h.has_nne) {
                term2 += /*sin (60)*/ R3_OVER_2 * (this->g[i][1][h.nne->vi] + this->g[i][1][h.vi]);
            } else {
                term2 += /*sin (60)*/ R3_OVER_2 * (this->g[i][1][h.vi]);
            }
            if (h.has_nnw) {
                term2 += /*sin (120)*/ R3_OVER_2 * (this->g[i][1][h.nnw->vi] + this->g[i][1][h.vi]);
            } else {
                term2 += /*sin (120)*/ R3_OVER_2 * (this->g[i][1][h.vi]);
            }
            //term2 += sin (180) * (this->g[i][1][h.nw->vi] + this->g[i][1][h.vi]);
            if (h.has_nsw) {
                term2 -= /*sin (240)*/ R3_OVER_2 * (this->g[i][1][h.nsw->vi] + this->g[i][1][h.vi]);
            } else {
                term2 -= /*sin (240)*/ R3_OVER_2 * (this->g[i][1][h.vi]);
            }
            if (h.has_nse) {
                term2 -= /*sin (300)*/ R3_OVER_2 * (this->g[i][1][h.nse->vi] + this->g[i][1][h.vi]);
            } else {
                term2 -= /*sin (300)*/ R3_OVER_2 * (this->g[i][1][h.vi]);
            }

            term2 /= (3.0 * this->d);
            term2 *= fa[h.vi];

            // 3. Third term is this->g . grad a_i. Should not
            // contribute to J, as g(x) decays towards boundary.
            double term3 = this->g[i][0][h.vi] * this->grad_a[i][0][h.vi]
                + this->g[i][1][h.vi] * this->grad_a[i][1][h.vi];

            this->divJ[i][h.vi] = term1 + term2 + term3;
        }
#endif

#ifdef VECTOR_CALCULUS_EXPANSION_METHOD_BOUNDARY_FORCED
        // Three terms to compute; see Eq. 14 in methods_notes.pdf

        // Compute gradient of a_i(x), for use computing the third term, below.
        this->spacegrad2D (fa, this->grad_a[i]);

        for (auto h : this->hg->hexen) {
            if (h.onBoundary() == true) {
                // Force divJ to 0 on boundary
                this->divJ[i][h.vi] = 0;

            } else {
                // 1. The D Del^2 a_i term
                // Compute the sum around the neighbours
                double thesum = -6 * fa[h.vi];
                thesum += fa[h.ne->vi];
                thesum += fa[h.nne->vi];
                thesum += fa[h.nnw->vi];
                thesum += fa[h.nw->vi];
                thesum += fa[h.nsw->vi];
                thesum += fa[h.nse->vi];
                // Multiply bu 2D/3d^2
                double term1 = (this->D * 2) / (3 * this->d * this->d) * thesum;

                // 2. The a div(g) term. Two sums for this.
                double term2 = 0.0;
                // First sum
                term2 += /*cos (0)*/ (this->g[i][0][h.ne->vi] + this->g[i][0][h.vi]);
                term2 += /*cos (60)*/ 0.5 * (this->g[i][0][h.nne->vi] + this->g[i][0][h.vi]);
                term2 -= /*cos (120)*/ 0.5 * (this->g[i][0][h.nnw->vi] + this->g[i][0][h.vi]);
                term2 -= /*cos (180)*/ (this->g[i][0][h.nw->vi] + this->g[i][0][h.vi]);
                term2 -= /*cos (240)*/ 0.5 * (this->g[i][0][h.nsw->vi] + this->g[i][0][h.vi]);
                term2 += /*cos (300)*/ 0.5 * (this->g[i][0][h.nse->vi] + this->g[i][0][h.vi]);
                // 2nd sum
                //term2 += sin (0) * (this->g[i][1][h.ne->vi] + this->g[i][1][h.vi]);
                term2 += /*sin (60)*/ R3_OVER_2 * (this->g[i][1][h.nne->vi] + this->g[i][1][h.vi]);
                term2 += /*sin (120)*/ R3_OVER_2 * (this->g[i][1][h.nnw->vi] + this->g[i][1][h.vi]);
                //term2 += sin (180) * (this->g[i][1][h.nw->vi] + this->g[i][1][h.vi]);
                term2 -= /*sin (240)*/ R3_OVER_2 * (this->g[i][1][h.nsw->vi] + this->g[i][1][h.vi]);
                term2 -= /*sin (300)*/ R3_OVER_2 * (this->g[i][1][h.nse->vi] + this->g[i][1][h.vi]);
                term2 /= 2.0;
                term2 *= fa[h.vi];

                // 3. Third term is this->g . grad a_i
                double term3 = this->g[i][0][h.vi] * this->grad_a[i][0][h.vi]
                    + this->g[i][1][h.vi] * this->grad_a[i][1][h.vi];

                this->divJ[i][h.vi] = term1 + term2 + term3;
            }
        }
#endif

#ifdef NAIVE_METHOD_WITH_BOUNDARY_TESTING
        // Compute gradient of a_i(x)
        this->spacegrad2D (fa, this->grad_a[i]);
        // Compute J. J blows up if grad_a blows up.
        for (auto h : this->hg->hexen) {
            if (h.onBoundary() == true) {
                // Force J to 0 on boundary. Is there a less brutal way?
                this->J[i][0][h.vi] = 0;
                this->J[i][1][h.vi] = 0;
            } else {
                this->J[i][0][h.vi] = this->D * this->grad_a[i][0][h.vi] - fa[h.vi] * this->g[i][0][h.vi];
                this->J[i][1][h.vi] = this->D * this->grad_a[i][1][h.vi] - fa[h.vi] * this->g[i][1][h.vi];
            }
        }
        // Compute divergence of J
        this->divergence_boundary (this->J[i], this->divJ[i]);
#endif

#ifdef NAIVE_METHOD
        // Compute gradient of a_i(x)
        this->spacegrad2D (fa, this->grad_a[i]);
        // Compute J. J blows up if grad_a blows up.
        for (unsigned int h = 0; h<this->nhex; ++h) {
            this->J[i][0][h] = this->D * this->grad_a[i][0][h] - fa[h] * this->g[i][0][h];
            this->J[i][1][h] = this->D * this->grad_a[i][1][h] - fa[h] * this->g[i][1][h];
        }
        // Compute divergence of J
        this->divergence_boundary (this->J[i], this->divJ[i]);
#endif
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
        DBG2 ("cosphi: " << cosphi);
        // Get minimum x and maximum x in the rotated co-ordinate system.
        double x_min_ = this->hg->getXmin (phi);
        DBG2 ("x_min_: " << x_min_);

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
        #pragma omp parallel for
        for (unsigned int t=0; t<100000; ++t) { // 300000 matches Stuart's 1D Karbowski model
            for (auto h : this->hg->hexen) {
                emx[h.vi] += tau_emx * (-emx[h.vi] + eta_emx[h.vi] / (1. + w2 * fgf[h.vi] + v2 * pax[h.vi]));
                pax[h.vi] += tau_pax * (-pax[h.vi] + eta_pax[h.vi] / (1. + v1 * emx[h.vi]));
                fgf[h.vi] += tau_fgf * (-fgf[h.vi] + eta_fgf[h.vi] / (1. + w1 * emx[h.vi]));
            }
            // Incompatible with parallel:
            //if (t%1000 == 0) {
            //    DBG ("Plot for t=" << t);
            //    this->plotexpression (displays);
            //}
        }
        this->plotexpression (displays);
    }

    /*!
     * Using this->emx, this->pax and this->fgf, populate rhoA/B/C
     */
    void populateChemoAttractants (vector<morph::Gdisplay>& displays) {
        // chemo-attraction gradient. cf Fig 1 of Karb 2004
        #pragma omp parallel for
        for (unsigned int h=0; h<this->nhex; ++h) {
            this->rhoA[h] = (kA/2.)*(1.+tanh((fgf[h]-theta1)/sigmaA));
            this->rhoB[h] = (kB/2.)*(1.+tanh((theta2-fgf[h])/sigmaB))*(kB/2.)*(1.+tanh((fgf[h]-theta3)/sigmaB));
            this->rhoC[h] = (kC/2.)*(1.+tanh((theta4-fgf[h])/sigmaC));

        }
        this->plotchemo (displays);
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

// Define this to run without the python server. Run as:
// process world00 logs/log00.txt 1
#define NO_NETWORKING 1

int main (int argc, char **argv)
{
    srand(atoi(argv[3]));

    // DISPLAYS
    vector<morph::Gdisplay> displays;
    vector<double> fix(3, 0.0);
    vector<double> eye(3, 0.0);
    eye[2] = -0.4;
    vector<double> rot(3, 0.0);

    double rhoInit = 1.5;
    string worldName(argv[1]);
    string winTitle = worldName + ": emx_pax_fgf";
    displays.push_back (morph::Gdisplay (1020, 300, winTitle.c_str(), rhoInit, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    winTitle = worldName + ": rhoA_rhoB_rhoC";
    displays.push_back (morph::Gdisplay (1020, 300, winTitle.c_str(), rhoInit, 0.0, 0.0, displays[0].win));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    winTitle = worldName + ": a[0] to a[4]";
    displays.push_back (morph::Gdisplay (1700, 300, winTitle.c_str(), rhoInit, 0.0, 0.0, displays[0].win));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    winTitle = worldName + ": c[0] to c[4]";
    displays.push_back (morph::Gdisplay (1700, 300, winTitle.c_str(), rhoInit, 0.0, 0.0, displays[0].win));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    // Instantiate the model object
    RD_2D_Karb M;
    try {
        M.init (displays);
    } catch (const exception& e) {
        cerr << "Exception initialising RD_2D_Karb object: " << e.what() << endl;
    }
#ifdef NO_NETWORKING
    morph::World W(argv[1], argv[2], atoi(argv[3]), M.dt);
#else
    morph::World W(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), M.dt);
#endif

    // Keep track of the frame number
    unsigned int frameN = 0;

    // Keep track of the time
    double TIME=0.;
    vector <double*> f;
    f.push_back(&TIME);

    // Start the loop
    bool doing = true;
    bool do_a_step = true;
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
#ifndef NO_NETWORKING
        string messageI = W.master.exchange (out.str().c_str());
        stringstream ss(messageI);
        while (ss.good()) {
            string substr;
            getline (ss,substr, ',');
            command.push_back (substr);
        }
        ss.clear();
#else
        command.push_back ("1");
#endif

#ifdef NO_NETWORKING
        // Force step/plot/step/plot etc.
        if (do_a_step == true) {
            command[0] = "1";
            do_a_step = false;
        } else {
            command[0] = "2";
            do_a_step = true;
        }
#endif
        // Interpret commands:
        switch (stoi(command[0])) {

        case 0: // *** QUIT ***
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 0=QUIT" << endl;

            W.logfile.close();
            for (unsigned int i=0; i<displays.size(); i++) {
                displays[i].closeDisplay();
            }
#ifndef NO_NETWORKING
            W.master.closeSocket();
            for (unsigned int i=0; i<W.ports.size(); i++) {
                W.ports[i].closeSocket();
            }
#endif
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
                //M.save (command[1].str());
            } else {
                W.logfile << "No output filename." << endl;
            }
            break;
        }

        case 6: // *** LOAD ***
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 6=LOAD" << endl;
            if (command.size() == 2) {
                //M.load(command[1]);
            } else {
                W.logfile << "No input filename." << endl;
            }
            break;
        }

        }
    }

    return 0;
};

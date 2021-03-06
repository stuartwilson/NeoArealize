#include "morph/display.h"
#include "morph/tools.h"
#include "morph/ReadCurves.h"
#include "morph/HexGrid.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>
#ifdef __GLN__ // Currently I've only tested OpenMP on Linux
#include <omp.h>
#endif
#include <hdf5.h>
#include <unistd.h>

#define DEBUG 1
#define DBGSTREAM std::cout
#include <morph/MorphDbg.h>

using std::vector;
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;

using morph::HexGrid;
using morph::ReadCurves;

/*!
 * Define this to use manufactured guidance molecules, A, B and C, in
 * contrast to what Karbowski et al did.
 */
//#define MANUFACTURE_GUIDANCE_MOLECULES 1

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
     * The logpath for this model. Used when saving data out.
     */
    string logpath = "logs";

    /*!
     * Setter which attempts to ensure the path exists.
     */
    void setLogpath (const string p) {
        this->logpath = p;
        // Ensure log directory exists
        morph::Tools::createDir (this->logpath);
    }

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
    /*!
     * Uncoupled concentrations of the factors - i.e. where they start.
     */
    //@{
    vector<double> eta_emx;
    vector<double> eta_pax;
    vector<double> eta_fgf;
    //@}

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
    double Aemx = 1; // 1.34 in paper
    double Apax = 1; // 1.4 in paper
    double Afgf = 1; // 0.9 in paper

    // These are scaled rougly in proportion with the values in
    // Karb2004. I have about a 1mm long cortex, so their Chis are
    // divided by 40 to get these values.
    double Chiemx = 0.64; // 25.6/40
    double Chipax = 0.68; // 27.3/40
    double Chifgf = 0.66; // 26.4/40

    // See Karb2004 Fig 8 - the arealization duplication figure
    //@{
    bool useSecondFgfSource = false;
    double Afgfprime = 1.5;
    double Chifgfprime = 0.3; // 12/40
    //@}

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
    float dirpax = 0;
    float dirfgf = 0;
    //@}

    //@} end factor expression dynamics parameters

    /*!
     * Params used in the calculation of rhoA, rhoB and rhoC from the
     * final eta, pax and fgf expression levels.
     */
    //@{
    double sigmaA = 0.2;
    double sigmaB = 0.3;
    double sigmaC = 0.2;

    double kA = 0.34;
    double kB = 0.9;
    double kC = 0.3;

    double theta1 = 0.4;  // 0.77 originally, in the paper
    double theta2 = 0.5;  // 0.5
    double theta3 = 0.39; // 0.39
    double theta4 = 0.3;  // 0.08
    //@}

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
     * Store Hex positions for saving.
     */
    vector<float> hgvx;
    vector<float> hgvy;

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
     * A frame number, incremented when an image is plotted to a PNG file.
     */
    unsigned int frameN = 0;

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
        for (unsigned int i=0; i<this->N; ++i) {
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
     *
     * I apply a sigmoid to the boundary hexes, so that the noise
     * drops away towards the edge of the domain.
     */
    void noiseify_vector_vector (vector<vector<double> >& vv) {
        double randNoiseOffset = 0.8;
        double randNoiseGain = 0.1;
        for (unsigned int i = 0; i<this->N; ++i) {
            for (auto h : this->hg->hexen) {
                // boundarySigmoid. Jumps sharply (100, larger is
                // sharper) over length scale 0.05 to 1. So if
                // distance from boundary > 0.05, noise has normal
                // value. Close to boundary, noise is less.
                vv[i][h.vi] = morph::Tools::randDouble() * randNoiseGain + randNoiseOffset;
                if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                    double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary-0.02)) );
                    vv[i][h.vi] = vv[i][h.vi] * bSig;
                }
            }
        }
    }

    /*!
     * Initialise HexGrid, variables and parameters. Carry out
     * one-time computations of the model.
     */
    void init (vector<morph::Gdisplay>& displays, bool useSavedGenetics = false) {

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
        // Save hex positions in vectors for datafile saving
        for (auto h : this->hg->hexen) {
            this->hgvx.push_back (h.x);
            this->hgvy.push_back (h.y);
        }

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

        // The gamma values - notice the symmetry here.
        // red
        this->gammaA[0] = -2.0;
        this->gammaB[0] =  0.5;
        this->gammaC[0] =  0.5;
        // yellow
        this->gammaA[1] = -2.0;
        this->gammaB[1] = -2.0;
        this->gammaC[1] =  0.5;
        // green
        this->gammaA[2] =  0.5;
        this->gammaB[2] = -2.0;
        this->gammaC[2] =  0.5;
        // blue
        this->gammaA[3] =  0.5;
        this->gammaB[3] = -2.0;
        this->gammaC[3] = -2.0;
        // magenta
        this->gammaA[4] =  0.5;
        this->gammaB[4] =  0.5;
        this->gammaC[4] = -2.0;

        for (unsigned int i=0; i<this->N; ++i) {
            this->alpha[i] = 3;
            this->beta[i] = 3;
        }

        if (useSavedGenetics == false) {
#ifdef MANUFACTURE_GUIDANCE_MOLECULES
            // Construct Gaussian-waves rather than doing the full-Karbowski shebang.
            this->makeupChemoAttractants();
#else
            // Generate the assumed uncoupled concentrations of growth/transcription factors
            this->createFactorInitialConc (this->diremx, this->Aemx, this->Chiemx, this->eta_emx);
            this->createFactorInitialConc (this->dirpax, this->Apax, this->Chipax, this->eta_pax);

            // Should we use a second Fgf source, as in the Karbowski paper (Fig 8)?
            if (this->useSecondFgfSource) {
                this->createFactorInitialConc (this->dirfgf, this->Afgf, this->Afgfprime,
                                               this->Chifgf, this->Chifgfprime, this->eta_fgf);
            } else {
                this->createFactorInitialConc (this->dirfgf, this->Afgf, this->Chifgf, this->eta_fgf);
            }
            // Run the expression dynamics, showing images as we go.
            this->runExpressionDynamics (displays);
            // Can now populate rhoA, rhoB and rhoC according to the paper.
            this->populateChemoAttractants();
#endif
            this->plotchemo (displays);

            // Save pngs of the factors and guidance expressions.
            displays[0].saveImage (this->logpath + "/factors.png");
            displays[1].saveImage (this->logpath + "/guidance.png");

            // Compute gradients of guidance molecule concentrations once only
            this->spacegrad2D (this->rhoA, this->grad_rhoA);
            this->spacegrad2D (this->rhoB, this->grad_rhoB);
            this->spacegrad2D (this->rhoC, this->grad_rhoC);

            // Having computed gradients, build this->g; has
            // to be done once only. Note that a sigmoid is applied so
            // that g(x) drops to zero around the boundary of the domain.
            for (unsigned int i=0; i<this->N; ++i) {
                for (auto h : this->hg->hexen) {
                    // Sigmoid/logistic fn params: 100 sharpness, 0.02 dist offset from boundary
                    double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary-0.02)) );
                    this->g[i][0][h.vi] = (this->gammaA[i] * this->grad_rhoA[0][h.vi]
                                           + this->gammaB[i] * this->grad_rhoB[0][h.vi]
                                           + this->gammaC[i] * this->grad_rhoC[0][h.vi]) * bSig;
                    this->g[i][1][h.vi] = (this->gammaA[i] * this->grad_rhoA[1][h.vi]
                                           + this->gammaB[i] * this->grad_rhoB[1][h.vi]
                                           + this->gammaC[i] * this->grad_rhoC[1][h.vi]) * bSig;
                }
            }

            // Save that data out
            this->saveFactorExpression();

        } else {
            // Load the data from files
            this->loadFactorExpression();
        }
    }

    /*!
     * HDF5 file saving/loading methods
     */
    //@{

    /*!
     * Save the c variable.
     */
    void saveC (void) {
        stringstream fname;
        fname << this->logpath << "/c_";
        fname.width(5);
        fname.fill('0');
        fname << this->stepCount << ".h5";
        hid_t file_id = H5Fcreate (fname.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        for (unsigned int i = 0; i<this->N; ++i) {
            stringstream path;
            path << "/c" << i;
            this->add_double_vector_to_hdf5file (file_id, path.str().c_str(), this->c[i]);
        }
        this->saveHexPositions (file_id);
        herr_t status = H5Fclose (file_id);
        if (status) {
            cerr << "status: " << status << endl;
        }
    }

    void saveHexPositions (const hid_t& file_id) {
        this->add_float_vector_to_hdf5file (file_id, "/x", this->hgvx);
        this->add_float_vector_to_hdf5file (file_id, "/y", this->hgvy);
        // And hex to hex distance
        this->add_double_to_hdf5file (file_id, "/d", this->d);
    }

    /*!
     * Save the results of running createFactorInitialConc(),
     * runExpressionDynamics() and populateChemoAttractants().
     */
    void saveFactorExpression (void) {
        string fname = this->logpath + "/factorexpression.h5";
        cout << "Saving to file " << fname << endl;
        hid_t file_id = H5Fcreate (fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // Save initial factor calculations - required in
        // createFactorInitialConc() and runExpressionDynamics()
        this->add_double_to_hdf5file (file_id, "/Aemx", this->Aemx);
        this->add_double_to_hdf5file (file_id, "/Apax", this->Apax);
        this->add_double_to_hdf5file (file_id, "/Afgf", this->Afgf);

        this->add_double_to_hdf5file (file_id, "/Chiemx", this->Chiemx);
        this->add_double_to_hdf5file (file_id, "/Chipax", this->Chipax);
        this->add_double_to_hdf5file (file_id, "/Chifgf", this->Chifgf);

        this->add_double_to_hdf5file (file_id, "/tau_emx", this->tau_emx);
        this->add_double_to_hdf5file (file_id, "/tau_pax", this->tau_pax);
        this->add_double_to_hdf5file (file_id, "/tau_fgf", this->tau_fgf);

        this->add_float_to_hdf5file (file_id, "/diremx", this->diremx);
        this->add_float_to_hdf5file (file_id, "/dirpax", this->dirpax);
        this->add_float_to_hdf5file (file_id, "/dirfgf", this->dirfgf);

        this->add_double_to_hdf5file (file_id, "/v1", this->v1);
        this->add_double_to_hdf5file (file_id, "/v2", this->v2);

        this->add_double_to_hdf5file (file_id, "/w1", this->w1);
        this->add_double_to_hdf5file (file_id, "/w2", this->w2);

        // Signalling molecule expression levels
        this->add_double_vector_to_hdf5file (file_id, "/emx", this->emx);
        this->add_double_vector_to_hdf5file (file_id, "/pax", this->pax);
        this->add_double_vector_to_hdf5file (file_id, "/fgf", this->fgf);

        this->add_double_vector_to_hdf5file (file_id, "/eta_emx", this->eta_emx);
        this->add_double_vector_to_hdf5file (file_id, "/eta_pax", this->eta_pax);
        this->add_double_vector_to_hdf5file (file_id, "/eta_fgf", this->eta_fgf);

        // parameters and vars for populateChemoAttractants
        this->add_double_to_hdf5file (file_id, "/sigmaA", this->sigmaA);
        this->add_double_to_hdf5file (file_id, "/sigmaB", this->sigmaB);
        this->add_double_to_hdf5file (file_id, "/sigmaC", this->sigmaC);

        this->add_double_to_hdf5file (file_id, "/kA", this->kA);
        this->add_double_to_hdf5file (file_id, "/kB", this->kB);
        this->add_double_to_hdf5file (file_id, "/kC", this->kC);

        this->add_double_to_hdf5file (file_id, "/theta1", this->theta1);
        this->add_double_to_hdf5file (file_id, "/theta2", this->theta2);
        this->add_double_to_hdf5file (file_id, "/theta3", this->theta3);
        this->add_double_to_hdf5file (file_id, "/theta4", this->theta4);

        // The axon guidance molecule expression levels
        this->add_double_vector_to_hdf5file (file_id, "/rhoA", this->rhoA);
        this->add_double_vector_to_hdf5file (file_id, "/rhoB", this->rhoB);
        this->add_double_vector_to_hdf5file (file_id, "/rhoC", this->rhoC);

        // And gradient thereof
        this->add_double_vector_to_hdf5file (file_id, "/grad_rhoA_x", this->grad_rhoA[0]);
        this->add_double_vector_to_hdf5file (file_id, "/grad_rhoA_y", this->grad_rhoA[1]);
        this->add_double_vector_to_hdf5file (file_id, "/grad_rhoB_x", this->grad_rhoB[0]);
        this->add_double_vector_to_hdf5file (file_id, "/grad_rhoB_y", this->grad_rhoB[1]);
        this->add_double_vector_to_hdf5file (file_id, "/grad_rhoC_x", this->grad_rhoC[0]);
        this->add_double_vector_to_hdf5file (file_id, "/grad_rhoC_y", this->grad_rhoC[1]);

        // g - the guidance molecular modifier on a.
        this->add_double_vector_to_hdf5file (file_id, "/g_0_x", this->g[0][0]);
        this->add_double_vector_to_hdf5file (file_id, "/g_0_y", this->g[0][1]);
        this->add_double_vector_to_hdf5file (file_id, "/g_1_x", this->g[1][0]);
        this->add_double_vector_to_hdf5file (file_id, "/g_1_y", this->g[1][1]);
        this->add_double_vector_to_hdf5file (file_id, "/g_2_x", this->g[2][0]);
        this->add_double_vector_to_hdf5file (file_id, "/g_2_y", this->g[2][1]);
        this->add_double_vector_to_hdf5file (file_id, "/g_3_x", this->g[3][0]);
        this->add_double_vector_to_hdf5file (file_id, "/g_3_y", this->g[3][1]);
        this->add_double_vector_to_hdf5file (file_id, "/g_4_x", this->g[4][0]);
        this->add_double_vector_to_hdf5file (file_id, "/g_4_y", this->g[4][1]);

        this->saveHexPositions (file_id);

        herr_t status = H5Fclose (file_id);
        if (status) {
            cerr << "status: " << status << endl;
        }
    }

    /*!
     * Load the results of running createFactorInitialConc(),
     * runExpressionDynamics() and populateChemoAttractants().
     */
    void loadFactorExpression (void) {
        throw runtime_error ("Implement loadFactorExpression().");
    }

    //@} // HDF5

    /*!
     * Computation methods
     */
    //@{

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
        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            Hex* h = this->hg->vhexen[hi];
            n[h->vi] = 0; // whoops forgot this!
            for (unsigned int i=0; i<N; ++i) {
                n[h->vi] += c[i][h->vi];
            }
            csum += c[0][h->vi];
            n[h->vi] = 1. - n[h->vi];
            nsum += n[h->vi];
        }
        if (this->stepCount % 20 == 0) {
            DBG ("sum of all n is " << nsum);
            DBG ("sum of all c for i=0 is " << csum);
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

            DBG2 ("(a) alpha_c_beta_na["<<i<<"][0] = " << this->alpha_c_beta_na[i][0]);

            // Runge-Kutta integration for A
            vector<double> q(this->nhex, 0.0);
            this->compute_divJ (a[i], i); // populates divJ[i]
            DBG2 ("Computing divJ, divJ[0][0]: " << divJ[0][0]);
            vector<double> k1(this->nhex, 0.0);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k1[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                q[h] = this->a[i][h] + k1[h] * halfdt;
            }
            DBG2 ("(a) After RK stage 1, q[0]: " << q[0]);

            vector<double> k2(this->nhex, 0.0);
            this->compute_divJ (q, i);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k2[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                q[h] = this->a[i][h] + k2[h] * halfdt; // Kaboom!
            }
            DBG2 ("(a) After RK stage 2, q[0]:" << q[0] << " from a["<<i<<"][0]:" << a[i][0] << " divj["<<i<<"][0]:" << divJ[i][0] << " k2[0]:" << k2[0]);

            vector<double> k3(this->nhex, 0.0);
            this->compute_divJ (q, i);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k3[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                q[h] = this->a[i][h] + k3[h] * dt;
            }
            DBG2 ("(a) After RK stage 3, q[0]: " << q[0]);

            vector<double> k4(this->nhex, 0.0);
            this->compute_divJ (q, i);
            for (unsigned int h=0; h<this->nhex; ++h) {
                k4[h] = this->divJ[i][h] + this->alpha_c_beta_na[i][h];
                a[i][h] += (k1[h] + 2.0 * (k2[h] + k3[h]) + k4[h]) * sixthdt;
            }
            DBG2 ("(a) After RK stage 4, a[" << i << "][0]: " << a[i][0]);

            DBG2("(a) Debug a["<<i<<"]");
        }

        // 3. Do integration of c
        #pragma omp parallel for
        for (unsigned int i=0; i<this->N; ++i) {

            for (unsigned int h=0; h<nhex; h++) {
                this->betaterm[i][h] = beta[i] * n[h] * pow (a[i][h], k);
            }
            DBG2 ("(c) betaterm[" << i << "][0]: " << betaterm[i][0]);

            // Runge-Kutta integration for C (or ci)
            vector<double> q(nhex,0.);
            vector<double> k1 = compute_dci_dt (c[i], i);
            for (unsigned int h=0; h<nhex; h++) {
                q[h] = c[i][h] + k1[h] * halfdt;
            }
            DBG2 ("(c) After RK stage 1, q[0]: " << q[0]);

            vector<double> k2 = compute_dci_dt (q, i);
            for (unsigned int h=0; h<nhex; h++) {
                q[h] = c[i][h] + k2[h] * halfdt;
            }
            DBG2 ("(c) After RK stage 2, q[0]: " << q[0]);

            vector<double> k3 = compute_dci_dt (q, i);
            for (unsigned int h=0; h<nhex; h++) {
                q[h] = c[i][h] + k3[h] * dt;
            }
            DBG2 ("(c) After RK stage 3, q[0]: " << q[0]);

            vector<double> k4 = compute_dci_dt (q, i);
            for (unsigned int h=0; h<nhex; h++) {
                c[i][h] += (k1[h]+2. * (k2[h] + k3[h]) + k4[h]) * sixthdt;
            }
            DBG2 ("(c) After RK stage 4, c["<<i<<"][0]: " << c[i][0]);

            DBG2("(c) Debug c["<<i<<"]");
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

        // Note - East is positive x; North is positive y. Does this match how it's drawn in the display??
        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            Hex* h = this->hg->vhexen[hi];

            gradf[0][h->vi] = 0.0;
            gradf[1][h->vi] = 0.0;

            DBG2 ("(h->ri,h->gi): (" << h->ri << "," << h->gi << ")");
            // Find x gradient
            if (h->has_ne && h->has_nw) {
                DBG2 ("x case 1 f[h->ne]: " << f[h->ne->vi] << " - f[h->nw]" << f[h->nw->vi] << "/ h->d*2: " << (double)h->d * 2.0);
                gradf[0][h->vi] = (f[h->ne->vi] - f[h->nw->vi]) / ((double)h->d * 2.0);
            } else if (h->has_ne) {
                DBG2 ("x case 2 f[h->ne]: " << f[h->ne->vi] << " - f[h]" << f[h->vi] << "/ h->d: " << (double)h->d);
                gradf[0][h->vi] = (f[h->ne->vi] - f[h->vi]) / (double)h->d;
            } else if (h->has_nw) {
                DBG2 ("x case 3 f[h]: " << f[h->vi] << " - f[h->nw]" << f[h->nw->vi] << "/ h->d: " << (double)h->d);
                gradf[0][h->vi] = (f[h->vi] - f[h->nw->vi]) / (double)h->d;
            } else {
                // zero gradient in x direction as no neighbours in
                // those directions? Or possibly use the average of
                // the gradient between the nw,ne and sw,se neighbours
            }

            // Find y gradient
            if (h->has_nnw && h->has_nne && h->has_nsw && h->has_nse) {
                // Full complement. Compute the mean of the nse->nne and nsw->nnw gradients
#ifdef DEBUG2
                if (h->vi == 0) {
                    double _d = (double)h->getV();
                    double _td = (double)h->getTwoV();
                    DBG2 ("y case 1. getV: " << _d << " getTwoV: " << _td);
                }
#endif
                gradf[1][h->vi] = ((f[h->nne->vi] - f[h->nse->vi]) + (f[h->nnw->vi] - f[h->nsw->vi])) / (double)h->getV();

            } else if (h->has_nnw && h->has_nne ) {
                //if (h->vi == 0) { DBG ("y case 2"); }
                gradf[1][h->vi] = ( (f[h->nne->vi] + f[h->nnw->vi]) / 2.0 - f[h->vi]) / (double)h->getV();

            } else if (h->has_nsw && h->has_nse) {
                //if (h->vi == 0) { DBG ("y case 3"); }
                gradf[1][h->vi] = (f[h->vi] - (f[h->nse->vi] + f[h->nsw->vi]) / 2.0) / (double)h->getV();

            } else if (h->has_nnw && h->has_nsw) {
                //if (h->vi == 0) { DBG ("y case 4"); }
                gradf[1][h->vi] = (f[h->nnw->vi] - f[h->nsw->vi]) / (double)h->getTwoV();

            } else if (h->has_nne && h->has_nse) {
                //if (h->vi == 0) { DBG ("y case 5"); }
                gradf[1][h->vi] = (f[h->nne->vi] - f[h->nse->vi]) / (double)h->getTwoV();
            } else {
                // Leave grady at 0
            }

            //if (h->vi == 0) {
            //    DBG ("gradf[0/1][0]: " << gradf[0][0] << "," << gradf[1][0]);
            //}
        }
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
     * spacegrad2D().  Output: this->divJ
     *
     * Stable with dt = 0.0001;
     */
    void compute_divJ (vector<double>& fa, unsigned int i) {

        // Three terms to compute; see Eq. 14 in methods_notes.pdf

        // Compute gradient of a_i(x), for use computing the third term, below.
        this->spacegrad2D (fa, this->grad_a[i]);

        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {

            Hex* h = this->hg->vhexen[hi];
            // 1. The D Del^2 a_i term
            // Compute the sum around the neighbours
            double thesum = -6 * fa[h->vi];
            if (h->has_ne) {
                thesum += fa[h->ne->vi];
            } else {
                // Apply boundary condition
            }
            if (h->has_nne) {
                thesum += fa[h->nne->vi];
            } else {
                thesum += fa[h->vi]; // A ghost neighbour-east with same value as Hex_0
            }
            if (h->has_nnw) {
                thesum += fa[h->nnw->vi];
            } else {
                thesum += fa[h->vi];
            }
            if (h->has_nw) {
                thesum += fa[h->nw->vi];
            } else {
                thesum += fa[h->vi];
            }
            if (h->has_nsw) {
                thesum += fa[h->nsw->vi];
            } else {
                thesum += fa[h->vi];
            }
            if (h->has_nse) {
                thesum += fa[h->nse->vi];
            } else {
                thesum += fa[h->vi];
            }
            // Multiply bu 2D/3d^2
            double term1 = (this->D * 2) / (3 * this->d * this->d) * thesum;

            // 2. The a div(g) term. Two sums for this.
            double term2 = 0.0;
            // First sum
            if (h->has_ne) {
                term2 += /*cos (0)*/ (this->g[i][0][h->ne->vi] + this->g[i][0][h->vi]);
            } else {
                // Boundary condition _should_ be satisfied by
                // sigmoidal roll-off of g towards the boundary, so
                // add only g[i][0][h->vi]
                term2 += /*cos (0)*/ (this->g[i][0][h->vi]);
            }
            if (h->has_nne) {
                term2 += /*cos (60)*/ 0.5 * (this->g[i][0][h->nne->vi] + this->g[i][0][h->vi]);
            } else {
                term2 += /*cos (60)*/ 0.5 * (this->g[i][0][h->vi]);
            }
            if (h->has_nnw) {
                term2 -= /*cos (120)*/ 0.5 * (this->g[i][0][h->nnw->vi] + this->g[i][0][h->vi]);
            } else {
                term2 -= /*cos (120)*/ 0.5 * (this->g[i][0][h->vi]);
            }
            if (h->has_nw) {
                term2 -= /*cos (180)*/ (this->g[i][0][h->nw->vi] + this->g[i][0][h->vi]);
            } else {
                term2 -= /*cos (180)*/ (this->g[i][0][h->vi]);
            }
            if (h->has_nsw) {
                term2 -= /*cos (240)*/ 0.5 * (this->g[i][0][h->nsw->vi] + this->g[i][0][h->vi]);
            } else {
                term2 -= /*cos (240)*/ 0.5 * (this->g[i][0][h->vi]);
            }
            if (h->has_nse) {
                term2 += /*cos (300)*/ 0.5 * (this->g[i][0][h->nse->vi] + this->g[i][0][h->vi]);
            } else {
                term2 += /*cos (300)*/ 0.5 * (this->g[i][0][h->vi]);
            }
            // 2nd sum
            //term2 += sin (0) * (this->g[i][1][h->ne->vi] + this->g[i][1][h->vi]);
            if (h->has_nne) {
                term2 += /*sin (60)*/ R3_OVER_2 * (this->g[i][1][h->nne->vi] + this->g[i][1][h->vi]);
            } else {
                term2 += /*sin (60)*/ R3_OVER_2 * (this->g[i][1][h->vi]);
            }
            if (h->has_nnw) {
                term2 += /*sin (120)*/ R3_OVER_2 * (this->g[i][1][h->nnw->vi] + this->g[i][1][h->vi]);
            } else {
                term2 += /*sin (120)*/ R3_OVER_2 * (this->g[i][1][h->vi]);
            }
            //term2 += sin (180) * (this->g[i][1][h->nw->vi] + this->g[i][1][h->vi]);
            if (h->has_nsw) {
                term2 -= /*sin (240)*/ R3_OVER_2 * (this->g[i][1][h->nsw->vi] + this->g[i][1][h->vi]);
            } else {
                term2 -= /*sin (240)*/ R3_OVER_2 * (this->g[i][1][h->vi]);
            }
            if (h->has_nse) {
                term2 -= /*sin (300)*/ R3_OVER_2 * (this->g[i][1][h->nse->vi] + this->g[i][1][h->vi]);
            } else {
                term2 -= /*sin (300)*/ R3_OVER_2 * (this->g[i][1][h->vi]);
            }

            term2 /= (3.0 * this->d);
            term2 *= fa[h->vi];

            // 3. Third term is this->g . grad a_i. Should not
            // contribute to J, as g(x) decays towards boundary.
            double term3 = this->g[i][0][h->vi] * this->grad_a[i][0][h->vi]
                + this->g[i][1][h->vi] * this->grad_a[i][1][h->vi];

            this->divJ[i][h->vi] = term1 + term2 + term3;
        }
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
            // x_ here is x from the Hex.
            result[h.vi] = Afac * exp (-(x_ * x_) / (chifac * chifac));
        }
    }

    /*!
     * A special version of createFactorInitialConc which introduces a
     * second source for the factor - see Karb2004 Fig 8.
     */
    void createFactorInitialConc (float phi,
                                  double Afac, double Afacprime,
                                  double chifac, double chifacprime,
                                  vector<double>& result) {

        // Work in a co-ordinate system rotated by phi radians, called x_, y_
        double x_ = 0.0;
        double x_rev = 0.0;

        double cosphi = (double) cos (phi);
        double sinphi = (double) sin (phi);

        // Get minimum x and maximum x in the rotated co-ordinate system.
        double x_min_ = this->hg->getXmin (phi);
        double x_max_ = this->hg->getXmax (phi);

        for (auto h : this->hg->hexen) {
            // Rotate x, then offset by the minimum along that line
            x_ = (h.x * cosphi) + (h.y * sinphi) - x_min_;
            x_rev = x_max_ - (h.x * cosphi) + (h.y * sinphi);
            // x_ here is x from the Hex.
            result[h.vi] = Afac * exp (-(x_ * x_) / (chifac * chifac)) + Afacprime * exp (-(x_rev * x_rev) / (chifacprime * chifacprime));
        }
    }

    /*!
     * Execute Eqs 5-7 of the Karbowski paper to find the steady state
     * of the growth/transcription factors after they have interacted
     * for a long time.
     */
    void runExpressionDynamics (vector<morph::Gdisplay>& displays) {
        #pragma omp parallel for
        for (unsigned int t=0; t<300000; ++t) { // 300000 matches Stuart's 1D Karbowski model
            for (auto h : this->hg->hexen) {
                emx[h.vi] += tau_emx * (-emx[h.vi] + eta_emx[h.vi] / (1. + w2 * fgf[h.vi] + v2 * pax[h.vi]));
                pax[h.vi] += tau_pax * (-pax[h.vi] + eta_pax[h.vi] / (1. + v1 * emx[h.vi]));
                fgf[h.vi] += tau_fgf * (-fgf[h.vi] + eta_fgf[h.vi] / (1. + w1 * emx[h.vi]));
            }
        }
        this->plotexpression (displays);
    }

    /*!
     * Using this->fgf and some hyperbolic tangents, populate rhoA/B/C
     */
    void populateChemoAttractants (void) {
        // chemo-attraction gradient. cf Fig 1 of Karb 2004
        #pragma omp parallel for
        for (unsigned int h=0; h<this->nhex; ++h) {
            this->rhoA[h] = (kA/2.)*(1.+tanh((fgf[h]-theta1)/sigmaA));
            this->rhoB[h] = (kB/2.)*(1.+tanh((theta2-fgf[h])/sigmaB))*(kB/2.)*(1.+tanh((fgf[h]-theta3)/sigmaB));
            this->rhoC[h] = (kC/2.)*(1.+tanh((theta4-fgf[h])/sigmaC));
        }
    }

#ifdef MANUFACTURE_GUIDANCE_MOLECULES
    /*!
     * Generate Gaussian profiles for the chemo-attractants.
     *
     * Instead of using the Karbowski equations, just make some
     * gaussian 'waves'
     */
    void makeupChemoAttractants (void) {

        // Potentially alter angle of Gaussian wave:
        double phi = M_PI*0.0;

        // Gaussian params:
        double sigma = 0.1;
        double gain = 1;

        // Positions of the bumps:
        double xoffA = -0.35;
        double xoffB = 0.0;
        double xoffC = 0.4;

        double cosphi = (double) cos (phi);
        double sinphi = (double) sin (phi);

        for (auto h : this->hg->hexen) {
            double x_ = (h.x * cosphi) + (h.y * sinphi);
            this->rhoA[h.vi] = gain * exp(-((x_-xoffA)*(x_-xoffA)) / sigma);
            this->rhoB[h.vi] = gain * exp(-((x_-xoffB)*(x_-xoffB)) / sigma);
            this->rhoC[h.vi] = gain * exp(-((x_-xoffC)*(x_-xoffC)) / sigma);
        }
    }
#endif

    /*!
     * Plotting code
     */
    //@{
    /*!
     * Plot the system on @a disps
     */
    void plot (vector<morph::Gdisplay>& disps, bool savePngs = false) {
        this->plot_f (this->a, disps[2]);
        this->plot_f (this->c, disps[3]);

        this->plot_contour (this->c, disps[4], 0.75);

        if (savePngs) {
            // a
            stringstream ff1;
            ff1 << this->logpath << "/a_";
            ff1 << std::setw(5) << std::setfill('0') << frameN;
            ff1 << ".png";
            disps[2].saveImage (ff1.str());
            // c
            stringstream ff2;
            ff2 << this->logpath << "/c_";
            ff2 << std::setw(5) << std::setfill('0') << frameN;
            ff2 << ".png";
            disps[3].saveImage (ff2.str());
            // contours
            stringstream ff3;
            ff3 << this->logpath << "/cntr_";
            ff3 << std::setw(5) << std::setfill('0') << frameN;
            ff3 << ".png";
            disps[4].saveImage (ff3.str());

            frameN++;
        }
    }

    /*!
     * Plot a or c
     */
    void plot_f (vector<vector<double> >& f, morph::Gdisplay& disp) {
        vector<double> fix(3, 0.0);
        vector<double> eye(3, 0.0);
        eye[2] = -0.4;
        vector<double> rot(3, 0.0);

#ifdef INDIVIDUAL_SCALING
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
#else
        // Copies data to plot out of the model
        double maxa = -1e7;
        double mina = +1e7;
        // Determines min and max

        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            Hex* h = this->hg->vhexen[hi];
            if (h->onBoundary() == false) {
                for (unsigned int i = 0; i<this->N; ++i) {
                    if (f[i][h->vi]>maxa) { maxa = f[i][h->vi]; }
                    if (f[i][h->vi]<mina) { mina = f[i][h->vi]; }
                }
            }
        }
        double scalea = 1.0 / (maxa-mina);

        // Determine a colour from min, max and current value
        vector<vector<double> > norm_a;
        this->resize_vector_vector (norm_a);
        for (unsigned int i = 0; i<this->N; ++i) {
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; h++) {
                norm_a[i][h] = fmin (fmax (((f[i][h]) - mina) * scalea, 0.0), 1.0);
            }
        }
#endif

        // Create an offset which we'll increment by the width of the
        // map, starting from the left-most map (f[0])
        float hgwidth = this->hg->getXmax()-this->hg->getXmin();
        array<float,3> offset = { 2*(-hgwidth-(hgwidth/20)), 0.0f, 0.0f };

        // Draw
        disp.resetDisplay (fix, eye, rot);
        for (unsigned int i = 0; i<this->N; ++i) {
            // Note: OpenGL isn't thread-safe, so no omp parallel for here.
            for (auto h : this->hg->hexen) {
                array<float,3> cl_a = morph::Tools::HSVtoRGB ((float)i/(float)this->N,
                                                              norm_a[i][h.vi], 1.0);
                disp.drawHex (h.position(), offset, (h.d/2.0f), cl_a);
            }
            offset[0] += hgwidth + (hgwidth/20);
        }
        disp.redrawDisplay();
    }

    void plot_contour (vector<vector<double> >& f, morph::Gdisplay& disp, double threshold) {

        vector<double> fix(3, 0.0);
        vector<double> eye(3, 0.0);
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

        // Re-normalize
        vector<vector<double> > norm_a;
        this->resize_vector_vector (norm_a);
        for (unsigned int i = 0; i<this->N; ++i) {
            for (unsigned int h=0; h<this->nhex; h++) {
                norm_a[i][h] = fmin (fmax (((f[i][h]) - mina[i]) * scalea[i], 0.0), 1.0);
            }
        }

        // Draw
        double c = threshold;
        disp.resetDisplay (fix, eye, rot);
        array<float,3> cl_blk = {0.0f, 0.0f, 0.0f};
        array<float,3> zero_offset = {0.0f, 0.0f, 0.0f};

        for (unsigned int i = 0; i<this->N; ++i) {
            array<float,3> cl_b = morph::Tools::HSVtoRGB ((double)i/(double)(this->N),1.,1.);
            for (auto h : this->hg->hexen) {
                if (h.onBoundary() == false) {
                    if (norm_a[i][h.vi]<c) {
                        if (norm_a[i][h.ne->vi]>c && h.has_ne) {
                            disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_b, 0);
                        }
                        if (norm_a[i][h.nne->vi]>c && h.has_nne) {
                            disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_b, 1);
                        }
                        if (norm_a[i][h.nnw->vi]>c && h.has_nnw) {
                            disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_b, 2);
                        }
                        if (norm_a[i][h.nw->vi]>c && h.has_nw) {
                            disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_b, 3);
                        }
                        if (norm_a[i][h.nsw->vi]>c && h.has_nsw) {
                            disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_b, 4);
                        }
                        if (norm_a[i][h.nse->vi]>c && h.has_nse) {
                            disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_b, 5);
                        }
                    }

                } else { // h.onBoundary() is true

                    if (!h.has_ne) {
                        disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_blk, 0);
                    }
                    if (!h.has_nne) {
                        disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_blk, 1);
                    }
                    if (!h.has_nnw) {
                        disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_blk, 2);
                    }
                    if (!h.has_nw) {
                        disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_blk, 3);
                    }
                    if (!h.has_nsw) {
                        disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_blk, 4);
                    }
                    if (!h.has_nse) {
                        disp.drawHexSeg (h.position(), zero_offset, (h.d/2.0f), cl_blk, 5);
                    }
                }
            }
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

        // Determines min and max
        double max = -1e7;
        double min = +1e7;
        for (auto h : this->hg->hexen) {
            if (h.onBoundary() == false) {
                if (this->emx[h.vi]>max) { max = this->emx[h.vi]; }
                if (this->emx[h.vi]<min) { min = this->emx[h.vi]; }

                if (this->pax[h.vi]>max) { max = this->pax[h.vi]; }
                if (this->pax[h.vi]<min) { min = this->pax[h.vi]; }

                if (this->fgf[h.vi]>max) { max = this->fgf[h.vi]; }
                if (this->fgf[h.vi]<min) { min = this->fgf[h.vi]; }
            }
        }

        double scale = 1.0 / (max-min);

        // Determine a colour from min, max and current value
        vector<double> norm_emx(this->nhex, 0.0);
        vector<double> norm_pax(this->nhex, 0.0);
        vector<double> norm_fgf(this->nhex, 0.0);
        for (unsigned int h=0; h<this->nhex; h++) {
            norm_emx[h] = fmin (fmax (((this->emx[h]) - min) * scale, 0.0), 1.0);
            norm_pax[h] = fmin (fmax (((this->pax[h]) - min) * scale, 0.0), 1.0);
            norm_fgf[h] = fmin (fmax (((this->fgf[h]) - min) * scale, 0.0), 1.0);
        }

        // Step through vectors or iterate through list? The latter should be just fine here.
        disps[0].resetDisplay (fix, eye, rot);

        // Set offsets for the three maps that we'll plot
        float hgwidth = this->hg->getXmax()-this->hg->getXmin();
        array<float,3> offset1 = { -hgwidth-(hgwidth/20), 0.0f, 0.0f };
        array<float,3> offset2 = { 0.0f, 0.0f, 0.0f };
        array<float,3> offset3 = { hgwidth+(hgwidth/20), 0.0f, 0.0f };

        for (auto h : this->hg->hexen) {
            array<float,3> cl_emx = morph::Tools::HSVtoRGB (0.0, norm_emx[h.vi], 1.0);
            disps[0].drawHex (h.position(), offset1, (h.d/2.0f), cl_emx);
        }
        for (auto h : this->hg->hexen) {
            array<float,3> cl_pax = morph::Tools::HSVtoRGB (0.33, norm_pax[h.vi], 1.0);
            disps[0].drawHex (h.position(), offset2, (h.d/2.0f), cl_pax);
        }
        for (auto h : this->hg->hexen) {
            array<float,3> cl_fgf = morph::Tools::HSVtoRGB (0.66, norm_fgf[h.vi], 1.0);
            disps[0].drawHex (h.position(), offset3, (h.d/2.0f), cl_fgf);
        }

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

        double max = -1e7;
        double min = +1e7;
        // Determines min and max
        for (auto h : this->hg->hexen) {
            if (h.onBoundary() == false) {
                if (this->rhoA[h.vi]>max) { max = this->rhoA[h.vi]; }
                if (this->rhoA[h.vi]<min) { min = this->rhoA[h.vi]; }

                if (this->rhoB[h.vi]>max) { max = this->rhoB[h.vi]; }
                if (this->rhoB[h.vi]<min) { min = this->rhoB[h.vi]; }

                if (this->rhoC[h.vi]>max) { max = this->rhoC[h.vi]; }
                if (this->rhoC[h.vi]<min) { min = this->rhoC[h.vi]; }
            }
        }
        double scale = 1.0 / (max-min);

        // Determine a colour from min, max and current value
        vector<double> norm_rhoA(this->nhex, 0.0);
        vector<double> norm_rhoB(this->nhex, 0.0);
        vector<double> norm_rhoC(this->nhex, 0.0);
        for (unsigned int h=0; h<this->nhex; h++) {
            norm_rhoA[h] = fmin (fmax (((this->rhoA[h]) - min) * scale, 0.0), 1.0);
            norm_rhoB[h] = fmin (fmax (((this->rhoB[h]) - min) * scale, 0.0), 1.0);
            norm_rhoC[h] = fmin (fmax (((this->rhoC[h]) - min) * scale, 0.0), 1.0);
        }

        // Set offsets for the three maps that we'll plot
        float hgwidth = this->hg->getXmax()-this->hg->getXmin();
        array<float,3> offset1 = { -hgwidth-(hgwidth/20), 0.0f, 0.0f };
        array<float,3> offset2 = { 0.0f, 0.0f, 0.0f };
        array<float,3> offset3 = { hgwidth+(hgwidth/20), 0.0f, 0.0f };

        // Step through vectors or iterate through list? The latter should be just fine here.
        disps[1].resetDisplay (fix, eye, rot);
        for (auto h : this->hg->hexen) {
            array<float,3> cl_rhoA = morph::Tools::HSVtoRGB (0.0, norm_rhoA[h.vi], 1.0);
            disps[1].drawHex (h.position(), offset1, (h.d/2.0f), cl_rhoA);
        }

        for (auto h : this->hg->hexen) {
            array<float,3> cl_rhoB = morph::Tools::HSVtoRGB (0.33, norm_rhoB[h.vi], 1.0);
            disps[1].drawHex (h.position(), offset2, (h.d/2.0f), cl_rhoB);
        }

        for (auto h : this->hg->hexen) {
            array<float,3> cl_rhoC = morph::Tools::HSVtoRGB (0.66, norm_rhoC[h.vi], 1.0);
            disps[1].drawHex (h.position(), offset3, (h.d/2.0f), cl_rhoC);
        }
        disps[1].redrawDisplay();
    }

    //@}

    /*!
     * tools for reading and writing to HDF5 files. To be transferred
     * into libmorpholgica in time.
     */
    //@{

    /*!
     * Makes necessary calls to add a double to an HDF5 file store,
     * using path as the name of the variable. Path could be /myvar or
     * /somegroup/myvar, though I think you'd have to have created the
     * group for the latter.
     */
    void add_double_to_hdf5file (hid_t file_id, const char* path, const double& val) {
        hsize_t dim_singleparam[1];
        dim_singleparam[0] = 1;
        hid_t dataspace_id = H5Screate_simple (1, dim_singleparam, NULL);
        hid_t dataset_id = H5Dcreate2 (file_id, path, H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        herr_t status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
        if (status) {
            cerr << "status after H5Dwrite: " << status << endl;
        }
        status = H5Dclose (dataset_id);
        if (status) {
            cerr << "status H5Dclose: " << status << endl;
        }
        status = H5Sclose (dataspace_id);
        if (status) {
            cerr << "status H5Sclose: " << status << endl;
        }
    }

    /*!
     * Makes necessary calls to add a float to an HDF5 file store,
     * using path as the name of the variable.
     */
    void add_float_to_hdf5file (hid_t file_id, const char* path, const float& val) {
        hsize_t dim_singleparam[1];
        herr_t status;
        dim_singleparam[0] = 1;
        hid_t dataspace_id = H5Screate_simple (1, dim_singleparam, NULL);
        hid_t dataset_id = H5Dcreate2 (file_id, path, H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
        if (status) {
            cerr << "status after H5Dwrite: " << status << endl;
        }
        status = H5Dclose (dataset_id);
        if (status) {
            cerr << "status H5Dclose: " << status << endl;
        }
        status = H5Sclose (dataspace_id);
        if (status) {
            cerr << "status H5Sclose: " << status << endl;
        }
    }

    /*!
     * Makes necessary calls to add a vector of doubles to an HDF5
     * file store, using path as the name of the variable.
     */
    void add_double_vector_to_hdf5file (hid_t file_id, const char* path, const vector<double>& vals) {
        hsize_t dim_singleparam[1];
        herr_t status;
        dim_singleparam[0] = vals.size();
        hid_t dataspace_id = H5Screate_simple (1, dim_singleparam, NULL);
        hid_t dataset_id = H5Dcreate2 (file_id, path, H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vals[0]));
        if (status) {
            cerr << "status after H5Dwrite: " << status << endl;
        }
        status = H5Dclose (dataset_id);
        if (status) {
            cerr << "status H5Dclose: " << status << endl;
        }
        status = H5Sclose (dataspace_id);
        if (status) {
            cerr << "status H5Sclose: " << status << endl;
        }
    }

    /*!
     * Makes necessary calls to add a vector of floats to an HDF5
     * file store, using path as the name of the variable.
     */
    void add_float_vector_to_hdf5file (hid_t file_id, const char* path, const vector<float>& vals) {
        hsize_t dim_singleparam[1];
        herr_t status;
        dim_singleparam[0] = vals.size();
        hid_t dataspace_id = H5Screate_simple (1, dim_singleparam, NULL);
        hid_t dataset_id = H5Dcreate2 (file_id, path, H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vals[0]));
        if (status) {
            cerr << "status after H5Dwrite: " << status << endl;
        }
        status = H5Dclose (dataset_id);
        if (status) {
            cerr << "status H5Dclose: " << status << endl;
        }
        status = H5Sclose (dataspace_id);
        if (status) {
            cerr << "status H5Sclose: " << status << endl;
        }
    }
    //@}

#ifdef GAUSSIAN_CODE_NEEDED
    /*!
     * Create a symmetric, 1D Gaussian hill centred at coordinate (x) with
     * width sigma and height gain. Place result into @a result.
     */
    void createGaussian1D (float x, float phi, double gain, double sigma, vector<double>& result) {

        // Once-only parts of the calculation of the Gaussian.
        double root_2_pi = 2.506628275;
        double one_over_sigma_root_2_pi = 1 / sigma * root_2_pi;
        double two_sigma_sq = 2 * sigma * sigma;

        // Gaussian dist. result, and a running sum of the results:
        double gauss = 0.0;

        double cosphi = (double) cos (phi);
        double sinphi = (double) sin (phi);

        // x and y components of the vector from (x,y) to any given Hex.
        float rx = 0.0f, ry = 0.0f;

        // Calculate each element of the kernel:
        for (auto h : this->hg->hexen) {
            rx = x - h.x;
            ry = 0 - h.y;
            double x_ = (rx * cosphi) + (ry * sinphi);
            gauss = gain * (one_over_sigma_root_2_pi
                            * exp ( static_cast<double>(-(x_*x_))
                                    / two_sigma_sq ));
            result[h.vi] = gauss;
            ++k;
        }
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
#endif

}; // RD_2D_Karb

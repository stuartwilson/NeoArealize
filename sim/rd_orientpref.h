#include "morph/display.h"
#include "morph/tools.h"
#include "morph/ReadCurves.h"
#include "morph/HexGrid.h"
#include "morph/HdfData.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>
#ifdef __GLN__
//#include <omp.h>
#endif
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
using std::exp;

using morph::HexGrid;
using morph::ReadCurves;
using morph::HdfData;

/*!
 * Reaction diffusion system; Orientation preference maps
 */
class RD_OrientPref
{
// This code implements my policy of keeping most parameters and
// variables public as this is research code.
public:
    /*!
     * Constants
     */
    //@{
    //! Square root of 3 over 2
    const double R3_OVER_2 = 0.866025403784439;
    //! Square root of 3
    const double ROOT3 = 1.73205080756888;
    //! 2 pi
    const double TWO_PI = 2 * M_PI;
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
     * Track the number of computational steps that we've carried
     * out, so that we can stop after some number of steps.
     */
    unsigned int stepCount = 0;

    /*!
     * A frame number, incremented when an image is plotted to a PNG
     * file.
     */
    unsigned int frameN = 0;

    /*!
     * Model variables
     */
    //@{

    /*!
     * The z variable - a vector field which represents orientation
     * selectivity of neurons. The magnitude of this field represents
     * the selectivity of the neuron, its direction represents the
     * orientation of edges for which it selects.
     */
    array<vector<double>, 2> z;
    /*!
     * z input pattern
     */
    array<vector<double>, 2> sz;

    /*!
     * Afferent input; scalar field.
     */
    vector<double> aff;

    /*!
     * exponential argument - a temporary result
     */
    vector<double> exparg;

    /*!
     * Receptive field centre coordinates for neurons.
     */
    array<vector<double>, 2> r;
    /*!
     * r input pattern
     */
    array<vector<double>, 2> sr;


    /*!
     * Fields used during Runge-Kutta integration (see integrate)
     */
    array<vector<double>, 2> q;
    array<vector<double>, 2> k1;
    array<vector<double>, 2> k2;
    array<vector<double>, 2> k3;
    array<vector<double>, 2> k4;
    array<vector<double>, 2> lap;
    array<vector<double>, 2> argu_r;
    array<vector<double>, 2> argu_z;
    array<vector<double>, 2> diff_r;
    array<vector<double>, 2> diff_z;
    //@} // model variables

    /*!
     * Model parameters
     */
    //@{

    /*!
     * Strength of lateral interactions.
     */
    double eta = 0.0006;

private:
    /*!
     * afferent activation selectivity
     */
    double sigma = 0.5;

    /*!
     * Used many times in code.
     */
    double oneOverTwoSigmaSquared = 1.0 / (2.0 * sigma * sigma);

public:
    /*!
     * Need a setter for sigma, as oneOverTwoSigmaSquared has to be
     * recomputed when sigma is changed.
     */
    void setSigma (double val) {
        this->sigma = val;
        this->oneOverTwoSigmaSquared = 1.0 / (2.0 * sigma * sigma);
    }

    //@} // model parameters

    /*!
     * Simple constructor; no arguments.
     */
    RD_OrientPref (void) {
        this->halfdt = this->dt/2.0;
        this->sixthdt = this->dt/6.0;
    }

    /*!
     * Destructor required to free up HexGrid memory
     */
    ~RD_OrientPref (void) {
        delete (this->hg);
    }

    /*!
     * Resize array (of length 2) of vectors to each be of nhex
     * length.
     */
    void resize_vector_field (array<vector<double>, 2>& av) {
        av[0].resize (this->nhex, 0.0);
        av[1].resize (this->nhex, 0.0);
    }

    void resize_scalar_field (vector<double>& vd) {
        vd.resize (this->nhex, 0.0);
    }

    /*!
     * Initialise this vector of vectors with noise. This is a
     * model-specific function.
     *
     * Optionally apply a sigmoid to the boundary hexes, so that the
     * noise drops away towards the edge of the domain.
     *
     * @param av Fixed size array of two vector<doubles> - this is the
     * input vector field.
     *
     * @param randNoiseOffset After randNoiseGain has been applied to
     * the uniformly sampled random number, apply this additive offset.
     *
     * @param randNoiseGain The uniform random number in range [0 1]
     * is multiplied by this value.
     *
     * @param useSigmoidFalloff If true, attenuate noise as boundary is approached
     *
     * @param sigmoidSharpness If useSigmoidFalloff is true, this is
     * the sharpness of the sigmoid fall-off.
     *
     * @param sigmoidOffset If useSigmoidFalloff is true, this is the
     * length scale over which the fall-off occurs.
     */
    void noiseify_vector_field (array<vector<double>, 2>& av,
                                double randNoiseOffset, double randNoiseGain,
                                bool useSigmoidFalloff,
                                double sigmoidSharpness /* 100.0 */, double sigmoidOffset /*0.02*/) {
        for (unsigned int i = 0; i<2; ++i) {
            #pragma omp parallel for
            for (unsigned int hi=0; hi<this->nhex; ++hi) {
                Hex* h = this->hg->vhexen[hi];
                av[i][h->vi] = morph::Tools::randDouble() * randNoiseGain + randNoiseOffset;
                if (useSigmoidFalloff == true) {
                    // boundarySigmoid. Jumps sharply (100, larger is
                    // sharper) over length scale 0.05 to 1. So if
                    // distance from boundary is greater than about
                    // 2*sigmoidOffset, noise has normal value. Close
                    // to boundary, noise is less.
                    if (h->distToBoundary > -0.5) { // It's possible that distToBoundary was initialised to -1.0
                        double bSig = 1.0 / ( 1.0 + exp (-sigmoidSharpness*(h->distToBoundary-sigmoidOffset)) );
                        av[i][h->vi] = av[i][h->vi] * bSig;
                    }
                }
            }
        }
    }

    /*!
     * Reproduce python code:
     *
     * z = np.exp(2.j*np.random.rand(N,N)*np.pi)
     *
     * Optionally with sigmoid style fall-off as the boundary is
     * approached.
     *
     * @param av Fixed size array of two vector<doubles> - this is the
     * input vector field.
     *
     * @param mult A multiplier to apply to the uniformly sampled
     * random number. This is 2pi for the original example, but can be
     * changed with this argument.
     *
     * @param useSigmoidFalloff If true, attenuate noise as boundary is approached
     *
     * @param sigmoidSharpness If useSigmoidFalloff is true, this is
     * the sharpness of the sigmoid fall-off.
     *
     * @param sigmoidOffset If useSigmoidFalloff is true, this is the
     * length scale over which the fall-off occurs.
     */
    void noiseify_z_field (array<vector<double>, 2>& av, double mult,
                           bool useSigmoidFalloff,
                           double sigmoidSharpness /* 100.0 */, double sigmoidOffset /*0.02*/) {
        double rn = 0;
        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            Hex* h = this->hg->vhexen[hi];

            rn = morph::Tools::randDouble();
            av[0][h->vi] = cos (rn*mult);
            av[1][h->vi] = sin (rn*mult);

            if (useSigmoidFalloff == true) {
                // boundarySigmoid. Jumps sharply (100, larger is
                // sharper) over length scale 0.05 to 1. So if
                // distance from boundary is greater than about
                // 2*sigmoidOffset, noise has normal value. Close
                // to boundary, noise is less.
                if (h->distToBoundary > -0.5) { // It's possible that distToBoundary was initialised to -1.0
                    double bSig = 1.0 / ( 1.0 + exp (-sigmoidSharpness*(h->distToBoundary-sigmoidOffset)) );
                    av[0][h->vi] = av[0][h->vi] * bSig;
                    av[1][h->vi] = av[1][h->vi] * bSig;
                }
            }
        }
    }

    /*!
     * Initialise HexGrid, variables and parameters. Carry out any
     * one-time computations required by the model.
     */
    void init (vector<morph::Gdisplay>& displays, bool useSavedGenetics = false) {

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
        this->resize_vector_field (this->z);
        this->resize_vector_field (this->sz);
        this->resize_vector_field (this->r);
        this->resize_vector_field (this->sr);

        this->resize_vector_field (this->q);
        this->resize_vector_field (this->k1);
        this->resize_vector_field (this->k2);
        this->resize_vector_field (this->k3);
        this->resize_vector_field (this->k4);
        this->resize_vector_field (this->lap);
        this->resize_vector_field (this->argu_r);
        this->resize_vector_field (this->argu_z);
        this->resize_vector_field (this->diff_r);
        this->resize_vector_field (this->diff_z);

        this->resize_scalar_field (this->aff);
        this->resize_scalar_field (this->exparg);

        // Initialise z with noise
        this->noiseify_z_field (this->z, TWO_PI, false, 100.0, 0.02);

        // Initialise r with noise
        this->noiseify_vector_field (this->r,
                                     -0.005, 0.01,
                                     false, 100.0, 0.02);
#if 0
        // debugging shows that initial values are sensible and small.
        string lp = this->logpath;
        this->logpath = "./logs/tmp";
        this->saveState();
        this->logpath = lp;
#endif
    }

    /*!
     * HDF5 file saving/loading methods
     */
    //@{
    void saveHexPositions (HdfData& dat) {
        dat.add_float_vector ("/x", this->hgvx);
        dat.add_float_vector ("/y", this->hgvy);
        // And hex to hex distance:
        dat.add_double ("/d", this->d);
    }

    /*!
     * Save some data like this.
     */
    void saveState (void) {
        string fname = this->logpath + "/pinwheel.h5";
        HdfData data (fname);
        // Save some variables
        data.add_double_vector ("/z_re", this->z[0]);
        data.add_double_vector ("/z_im", this->z[1]);
        data.add_double_vector ("/r_re", this->r[0]);
        data.add_double_vector ("/r_im", this->r[1]);
        this->saveHexPositions (data);
    }
    //@} // HDF5

    /*!
     * Computation methods
     */
    //@{

    /*!
     * Do a step through the model.
     */
    void step (void) {

        stringstream fss;
        fss << this->logpath << "/pinwheel_step_" << stepCount << ".h5";
        HdfData data (fss.str());

        if (this->stepCount % 100 == 0) {
            DBG ("System computed " << this->stepCount << " times so far...");
        }

        // Generate input patter
        this->noiseify_z_field (this->sz, 2.0, false, 100.0, 0.02); // here?
        this->noiseify_vector_field (this->sr,
                                     //-0.5, 1.0,
                                     -0.005, 0.01,
                                     false, 100.0, 0.02);

        data.add_double_vector ("/z_re", this->z[0]);
        data.add_double_vector ("/z_im", this->z[1]);
        data.add_double_vector ("/r_re", this->r[0]);
        data.add_double_vector ("/r_im", this->r[1]);
        data.add_double_vector ("/sz_re", this->sz[0]);
        data.add_double_vector ("/sz_im", this->sz[1]);
        data.add_double_vector ("/sr_re", this->sr[0]);
        data.add_double_vector ("/sr_im", this->sr[1]);

        // Compute afferent response
        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {

            this->diff_r[0][hi] = this->sr[0][hi] - this->r[0][hi];
            this->diff_r[1][hi] = this->sr[1][hi] - this->r[1][hi];
            this->diff_z[0][hi] = this->sz[0][hi] - this->z[0][hi];
            this->diff_z[1][hi] = this->sz[1][hi] - this->z[1][hi];

            // Minimal computation version:
            this->exparg[hi] = - diff_r[0][hi]*diff_r[0][hi] - diff_r[1][hi]*diff_r[1][hi] - diff_z[0][hi]*diff_z[0][hi] - diff_z[1][hi]*diff_z[1][hi];

            if (this->oneOverTwoSigmaSquared * this->exparg[hi] < -708.4 && hi < 5) {
                cerr << "Should get an underflow into this->aff[" << hi << "]" << endl;
            }
            this->aff[hi] = exp (this->oneOverTwoSigmaSquared * this->exparg[hi]);
        }
        data.add_double_vector ("/exparg", this->exparg);
        data.add_double_vector ("/aff", this->aff);

        data.add_double_vector ("/diff_z_re", this->diff_z[0]);
        data.add_double_vector ("/diff_z_im", this->diff_z[1]);
        data.add_double_vector ("/diff_r_re", this->diff_r[0]);
        data.add_double_vector ("/diff_r_im", this->diff_r[1]);

        // Serial sum. Vector sum candidate?
        double affsum = 0;
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            affsum += this->aff[hi];
        }
        data.add_double ("/affsum", affsum);

        // Parallel division-by
        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->aff[hi] /= affsum; // hmm
            // Prepare args for integrate. Vector multiplication candidate!
            this->argu_z[0][hi] = this->aff[hi] * diff_z[0][hi];
            this->argu_z[1][hi] = this->aff[hi] * diff_z[1][hi];
            this->argu_r[0][hi] = this->aff[hi] * diff_r[0][hi];
            this->argu_r[1][hi] = this->aff[hi] * diff_r[1][hi];
        }

        data.add_double_vector ("/argu_z_re", this->argu_z[0]);
        data.add_double_vector ("/argu_z_im", this->argu_z[1]);
        data.add_double_vector ("/argu_r_re", this->argu_r[0]);
        data.add_double_vector ("/argu_r_im", this->argu_r[1]);

        double norm  = 2.0 / (3.0 * this->d * this->d);
        data.add_double ("/norm", norm);
        data.add_double ("/d", this->d);

        this->integrate (this->argu_z, this->z, 0.05, true);
        this->integrate (this->argu_r, this->r, 0.05, false);

        // pi for post integrate
        data.add_double_vector ("/pi_z_re", this->z[0]);
        data.add_double_vector ("/pi_z_im", this->z[1]);
        data.add_double_vector ("/pi_r_re", this->r[0]);
        data.add_double_vector ("/pi_r_im", this->r[1]);

        this->stepCount++;
    }

    void integrate (array<vector<double>, 2>& E, array<vector<double>, 2>& x, double h_ = 0.05, bool first=false) {

        stringstream fss;
        // Hacky way of having two different file names for numerical debugging
        if (first) {
            fss << this->logpath << "/pinwheel_intz_step_" << stepCount << ".h5";
        } else {
            fss << this->logpath << "/pinwheel_intr_step_" << stepCount << ".h5";
        }
        HdfData data (fss.str());
        DBG("stepCount:" << stepCount);
        // Runge-Kutta integration:

        data.add_double_vector ("/x_re", x[0]);
        data.add_double_vector ("/x_im", x[1]);
        data.add_double_vector ("/E_re", E[0]);
        data.add_double_vector ("/E_im", E[1]);

        DBG ("E_re[0]:" << E[0][0]);

        double h2 = h_ * h_;

        this->compute_lapl_cmplx (x, lap);
        DBG("lap_1_re[0]:" << lap[0][0]);
        data.add_double_vector ("/lap_1_re", lap[0]);
        data.add_double_vector ("/lap_1_im", lap[1]);
        #pragma omp parallel for
        for (unsigned int h=0; h<this->nhex; ++h) {
            k1[0][h] = h_ * (E[0][h] + this->eta * lap[0][h]);
            k1[1][h] = h_ * (E[1][h] + this->eta * lap[1][h]);
            q[0][h] = x[0][h] + k1[0][h] * this->halfdt;
            q[1][h] = x[1][h] + k1[1][h] * this->halfdt;
        }

        this->compute_lapl_cmplx (q, lap);
        data.add_double_vector ("/lap_2_re", lap[0]);
        data.add_double_vector ("/lap_2_im", lap[1]);
        #pragma omp parallel for
        for (unsigned int h=0; h<this->nhex; ++h) {
            k2[0][h] = this->halfdt * h2 * (E[0][h] + this->eta * lap[0][h]);
            k2[1][h] = this->halfdt * h2 * (E[1][h] + this->eta * lap[1][h]);
            q[0][h] = x[0][h] + k2[0][h] * this->halfdt;
            q[1][h] = x[1][h] + k2[1][h] * this->halfdt;
        }

        this->compute_lapl_cmplx (q, lap);
        data.add_double_vector ("/lap_3_re", lap[0]);
        data.add_double_vector ("/lap_3_im", lap[1]);
        #pragma omp parallel for
        for (unsigned int h=0; h<this->nhex; ++h) {
            k3[0][h] = this->halfdt * h2 * (E[0][h] + this->eta * lap[0][h]);
            k3[1][h] = this->halfdt * h2 * (E[1][h] + this->eta * lap[1][h]);
            q[0][h] = x[0][h] + k3[0][h];
            q[1][h] = x[1][h] + k3[1][h];
        }

        this->compute_lapl_cmplx (q, lap);
        data.add_double_vector ("/lap_4_re", lap[0]);
        data.add_double_vector ("/lap_4_im", lap[1]);
        #pragma omp parallel for
        for (unsigned int h=0; h<this->nhex; ++h) {
            k4[0][h] = h2 * (E[0][h] + this->eta * lap[0][h]);
            k4[1][h] = h2 * (E[1][h] + this->eta * lap[1][h]);

            // Write final result on this loop back into x
            x[0][h] = x[0][h] + (k1[0][h] + 2.0*(k2[0][h]+k3[0][h]) + k4[0][h]) / 6.0;
            x[1][h] = x[1][h] + (k1[1][h] + 2.0*(k2[1][h]+k3[1][h]) + k4[1][h]) / 6.0;
            if (h < 3) {
                DBG("x[0][" << h << "]=" << x[0][h] << " x[1][" << h << "]=" << x[1][h]);
            }
        }
        data.add_double_vector ("/k1_re", k1[0]);
        data.add_double_vector ("/k1_im", k1[1]);
        data.add_double_vector ("/k2_re", k2[0]);
        data.add_double_vector ("/k2_im", k2[1]);
        data.add_double_vector ("/k3_re", k3[0]);
        data.add_double_vector ("/k3_im", k3[1]);
        data.add_double_vector ("/k4_re", k4[0]);
        data.add_double_vector ("/k4_im", k4[1]);
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

            // Find x gradient
            if (h->has_ne && h->has_nw) {
                gradf[0][h->vi] = (f[h->ne->vi] - f[h->nw->vi]) / ((double)h->d * 2.0);
            } else if (h->has_ne) {
                gradf[0][h->vi] = (f[h->ne->vi] - f[h->vi]) / (double)h->d;
            } else if (h->has_nw) {
                gradf[0][h->vi] = (f[h->vi] - f[h->nw->vi]) / (double)h->d;
            } /*
               * else zero gradient in x direction as no neighbours in
               * those directions? Or possibly use the average of
               * the gradient between the nw,ne and sw,se neighbours
               */

            // Find y gradient
            if (h->has_nnw && h->has_nne && h->has_nsw && h->has_nse) {
                // Full complement. Compute the mean of the nse->nne and nsw->nnw gradients
                gradf[1][h->vi] = ((f[h->nne->vi] - f[h->nse->vi]) + (f[h->nnw->vi] - f[h->nsw->vi])) / (double)h->getV();

            } else if (h->has_nnw && h->has_nne ) {
                gradf[1][h->vi] = ( (f[h->nne->vi] + f[h->nnw->vi]) / 2.0 - f[h->vi]) / (double)h->getV();

            } else if (h->has_nsw && h->has_nse) {
                gradf[1][h->vi] = (f[h->vi] - (f[h->nse->vi] + f[h->nsw->vi]) / 2.0) / (double)h->getV();

            } else if (h->has_nnw && h->has_nsw) {
                gradf[1][h->vi] = (f[h->nnw->vi] - f[h->nsw->vi]) / (double)h->getTwoV();

            } else if (h->has_nne && h->has_nse) {
                gradf[1][h->vi] = (f[h->nne->vi] - f[h->nse->vi]) / (double)h->getTwoV();
            } // else leave grady at 0
        }
    }

    /*!
     * Computes the Laplacian (the divergence of the gradient) of the
     * scalar field fa, placing the result in the scalar field
     * laplace.
     */
    void compute_lapl (vector<double>& fa, vector<double>& laplace) {

        // The normalisation comes from the area of the hex. See methods_notes.pdf.
        double norm  = 2.0 / (3.0 * this->d * this->d);

        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {

            Hex* h = this->hg->vhexen[hi];

            // Compute the sum around the neighbours
            double thesum = -6 * fa[h->vi];
            if (h->has_ne) {
                thesum += fa[h->ne->vi];
            } else {
                thesum += fa[h->vi]; // A ghost neighbour-east with same value as Hex_0
            }
            if (h->has_nne) {
                thesum += fa[h->nne->vi];
            } else {
                thesum += fa[h->vi];
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

            laplace[h->vi] = norm * thesum;
        }
    }

    /*!
     * Like compute_lapl, but with the input and output treated as
     * "scalar fields of complex numbers".
     */
    void compute_lapl_cmplx (array<vector<double>, 2>& fa, array<vector<double>, 2>& laplace) {

        double norm  = 2.0 / (3.0 * this->d * this->d);
        // DBG ("norm: " << norm); it's 6666.67

        double sum_re_tot = 0.0;
        double sum_im_tot = 0.0;

        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {

            Hex* h = this->hg->vhexen[hi];

            // Compute the sum around the neighbours
            double sum_re = -6 * fa[0][h->vi];
            double sum_im = -6 * fa[1][h->vi];
            if (h->has_ne) {
                sum_re += fa[0][h->ne->vi];
                sum_im += fa[1][h->ne->vi];
            } else {
                sum_re += fa[0][h->vi]; // A ghost neighbour-east with
                                        // same value as Hex_0 - this
                                        // is therefore no-flux
                                        // boundary condition. Not
                                        // quite the same as the
                                        // original model. FIXME
                                        // perhaps.
                sum_im += fa[1][h->vi];
            }
            if (h->has_nne) {
                sum_re += fa[0][h->nne->vi];
                sum_im += fa[1][h->nne->vi];
            } else {
                sum_re += fa[0][h->vi];
                sum_im += fa[1][h->vi];
            }
            if (h->has_nnw) {
                sum_re += fa[0][h->nnw->vi];
                sum_im += fa[1][h->nnw->vi];
            } else {
                sum_re += fa[0][h->vi];
                sum_im += fa[1][h->vi];
            }
            if (h->has_nw) {
                sum_re += fa[0][h->nw->vi];
                sum_im += fa[1][h->nw->vi];
            } else {
                sum_re += fa[0][h->vi];
                sum_im += fa[1][h->vi];
            }
            if (h->has_nsw) {
                sum_re += fa[0][h->nsw->vi];
                sum_im += fa[1][h->nsw->vi];
            } else {
                sum_re += fa[0][h->vi];
                sum_im += fa[1][h->vi];
            }
            if (h->has_nse) {
                sum_re += fa[0][h->nse->vi];
                sum_im += fa[1][h->nse->vi];
            } else {
                sum_re += fa[0][h->vi];
                sum_im += fa[1][h->vi];
            }

            laplace[0][h->vi] = norm * sum_re;
            laplace[1][h->vi] = norm * sum_im;

            // Just for debugging, will break in parallel
            sum_re_tot += laplace[0][h->vi];
            sum_im_tot += laplace[1][h->vi];
        }

        DBG ("Real laplacian total: " << sum_re_tot);
        DBG ("Imag laplacian total: " << sum_im_tot);
    }

    /*!
     * Plotting code
     */
    //@{
    /*!
     * Plot the system on @a disps
     */
    void plot (vector<morph::Gdisplay>& disps, bool savePngs = false) {

        //this->plot_f (this->a, disps[2]);
        //this->plot_contour (this->c, disps[4], 0.75);
#if 0
        if (savePngs) {
            // a
            stringstream ff1;
            ff1 << this->logpath << "/a_";
            ff1 << std::setw(5) << std::setfill('0') << frameN;
            ff1 << ".png";
            disps[2].saveImage (ff1.str());
            // contours
            stringstream ff3;
            ff3 << this->logpath << "/cntr_";
            ff3 << std::setw(5) << std::setfill('0') << frameN;
            ff3 << ".png";
            disps[4].saveImage (ff3.str());

            frameN++;
        }
#endif
    }

#if 0
    /*!
     * Plot a field. If individual_scaling is true, then each map is
     * normalised individually, otherwise the group is normalised so
     * that maps can be compared alongside each other.
     */
    void plot_f (vector<vector<double> >& f, morph::Gdisplay& disp, bool individual_scaling=false) {
        vector<double> fix(3, 0.0);
        vector<double> eye(3, 0.0);
        eye[2] = -0.4;
        vector<double> rot(3, 0.0);

        vector<vector<double> > norm_a;
        this->resize_vector_vector (norm_a);
        if (individual_scaling) {
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
            for (unsigned int i = 0; i<this->N; ++i) {
                for (unsigned int h=0; h<this->nhex; h++) {
                    norm_a[i][h] = fmin (fmax (((f[i][h]) - mina[i]) * scalea[i], 0.0), 1.0);
                }
            }
        } else {
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
            for (unsigned int i = 0; i<this->N; ++i) {
                #pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; h++) {
                    norm_a[i][h] = fmin (fmax (((f[i][h]) - mina) * scalea, 0.0), 1.0);
                }
            }
        }

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
#endif
    //@} // Plotting code

}; // RD_OrientPref

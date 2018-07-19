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
     * n(x,t) variable from the Karb2004 paper.
     */
    vector<double> n;

    /*!
     * J_i(x,t) variables - the "flux current of axonal branches of
     * type i"
     */
    vector<vector<double> > J;

    /*!
     * Our choice of dt.
     */
    double dt = 0.001;

    /*!
     * The power to which a_i(x,t) is raised in Eqs 1 and 2 in the
     * paper.
     */
    double k = 3.0;

    /*!
     * The diffusion parameter.
     */
    array<double,2> D = { { 1.0, 1.0 } };

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
     * Simple constructor; no arguments
     */
    RD_2D_Karb (void) {
        this->init();
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
    void init (void) {

        cout << "init() called" << endl;

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
        this->resize_vector_vector (this->J);
        this->resize_vector_vector (this->betaterm);

        this->resize_vector_variable (this->n);
        this->resize_vector_variable (this->rhoA);
        this->resize_vector_variable (this->rhoB);
        this->resize_vector_variable (this->rhoC);

        this->resize_vector_param (this->alpha);
        this->resize_vector_param (this->beta);
        this->resize_vector_param (this->gammaA);
        this->resize_vector_param (this->gammaB);
        this->resize_vector_param (this->gammaC);

        this->resize_gradient_field (this->grad_rhoA);
        this->resize_gradient_field (this->grad_rhoB);
        this->resize_gradient_field (this->grad_rhoC);

        // Resize grad_a
        this->grad_a.resize (this->N);
        for (unsigned int i = 0; i<this->N; ++i) {
            this->resize_gradient_field (this->grad_a[i]);
        }

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

        // Put some hand-crafted 2D Gaussians in rhoA, rhoB and rhoC.
        //    createGaussian ( x      y    gain  sigma    result)
        this->createGaussian (0.2f, 0.05f, 1.0f, 0.2f, this->rhoA);
        this->createGaussian (-0.2f, 0.05f, 1.0f, 0.3f, this->rhoB);
        this->createGaussian (0.0f, -0.05f, 0.8f, 0.25f, this->rhoC);

        // Compute gradients of guidance molecule concentrations once only
        this->spacegrad2D (this->rhoA, this->grad_rhoA);
        this->spacegrad2D (this->rhoB, this->grad_rhoB);
        this->spacegrad2D (this->rhoC, this->grad_rhoC);
#if 0
        for (unsigned int i=0; i<this->N; ++i) {
            this->Gsum[i].resize (nhex, 0.0);
            for (unsigned int h=0; h<nhex; ++h) {
                // Here I'm adding, but does Ji have components? Revisit this question.
                this->Gsum[i][h] = this->gammaA[i] * (this->grad_rhoA[0][h]  + this->grad_rhoA[1][h])
                    + this->gammaB[i] * (this->grad_rhoB[0][h] + this->grad_rhoB[1][h])
                    + this->gammaC[i] * (this->grad_rhoC[0][h] + this->grad_rhoC[1][h]);
            }
        }
#endif
    }

    //! 2D spatial integration of the function f. Result placed in gradf.
    void spacegrad2D (vector<double> f, array<vector<double>, 2>& gradf) {

        // For each Hex, work out the gradient in x and y directions
        // using whatever neighbours can contribute to an estimate.

        // Note - East is positive x; North is positive y.
        for (auto h : this->hg->hexen) {

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

    void step (void) {
        // Do some work. The order of proceedings is:
        // 1. Compute Karb2004 Eq 3.
        // 2. Do integration of a (RK in the 1D model). Involves computing axon branching flux.
        // 3. Do integration of c
    }

    /*!
     * Plot the system on @a disp
     */
    void plot (morph::Gdisplay& disp) {

        // Copies data to plot out of the model
        vector<double> plt = this->rhoA;
        double maxV = -1e7;
        double minV = +1e7;
        // Determines min and max
        for (auto h : this->hg->hexen) {
            if (h.onBoundary() == false) {
                if (plt[h.vi]>maxV) { maxV = plt[h.vi]; }
                if (plt[h.vi]<minV) { minV = plt[h.vi]; }
            }
        }
        double scaleV = 1.0 / (maxV-minV);

        // Determine a colour from min, max and current value
        vector<double> P(this->nhex, 0.0);
        for (unsigned int h=0; h<this->nhex; h++) {
            P[h] = fmin (fmax (((plt[h]) - minV) * scaleV, 0.0), 1.0);
        }

        // Step through vectors or iterate through list? The latter should be just fine here.
        for (auto h : this->hg->hexen) {
            array<float,3> cl = morph::Tools::getJetColorF (P[h.vi]);
            // Reduce the boundary artificially, so it shows up.
            array<float,3> cl_bound = morph::Tools::getJetColorF (P[h.vi] * 0.8);
            if (h.boundaryHex == false) {
                disp.drawHex (h.position(), (h.d/2.0f), cl);
            } else {
                disp.drawHex (h.position(), (h.d/2.0f), cl_bound);
            }
        }
        disp.redrawDisplay();

        return;
    }

    //! Does: f = (alpha * ci) + betaterm. c.f. Karb2004, Eq 1
    vector<double> compute_dci_dt (unsigned int i) { //vector<double> ci, double alpha, vector<double> betaterm) {
        vector<double> dci_dt;
        for (unsigned int xi=0; xi<this->nhex; xi++) {
            dci_dt.push_back (this->alpha[i] * this->c[i][xi] + this->betaterm[i][xi]);
        }
        // Probably best to use another intermediate member variable for dci_dt.
        return dci_dt;
    }

#if 0
    /*!
     * Computes the "flux of axonal branches" term, and adds "addon",
     * which is alpha * ci - betaterm. Computes dJdX, then multiplies by
     * 2.dx as the sum is carried out over two distance increments.
     */
    vector<double> compute_axonalbranchflux (vector<double> ai,
                                             vector<double> gi,
                                             vector<double> addon) {
        double dJdX;
        vector<double> output(this->nhex, 0.0);
        // This all needs a rewrite to work with adjacent hexes, or to
        // work on a single hex and all its neighbours or something
        // like this.
        for (int j=0; j<this->nhex; j++) {
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
            output[j] = dJdX*2.*this->dx + addon[j];
        }
        return output;
    }
#endif

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

#if 0
/*!
 * c is threshold. p is value of the thing you're contouring; P is a
 * vector of the values in the adjacent hexes, with P[0] being the ___
 * hex.
 */
vector<int> getEdges (double p, vector<double>  P, double c)
{
    vector<int> k;
    if (p<c) {
        if (P[0]>c) { k.push_back(0); }
        if (P[1]>c) { k.push_back(1); }
        if (P[2]>c) { k.push_back(2); }
        if (P[3]>c) { k.push_back(3); }
        if (P[4]>c) { k.push_back(4); }
        if (P[5]>c) { k.push_back(5); }
    }
    return k;
}
#endif

int main (int argc, char **argv)
{
    srand(atoi(argv[3]));

    // DISPLAYS
    vector<morph::Gdisplay> displays;
    vector<double> fix(3, 0.0);
    vector<double> eye(3, 0.0);
    eye[2] = -0.4;
    vector<double> rot(3, 0.0);
    displays.push_back (morph::Gdisplay (600, "morphologica", 0.0, 0.0, 0.0));
    displays[0].resetDisplay (fix, eye, rot);
    displays[0].redrawDisplay();

    // Instatiate the model object
    RD_2D_Karb M;
    try {
        M.init();
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
                M.plot (displays[0]);
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

        case 7:
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 7=DISA" << endl;
            displays[0].resetDisplay (fix, eye, rot);
#if 0
            for (int I=0; I<M.nFields; I++) {

                vector<double> contourcol = morph::Tools::getJetColor((double)I / (double)M.nFields);

                vector<double> plt = M.NN[I];
                double maxV = -1e7;
                double minV = +1e7;
                for(int i=0;i<M.nHexes;i++){
                    if(M.C[i]==6){
                        if(plt[i]>maxV){maxV = plt[i];}
                        if(plt[i]<minV){minV = plt[i];}
                    }
                }
                double scaleV = 1./(maxV-minV);
                vector<double> P(M.nHexes,0.);
                for (int i=0; i<M.nHexes; i++) {
                    P[i] = fmin(fmax( ((plt[i])-minV)*scaleV,0.),1.);
                    // M.X[i][2] = P[i];
                }

                double dh = fabs(M.H[0][0]-M.H[0][1])*0.5;

                for(int i=0; i<M.nHexes; i++) {

                    vector <double> cl = morph::Tools::getJetColor(P[i]);

                    vector<double> p(6,0);
                    p[0] = P[M.N[i][0]];
                    p[1] = P[M.N[i][1]];
                    p[2] = P[M.N[i][2]];
                    p[3] = P[M.N[i][3]];
                    p[4] = P[M.N[i][4]];
                    p[5] = P[M.N[i][5]];
                    vector<int> k1 = getEdges(P[i],p,0.5);


                    //k=1;
                    vector <double> cl2(3,0);
                    if(I==0){
                        cl2[0] = 1.;
                    }

                    for(unsigned int j=0;j<k1.size();j++){
                        displays[0].drawHexSeg(M.H[0][i],M.H[1][i],0.01,dh,contourcol[0],contourcol[1],contourcol[2],k1[j]);
                    }

                }
            }
#endif
            displays[0].redrawDisplay();

            break;
        }

        case 8:
        {
            //W.logfile << W.processName << "@" << TIMEcs << ": 8=DISB" << endl;
            displays[0].resetDisplay(fix,eye,rot);
#if 0
            double dh = fabs(M.H[0][0]-M.H[0][1])*0.5;
            for (int i=0; i<M.nHexes; i++) {
                double maxV = -1e7;
                int maxI = 0;
                double sel = 0.;
                for(int I=0;I<M.nFields;I++){
                    sel += M.NN[I][i];
                    if(M.NN[I][i]>maxV){
                        maxV = M.NN[I][i];
                        maxI = I;
                    }
                }
                sel = maxV/sel;
                vector<double> cl = morph::Tools::getJetColor((double)maxI/(double)M.nFields);
                cl[0] *= sel;
                cl[1] *= sel;
                cl[2] *= sel;

                displays[0].drawHex(M.H[0][i],M.H[1][i],0.,dh,cl[0],cl[1],cl[2]);

            }
#endif
            displays[0].redrawDisplay();
            break;
        }

        }
    }

    return 0;
};

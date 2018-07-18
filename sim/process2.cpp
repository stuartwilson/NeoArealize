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
     * Default constructor.
     */
    RD_2D_Karb (void) {
        this->init();
    }

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

    void resize_gradient_field (array<vector<double>, 2>& g) {
        g[0].resize (this->nhex, 0.0);
        g[1].resize (this->nhex, 0.0);
    }

    /*!
     * Initialise HexGrid, variables and parameters.
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

        // Do some work.

#ifdef EXAMPLE_STUFF
        // Enforce no-flux boundary conditions
        vector<double> CCp = CC[I];
        vector<double> NNp = NN[I];
        for (int i=0; i<nHexes; i++) {
            if (C[i]<6) { // then on boundary
                double nfC = 0.0, nfN = 0.0;
                for (int j=0; j<6; j++) {
                    if (N[i][j] != i) {
                        nfC += CCp[N[i][j]];
                        nfN += NNp[N[i][j]];
                    }
                }
                // Sets value of
                CC[I][i] = nfC / (double)C[i]; // value in C[i] used here for normalization
                NN[I][i] = nfN / (double)C[i];
            }
        }
#endif
    }

    /*!
     * Plot the system on @a disp
     */
    string plot (morph::Gdisplay& disp) {
        string rtnMsg("");
#if 0
        // Copies data to plot out of the model
        vector<double> plt = M.NN[0];
        double maxV = -1e7;
        double minV = +1e7;
        // Determines min and max
        for (int i=0; i<M.nHexes; i++) {
            if (M.C[i] == 6) {
                if (plt[i]>maxV) { maxV = plt[i]; }
                if (plt[i]<minV) { minV = plt[i]; }
            }
        }
        double scaleV = 1.0 / (maxV-minV);

        // Determine a colour from min, max and current value
        vector<double> P(M.nhex, 0.0);
        for (int i=0; i<M.nhex; i++) {
            P[i] = fmin (fmax (((plt[i]) - minV) * scaleV, 0.0), 1.0);
            // M.X[i][2] = P[i];
        }
#endif

        // Step through vectors or iterate through list? The latter should be just fine here.
        for (auto h : this->hg->hexen) {
            array<float,3> cl = morph::Tools::getJetColorF (/*P[i]*/ 0.1);
            array<float,3> cl2 = morph::Tools::getJetColorF (/*P[i]*/ 0.3);
            array<float,3> cl3 = morph::Tools::getJetColorF (/*P[i]*/ 0.8);
            array<float,3> cl4 = morph::Tools::getJetColorF (/*P[i]*/ 0.85);
            array<float,3> cl5 = morph::Tools::getJetColorF (/*P[i]*/ 0.95);
            if (h.has_ne && h.has_nse) {
                if (h.boundaryHex == true) {
                    if (h.x == 0.0f && h.y == 0.0f) {
                        disp.drawTriFill (h.position(), h.ne->position(), h.nse->position(), cl3);
                    } else if (h.ne->x == 0.0f && h.ne->y == 0.0f) {
                        disp.drawTriFill (h.position(), h.ne->position(), h.nse->position(), cl4);
                    } else if (h.nse->x == 0.0f && h.nse->y == 0.0f) {
                        disp.drawTriFill (h.position(), h.ne->position(), h.nse->position(), cl5);
                    } else {
                        disp.drawTriFill (h.position(), h.ne->position(), h.nse->position(), cl2);
                    }
                } else {
                    disp.drawTriFill (h.position(), h.ne->position(), h.nse->position(), cl);
                }
            }
            if (h.has_nw && h.has_nnw) {
                disp.drawTriFill (h.position(), h.nw->position(), h.nnw->position(), cl3);
            }
        }
        disp.redrawDisplay();

        return rtnMsg;
    }

    //! Spatial integration
    vector<double> space (vector<double> x, double dx) {
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

    /*!
     * Computes the "flux of axonal branches" term, and adds "addon",
     * which is alpha * ci - betaterm. Computes dJdX, then multiplies by
     * 2.dx as the sum is carried out over two distance increments.
     */
#if 0
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
        vector <string> command;
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
                W.logfile << M.plot (displays[0]); // would be the ideal API
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

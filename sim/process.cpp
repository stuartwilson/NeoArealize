#include "morph/world.h"
#include "morph/sockserve.h"
#include "morph/display.h"
#include "morph/tools.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

#define CIRCLE_RAD 0.75

class ReactDiff
{
public:
    /*!
     * This is the "scale of the map" but might be better named the
     * resolution. If scale is 7 you get 12481 hexes as an upper limit
     * (drawing a circle or other boundary may reduce this). If scale
     * is 6, you get 3169 hexes to work with. If scale is 8, then it's
     * 49537.
     *
     * This determines the number of hexagonal rings of hexagons in
     * the field.
     */
    int scale = 7;

    /*!
     * I also need an actual scaling factor, which says how far in
     * real measurements the disance from one hex to another
     * represents.
     */

    /*!
     * 0 or 1. Determines how grids at different resolutions are
     * overlaid, with respect to their centres and vertices.
     */
    int offset = 0;

    /*!
     * Number of hexes
     */
    int nHexes = 0;

#if 0
    /*!
     * Appears to be unused. populated as 1/maxR.
     */
    double DR;
#endif

    int nFields = 1; // Number of competing reaction-diffusion systems

    /*!
     * Populated with a subset of the data in H.
     */
    vector<vector<double> > X;

    /*!
     * HexData? hex-grid info
     */
    vector<vector<double> > H;

    /*!
     * Neighbourhood information for each hex element.
     */
    vector<vector<double> > N;

    /*!
     * C for "Count" of number of neighbours. if C[i]<6 then you're on
     * a boundary.
     */
    vector<int> C;

    /*!
     * Contains an N for each of the nFields fields.
     */
    vector<vector<double> > NN;

    /*!
     * Contains a C for each of the nFields fields.
     */
    vector<vector<double> > CC;

    ReactDiff (int scale, int offset, double z, int numFields) {

        cout << "ReactDiff constructor called with scale: " << scale
             << ", offset: " << offset << " and numFields: " << numFields << endl;

        this->scale = scale;
        this->offset = offset;
        this->nFields = numFields;

        H.resize(6);

        // This loop builds hexagonal rings of hexagons. S is the
        // number of hexagonal rings.

        int S = pow(2.0, scale-1) - offset;
        double s = S - offset;
        unsigned int uncounted = 0;
        for (int d = 0; d<=S; d++) { // Hex ring number out from centre. Index of central hex is 0.
            for (int r=-S; r<=S; r++) {
                for (int g=-S; g<=S; g++) {
                    for (int b=-S; b<=S; b++) {
                        if ((abs(r) + abs(g) + abs(b) == d * 2)
                            && (r + g + b == 0)) {
                            double x = (g / 2.0 + r) / s;
                            double y = g * (sqrt(3.0) / 2.0) / s;
                            if (x * x + y * y <= CIRCLE_RAD) { // circular
                                //if (x*x/0.1 + y*y/0.8 <= 1) { // ellipse 1
                                //if (x*x/0.8 + y*y/0.1 <= 1) { // ellipse 2
                                //if (x*x/0.8 + y*y/0.01 <= 0.9) { // ellipse 3
                                //if (x*x/0.01 + y*y/0.8 <= 0.9) { // ellipse 4
                                //if (x*x < 0.4 && y*y < 0.4) { // square 1
                                //if (x*x < 0.2 && y*y < 0.2) { // square 2
                                //if (x*x/0.05 - y*y/0.05 <= 0.7) { // hyperbolic 1
                                H[0].push_back(x);                      //X
                                H[1].push_back(y);                      //Y
                                H[2].push_back(r);                      //R
                                H[3].push_back(g);                      //G
                                H[4].push_back(b);                      //B
                                H[5].push_back(d);                      //D
                                nHexes++;
                            } // else outside region
                        } // else do nothing.
                        else {
                            ++uncounted;
                        }
                    }
                }
            }
        }

        cout << "Hexes counted: " << nHexes << " and not counted: " << uncounted << endl;

#if 0
        int maxR = -1;
        for (int i=0; i<nHexes; i++) {
            if (H[2][i] > maxR ) {
                maxR = H[2][i];
            }
        }
        DR = 1./(double)(maxR);
#endif

        // get neighbours
        N.resize (nHexes);
        C.resize (nHexes, 0);
        for (int i=0; i<nHexes; i++) {
            N[i].resize (8, i); // CONNECT ALL TO 'boundary' UNIT AT N+1
            for (int j=0; j<nHexes; j++) {
                int dr = H[2][j] - H[2][i];
                int dg = H[3][j] - H[3][i];
                int db = H[4][j] - H[4][i];
                if (max (max (abs(dr), abs(dg)), abs(db)) == 1) {
                    // anticlockwise from east
                    if (db==-1 && dr==+1) { N[i][0]=j; C[i]++; } // May be out of place (should be 5?)
                    if (db==-1 && dr== 0) { N[i][1]=j; C[i]++; } // East? (should be 0?)
                    if (db== 0 && dr==-1) { N[i][2]=j; C[i]++; }
                    if (db==+1 && dr==-1) { N[i][3]=j; C[i]++; }
                    if (db==+1 && dr== 0) { N[i][4]=j; C[i]++; }
                    if (db== 0 && dr==+1) { N[i][5]=j; C[i]++; }
                }
            }
        }
        X.resize (nHexes);

        // This is a solution to having H, which is 6 long vectors,
        // each nHexes long, and wanting to pass something with
        // Cartesian position information to the drawing code which is
        // nHexes times a vector of 3 coordinates. It may be a hack
        // which can be recoded out...
        for (int i=0; i<nHexes; i++) {
            X[i].resize (3, 0.0);
            X[i][0] = H[0][i];
            X[i][1] = H[1][i];
        }

        NN.resize (nFields);
        CC.resize (nFields);
        for (int i=0; i<nFields; i++) {
            NN[i].resize (nHexes);
            CC[i].resize (nHexes);
        }
    };

    vector<double> getLaplacian (vector<double> Q, double dx) {
        double overdxSquare = 1.0/(dx*dx);
        vector<double> L(nHexes, 0.0);
        for (int i=0; i<nHexes; i++) {
            L[i] = (Q[N[i][0]] + Q[N[i][1]] + Q[N[i][2]] + Q[N[i][3]] + Q[N[i][4]] + Q[N[i][5]] - 6.0 * Q[i]) * overdxSquare;
        }
        return L;
    }

#define FUDGEFAC 50
    void step (double dt, vector<double> Dn, vector<double> Dc) {

        // ds is a grid-resolution independent length scale. If we
        // make a grid with a known length scale, this will be set by
        // choice.
        double ds = FUDGEFAC * fabs (H[0][0] - H[0][1]) / CIRCLE_RAD;

        double beta = 5.0;
        double b = 1.0, mu = 1.0;
        vector<double> chi = Dn;

        // coupling of fields
        vector<double> a(nHexes,0.0);
        for (int j=0; j<nHexes; j++) {
            for (int i=0; i<nFields; i++) {
                a[j] += CC[i][j];
            }
            a[j] = 1.0 - a[j];
        }

        for (int I=0; I<nFields; I++) {
            vector<double> lapN = getLaplacian (NN[I], ds);
            vector<double> lapC = getLaplacian (CC[I], ds);

            // 1. https://pdfs.semanticscholar.org/a103/b0ab83e2553bca7db069e4962049c4f3e966.pdf
            // 2. http://systems-sciences.uni-graz.at/etextbook/sw3/continuousfield.html
            // 3. http://textbooks.opensuny.org/introduction-to-the-modeling-and-analysis-of-complex-systems/

            // G stands for gradient, but G is really gradient of
            // divergence term in eq 1 in the Ermentrout paper.
            vector<double> G(nHexes, 0.0);

            for (int i=0; i<nHexes; i++) {
                // METHOD 3
                double drN = NN[I][N[i][0]] - NN[I][N[i][3]];
                double dgN = NN[I][N[i][1]] - NN[I][N[i][4]];
                double dbN = NN[I][N[i][2]] - NN[I][N[i][5]];

                double drC = CC[I][N[i][0]] - CC[I][N[i][3]];
                double dgC = CC[I][N[i][1]] - CC[I][N[i][4]];
                double dbC = CC[I][N[i][2]] - CC[I][N[i][5]];

                // Possibly one ds too many below (ds^2 not ds^3)
                G[i] = (drN * drC + dgN * dgC + dbN * dbC) / (6.0 * ds * ds * ds) + NN[I][i] * lapC[i];
            }

            // step N
            for (int i=0; i<nHexes; i++) {
                NN[I][i] += dt * (a[i] - b * NN[I][i] + Dn[I] * lapN[i] - chi[I] * G[i]);
            }

            // step C
            double N2;
            for (int i=0; i<nHexes; i++) {
                N2 = NN[I][i] * NN[I][i];
                CC[I][i] += dt * (beta * N2 / (1.0 + N2) - mu * CC[I][i] + Dc[I] * lapC[i]);
            }

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
        }
    }
}; // ReactDiff

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

int main (int argc, char **argv)
{
    srand(atoi(argv[3]));

    // integration timestep; should be 0.001 to match Ermentrout 2009.
    double dt = .0001;

    // INITIALIZATION
    morph::World W(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), dt);

    // DISPLAYS
    vector<morph::Gdisplay> displays;
    vector<double> fix(3, 0.0);
    vector<double> eye(3, 0.0);
    eye[2] = -0.4;
    vector<double>rot(3, 0.0);
    displays.push_back (morph::Gdisplay (600, "morphologica", 0.0, 0.0, 0.0));
    displays[0].resetDisplay (fix, eye, rot);
    displays[0].redrawDisplay();

    ReactDiff M(7, 0, 0.0, 5);
    for (int I=0; I<M.nFields; I++) {
        for (int i=0; i<M.nHexes; i++) {
            M.NN[I][i] = (morph::Tools::randFloat()) * 0.1;
            M.CC[I][i] = (morph::Tools::randFloat()) * 0.1;
        }
    }

    unsigned int frameN = 0;
    vector<double> Dn (M.nFields, 100.0);
    vector<double> Dc (M.nFields, 100.0 * 0.3);
    double TIME=0.;
    vector <double*> f;
    f.push_back(&TIME);

    vector<vector<double*> > S;
    {
        vector<double*> s;
        for (unsigned int I=0; I<M.CC.size(); I++) {
            for (unsigned int i=0; i<M.CC[I].size(); i++) {
                s.push_back (&M.CC[I][i]);
            }
            S.push_back (s);
        }
    }
    {
        vector<double*> s;
        for (unsigned int I=0; I<M.NN.size(); I++) {
            for (unsigned int i=0; i<M.NN[I].size(); i++) {
                s.push_back (&M.NN[I][i]);
            }
            S.push_back (s);
        }
    }

    bool doing = true;
    while (doing) {

        std::stringstream TIMEss;
        TIMEss << setw(10) << setfill('0') << TIME;
        const char* TIMEcs = TIMEss.str().c_str();

        std::stringstream out;
        out.clear();
        out.setf (ios::fixed, ios::floatfield);

        // DEFINE OUTPUT MESSAGE
        out << TIME << ",";

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
            //W.logfile << W.processName << "@" << TIMEcs << ": 1=STEP" << endl;

            // Exchange comms
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

            // DO STUFF HERE

            M.step (dt, Dn, Dc);
            // END STEP
            TIME++;
            break;
        }

        case 2: // *** PLOT ***
        {
            //W.logfile << W.processName << "@" << TIMEcs << ": 8=DISA" << endl;
            displays[0].resetDisplay (fix, eye, rot);

            vector<double> plt = M.NN[0];
            double maxV = -1e7;
            double minV = +1e7;
            for (int i=0; i<M.nHexes; i++) {
                if (M.C[i] == 6) {
                    if (plt[i]>maxV) { maxV = plt[i]; }
                    if (plt[i]<minV) { minV = plt[i]; }
                }
            }
            double scaleV = 1.0 / (maxV-minV);
            vector<double> P(M.nHexes, 0.0);
            for (int i=0; i<M.nHexes; i++) {
                P[i] = fmin (fmax (((plt[i]) - minV) * scaleV, 0.0), 1.0);
                // M.X[i][2] = P[i];
            }

            for (int i=0; i<M.nHexes; i++) {
                vector<double> cl = morph::Tools::getJetColor (P[i]);
                displays[0].drawTriFill (M.X[i], M.X[M.N[i][0]], M.X[M.N[i][1]], cl);
                displays[0].drawTriFill (M.X[i], M.X[M.N[i][3]], M.X[M.N[i][4]], cl);
            }
            displays[0].redrawDisplay();
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

            switch(stoi(command[1])){
            case 0:
            {
                Dn[stoi(command[2])] = stod(command[3]);
                break;
            }
            case 1:
            {
                Dc[stoi(command[2])] = stod(command[3]);
                break;
            }
            }
            break;
        }

        case 5: // *** SAVE ***
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 5=SAVE" << endl;

            if (command.size() == 2) {
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
            } else {
                W.logfile << "No output filename." << endl;
            }
            break;
        }

        case 6: // *** LOAD ***
        {
            W.logfile << W.processName << "@" << TIMEcs << ": 6=LOAD" << endl;

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
            break;
        }

        case 7:
        {
            //W.logfile << W.processName << "@" << TIMEcs << ": 8=DISA" << endl;
            displays[0].resetDisplay (fix, eye, rot);

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
                    //displays[0].drawHex(M.H[0][i],M.H[1][i],0.,dh,cl[0],cl[1],cl[2]);
                    //displays[0].drawHex(M.H[0][i],M.H[1][i],0.,dh,cl[0]*0.8,cl[1]*0.8,cl[2]*0.8);

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
            /*
              for(int i=0;i<M.n;i++){
              vector <double> cl = tools::getJetColor(P[i]);
              displays[0].drawTriFill(M.X[i],M.X[M.N[i][0]],M.X[M.N[i][1]],cl);
              displays[0].drawTriFill(M.X[i],M.X[M.N[i][3]],M.X[M.N[i][4]],cl);
              }
            */
            displays[0].redrawDisplay();

            break;
        }

        case 8:
        {
            //W.logfile << W.processName << "@" << TIMEcs << ": 8=DISA" << endl;

            double dh = fabs(M.H[0][0]-M.H[0][1])*0.5;

            displays[0].resetDisplay(fix,eye,rot);
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

            displays[0].redrawDisplay();
            break;
        }

        }
    }

    return 0;
};

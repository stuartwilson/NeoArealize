#include "morph/world.h"
#include "morph/sockserve.h"
#include "morph/display.h"
#include "morph/tools.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

class Erm2009
{
public:
    int scale;
    int offset;
    int n;
    double DR;

    vector<vector<double> > X; //
    vector<vector<double> > H; // hex-grid info
    vector<vector<double> > N; // hex neighbourhood
    vector<int> C;

    vector<double> NN, CC;
    Erm2009 (int scale, int offset, double z) {

        this->scale = scale;
        this->offset = offset;

        H.resize(6);
        int S = pow(2.0, scale-1) - offset;
        double s = S - offset;
        n = 0;
        for (int d=0; d<=S; d++) {
            for (int r=-S; r<=S; r++) {
                for (int g=-S; g<=S; g++) {
                    for (int b=-S;b<=S;b++) {
                        if (abs(r)+abs(g)+abs(b)==d*2 && (r+g+b==0)) {
                            double x = (g/2.+r)/s;
                            double y = g*(sqrt(3.)/2.)/s;
                            if (x*x+y*y <= 0.75) { // circular
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
                                n++;
                            }
                        }
                    }
                }
            }
        }

        int maxR = -1;
        for(int i=0;i<n;i++){
            if(H[2][i]>maxR){
                maxR = H[2][i];
            }
        }
        DR = 1./(double)(maxR);

        // get neighbours
        N.resize(n);
        C.resize(n,0);
        for(int i=0;i<n;i++){
            N[i].resize(6,i); // CONNECT ALL TO 'boundary' UNIT AT N+1
            for(int j=0;j<n;j++){
                int dr = H[2][j]-H[2][i];
                int dg = H[3][j]-H[3][i];
                int db = H[4][j]-H[4][i];
                if(max(max(abs(dr),abs(dg)),abs(db))==1){
                    //anticlockwise from east
                    if(db==-1&&dr==+1){N[i][0]=j;C[i]++;}
                    if(db==-1&&dr== 0){N[i][1]=j;C[i]++;}
                    if(db== 0&&dr==-1){N[i][2]=j;C[i]++;}
                    if(db==+1&&dr==-1){N[i][3]=j;C[i]++;}
                    if(db==+1&&dr== 0){N[i][4]=j;C[i]++;}
                    if(db== 0&&dr==+1){N[i][5]=j;C[i]++;}
                }
            }
        }
        X.resize(n);
        for(int i=0;i<n;i++){
            X[i].resize(3,0.);
            X[i][0] = H[0][i];
            X[i][1] = H[1][i];
        }
        NN.resize(n);
        CC.resize(n);
    };

    vector<double> getLaplacian(vector<double> Q, double dx) {
        double overdxSquare = 1./(dx*dx);
        vector<double> L(n,0.);
        for(int i=0;i<n;i++){
            L[i]=(Q[N[i][0]]+Q[N[i][1]]+Q[N[i][2]]+Q[N[i][3]]+Q[N[i][4]]+Q[N[i][5]]-6.*Q[i])*overdxSquare;
        }
        return L;
    }

    void step(double dt, double Dn, double Dc) {

        // ORIGINAL MASK (a,b,c) OR dr, dg, bd METHOD
        // 100,0.3* coffeebean
        // 80, 0.3* bullseye
        // 60,0.3*  mercedes
        // 40,0.3*  baseball
        // 10, 0.3* hexagons (10 counted)
        // 5, 0.3* hexagons (17 counted)

        // double Dn = 80.;
        // double Dc = 0.3*Dn;

        double ds = 50.*fabs(H[0][0]-H[0][1]) / 0.75; // SHOULD BE 25.*DR
#ifdef METHOD_1_OR_2
        double over4dsSquared = 1./(4.*ds*ds);
#endif
        double beta = 5.;
        double a = 1., b = 1., mu = 1., chi = Dn;

        vector<double> lapN = getLaplacian(NN,ds);
        vector<double> lapC = getLaplacian(CC,ds);

        //
        // 1. https://pdfs.semanticscholar.org/a103/b0ab83e2553bca7db069e4962049c4f3e966.pdf
        // 2. http://systems-sciences.uni-graz.at/etextbook/sw3/continuousfield.html
        // 3. http://textbooks.opensuny.org/introduction-to-the-modeling-and-analysis-of-complex-systems/
#ifdef METHOD_1_OR_2
        double dxN, dyN, dxC, dyC;
        double ma= 0.22061, mb=0.43174, mc= 0.37665; // mask for gradient on hex (ref 1 above) (ref 2/3 above)
#endif
        vector<double> G(n,0.);

        for (int i=0;i<n;i++) {
#ifdef METHOD_1
            dxN = mb*NN[N[i][0]]+ma*NN[N[i][1]]-ma*NN[N[i][2]]-mb*NN[N[i][3]]-ma*NN[N[i][4]]+ma*NN[N[i][5]];
            dyN = mc*NN[N[i][1]]+mc*NN[N[i][2]]-mc*NN[N[i][4]]-mc*NN[N[i][5]];
            dxC = mb*CC[N[i][0]]+ma*CC[N[i][1]]-ma*CC[N[i][2]]-mb*CC[N[i][3]]-ma*CC[N[i][4]]+ma*CC[N[i][5]];
            dyC = mc*CC[N[i][1]]+mc*CC[N[i][2]]-mc*CC[N[i][4]]-mc*CC[N[i][5]];
            G[i] = (dxN*dxC+dyN*dyC)*over4dsSquared + NN[i]*lapC[i];
#endif
#ifdef METHOD_2
            double h1N = NN[N[i][1]]+NN[N[i][2]]-NN[N[i][4]]-NN[N[i][5]];
            double h2N = NN[N[i][0]]+NN[N[i][1]]-NN[N[i][3]]-NN[N[i][4]];
            double h3N = NN[N[i][2]]+NN[N[i][3]]-NN[N[i][5]]-NN[N[i][0]];
            double h1C = CC[N[i][1]]+CC[N[i][2]]-CC[N[i][4]]-CC[N[i][5]];
            double h2C = CC[N[i][0]]+CC[N[i][1]]-CC[N[i][3]]-CC[N[i][4]];
            double h3C = CC[N[i][2]]+CC[N[i][3]]-CC[N[i][5]]-CC[N[i][0]];
            dxN = h1N;
            dyN = h2N-h3N;
            dxC = h1C;
            dyC = h2C-h3C;
            G[i] = (dxN*dxC+dyN*dyC)*over4dsSquared + NN[i]*lapC[i];
#endif

            // METHOD 3
            double drN = NN[N[i][0]]-NN[N[i][3]];
            double dgN = NN[N[i][1]]-NN[N[i][4]];
            double dbN = NN[N[i][2]]-NN[N[i][5]];

            double drC = CC[N[i][0]]-CC[N[i][3]];
            double dgC = CC[N[i][1]]-CC[N[i][4]];
            double dbC = CC[N[i][2]]-CC[N[i][5]];

            G[i] = (drN*drC+dgN*dgC+dbN*dbC)/(6.*ds*ds*ds) + NN[i]*lapC[i];
        }

        // step N
        for(int i=0;i<n;i++){
            NN[i]+=dt*( a-b*NN[i] + Dn*lapN[i] - chi*G[i]);
        }

        // step C
        double N2;
        for(int i=0;i<n;i++){
            N2 = NN[i]*NN[i];
            CC[i]+=dt*( beta*N2/(1.+N2) - mu*CC[i] + Dc*lapC[i] );
        }

        // Enforce no-flux boundary
        vector<double> CCp = CC;
        vector<double> NNp = NN;
        for(int i=0;i<n;i++){
            if(C[i]<6){
                double nfC = 0., nfN = 0.;
                for(int j=0;j<6;j++){
                    if(N[i][j] !=i){
                        nfC += CCp[N[i][j]];
                        nfN += NNp[N[i][j]];
                    }
                }
                CC[i] = nfC / (double)C[i];
                NN[i] = nfN / (double)C[i];
            }
        }
    }
}; // Erm2009


int main (int argc, char **argv)
{
    srand(atoi(argv[3]));

    double dt = .0001;            // integration timestep // SHOULD BE 0.001

    // INITIALIZATION
    morph::World W(argv[1],argv[2],atoi(argv[3]),atoi(argv[4]),dt);
    //W.waitForConnection();

    // DISPLAYS
    vector<morph::Gdisplay> displays;
    vector<double>fix(3,0.);
    vector<double>eye(3,0.);
    eye[2] = -0.4;
    vector<double>rot(3,0.);
    displays.push_back (morph::Gdisplay (600, "morphologica", 0., 0., 0.));
    displays[0].resetDisplay(fix,eye,rot);
    displays[0].redrawDisplay();

    Erm2009 M(7,0,0.);
    for (int i=0;i<M.n;i++) {
        M.NN[i]=(morph::Tools::randFloat())*0.1;
        M.CC[i]=(morph::Tools::randFloat())*0.1;
    }

    unsigned int frameN = 0;
    double Dn = 100.;
    double Dc = Dn*0.3;
    double TIME=0.;
    vector <double*> f;
    f.push_back(&TIME);

    vector<vector<double*> > S;
    {
        vector<double*> s;
        for(unsigned int i=0;i<M.CC.size();i++){
            s.push_back(&M.CC[i]);
        } S.push_back(s);
    }
    {
        vector<double*> s;
        for(unsigned int i=0;i<M.NN.size();i++){
            s.push_back(&M.NN[i]);
        } S.push_back(s);
    }

    bool doing = true;
    while (doing) {

        std::stringstream TIMEss;
        TIMEss<<setw(10)<<setfill('0')<<TIME;
        const char* TIMEcs = TIMEss.str().c_str();

        std::stringstream out;
        out.clear();
        out.setf(ios::fixed,ios::floatfield);

        // DEFINE OUTPUT MESSAGE
        out<<TIME<<",";

        vector <string> command;
        string messageI=W.master.exchange(out.str().c_str());
        stringstream ss(messageI);
        while (ss.good()){
            string substr;
            getline(ss,substr,',');
            command.push_back(substr);
        } ss.clear();

        // Interpret commands:
        switch (stoi(command[0])) {

        case 0: // *** QUIT ***
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 0=QUIT"<<endl<<flush;

            W.logfile.close();
            for(unsigned int i=0;i<displays.size();i++){
                displays[i].closeDisplay();
            }
            W.master.closeSocket();
            for(unsigned int i=0;i<W.ports.size();i++){
                W.ports[i].closeSocket();
            }
            doing=false;
            break;
        }

        case 1: // *** STEP ***
        {
            //W.logfile<<W.processName<<"@"<<TIMEcs<<": 1=STEP"<<endl<<flush;

            // Exchange comms
            vector< vector <string> > commands(W.ports.size());
            for(unsigned int i=0;i<W.ports.size();i++){
                string messageI=W.ports[i].exchange(out.str().c_str());
                stringstream ss(messageI);
                while (ss.good()){
                    string substr;
                    getline(ss,substr,',');
                    commands[i].push_back(substr);
                } ss.clear();
            }

            // DO STUFF HERE

            M.step(dt, Dn, Dc);
            // END STEP
            TIME++;
            break;
        }

        case 2: // *** PLOT ***
        {
            //W.logfile<<W.processName<<"@"<<TIMEcs<<": 8=DISA"<<endl<<flush;
            displays[0].resetDisplay(fix,eye,rot);

            vector<double> plt = M.NN;
            double maxV = -1e7;
            double minV = +1e7;
            for (int i=0;i<M.n;i++) {
                if (M.C[i]==6) {
                    if (plt[i]>maxV) { maxV = plt[i]; }
                    if (plt[i]<minV) { minV = plt[i]; }
                }
            }
            double scaleV = 1./(maxV-minV);
            vector<double> P(M.n,0.);
            for (int i=0;i<M.n;i++) {
                P[i] = fmin (fmax (((plt[i])-minV)*scaleV,0.),1.);
                // M.X[i][2] = P[i];
            }

            for (int i=0;i<M.n;i++) {
                vector <double> cl = morph::Tools::getJetColor(P[i]);
                displays[0].drawTriFill(M.X[i],M.X[M.N[i][0]],M.X[M.N[i][1]],cl);
                displays[0].drawTriFill(M.X[i],M.X[M.N[i][3]],M.X[M.N[i][4]],cl);
            }
            displays[0].redrawDisplay();
            break;
        }

        case 3:
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 3=RECD"<<endl<<flush;
            std::stringstream frameFile1;
            frameFile1<<"logs/"<<W.processName<<"frameA";
            frameFile1<<setw(5)<<setfill('0')<<frameN;
            frameFile1<<".png";
            displays[0].saveImage(frameFile1.str());
            frameN++;
            break;
        }

        case 4:
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 4=RECD"<<endl<<flush;
            switch(stoi(command[1]))
            {
            case 0: Dn = stod(command[2]); break;
            case 1: Dc = stod(command[2]); break;
            }
            break;
        }

        case 5: // *** SAVE ***
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 5=SAVE"<<endl<<flush;

            if (command.size()==2) {
                ofstream outFile;
                std::stringstream oFile; oFile<<command[1];
                outFile.open(oFile.str().c_str(),ios::out|ios::binary);
                double n=(double)S.size();double* N=&n;
                outFile.write((char*)N,sizeof(double));
                for (unsigned int i=0;i<S.size();i++) {
                    double m=(double)S[i].size();double* M=&m;
                    outFile.write((char*)M,sizeof(double));
                    for (unsigned int j=0;j<S[i].size();j++) {
                        outFile.write((char*)S[i][j],sizeof(double));
                    }
                }
                outFile.close();
            } else {
                W.logfile<<"No output filename."<<endl<<flush;
            }
            break;
        }

        case 6: // *** LOAD ***
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 6=LOAD"<<endl<<flush;

            if (command.size()==2) {
                double dummy;
                ifstream inFile;
                std::stringstream iFile;
                iFile<<command[1];
                inFile.open(iFile.str().c_str(),ios::in|ios::binary);
                inFile.read((char*)&dummy,sizeof(double));
                int I=(int)dummy;
                if (I==static_cast<int>(S.size())) {
                    for (int i=0;i<I;i++) {
                        inFile.read((char*)&dummy,sizeof(double));
                        int J=(int)dummy;
                        if (J==static_cast<int>(S[i].size())) {
                            for (int j=0;j<J;j++) {
                                inFile.read((char*)S[i][j],sizeof(double));
                            }
                        } else {
                            W.logfile<<"Wrong dims I."<<endl<<flush;
                        }
                    }
                } else {
                    W.logfile<<"Wrong dims J."<<endl<<flush;
                }
                inFile.close();
            } else {
                W.logfile<<"No input filename."<<endl<<flush;
            }
            break;
        }
        }
    }

    return 0;
};

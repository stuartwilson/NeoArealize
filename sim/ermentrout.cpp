#include "rd_2d_erm.h"

#include "morph/display.h"
#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main (int argc, char **argv){

    if (argc < 2) {
        cerr << "\nUsage: ./build/sim/process w0\n\n";
        cerr << "Be sure to run from the base source directory.\n";
        return -1;
    }
    vector<morph::Gdisplay> displays;
    vector<double> fix(3, 0.0);
    vector<double> eye(3, 0.0);
    vector<double> rot(3, 0.0);

    double rhoInit = 1.4;
    string worldName(argv[1]);
    string winTitle = worldName + ": n";
    displays.push_back (morph::Gdisplay (500, 500, 100, 0, winTitle.c_str(), rhoInit, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    rhoInit = 1.4;
    winTitle = worldName + ": c";
    displays.push_back (morph::Gdisplay (500, 500, 100, 0, winTitle.c_str(), rhoInit, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    // Instantiate the model object
    RD_2D_Erm M;
    try {
        M.init (displays);
    } catch (const exception& e) {
        cerr << "Exception initialising RD_2D_Karb object: " << e.what() << endl;
    }

    M.Dn = stod(argv[3]);
    M.chi = M.Dn;
    M.Dc = 0.3*M.Dn;

    // Start the loop
    bool doing = true;
    while (doing) {
        // Step the model
        try {
            M.step();
        } catch (const exception& e) {
            cerr << "Caught exception calling M.step(): " << e.what() << endl;
            doing = false;
        }

        displays[0].resetDisplay (fix, eye, rot);
        try {
            M.plot (displays);
        } catch (const exception& e) {
            cerr << "Caught exception calling M.plot(): " << e.what() << endl;
            doing = false;
        }
    }

    return 0;
};

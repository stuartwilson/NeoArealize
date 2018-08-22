#include "rd_orientpref.h"

#include "morph/display.h"
#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main (int argc, char **argv)
{
    if (argc < 1) {
        cerr << "\nUsage: ./build/sim/orient";
        cerr << "Be sure to run from the base NeoArealize source directory.\n";
        return -1;
    }

    // Set RNG seed
    int rseed = 1;
    srand(rseed);

    cout << "First rand num for seed " << rseed << " is: " << morph::Tools::randDouble() << endl;

    // Instantiate the model object
    RD_OrientPref RD;

    // Do any modifications to RD parameters here.
    // RD.myparam = 3.4;

    // Call the init function, which can allocate variables and run
    // any pre-stepping computations.
    try {
        RD.init();
    } catch (const exception& e) {
        cerr << "Exception initialising RD object: " << e.what() << endl;
    }

    // How many iterations to compute?
    unsigned int T = 100; // Care - while debugging h5 files are being created, with T=10000 30 GB of logs will be created!

    // Start the loop
    bool finished = false;
    while (!finished) {
        // Step the model
        try {
            RD.step();
        } catch (const exception& e) {
            cerr << "Caught exception calling RD.step(): " << e.what() << endl;
            finished = true;
        }
        if (RD.stepCount > T) {
            finished = true;
        }
    }

    // Save final state of RD system
    RD.saveState();

    return 0;
};

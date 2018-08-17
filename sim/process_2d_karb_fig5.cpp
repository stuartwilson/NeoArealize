/*
 * Reproduce approximate behaviour of Fig 5 of the Karbowski 2004
 * paper in our 2D system of arbitrary boundary shape.
 */

#include "rd_2d_karb.h"

#include "morph/display.h"
#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main (void)
{
    // Set RNG seed
    int rseed = 1;
    srand(rseed);

    // Create some displays
    vector<morph::Gdisplay> displays;
    vector<double> fix(3, 0.0);
    vector<double> eye(3, 0.0);
    eye[2] = -0.4;
    vector<double> rot(3, 0.0);

    double rhoInit = 1.5;
    string worldName("2dkarb_f5");
    string winTitle = worldName + ": emx_pax_fgf";
    displays.push_back (morph::Gdisplay (1020, 300, 100, 0, winTitle.c_str(), rhoInit, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    winTitle = worldName + ": rhoA_rhoB_rhoC";
    displays.push_back (morph::Gdisplay (1020, 300, 100, 300, winTitle.c_str(), rhoInit, 0.0, 0.0, displays[0].win));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    winTitle = worldName + ": a[0] to a[4]";
    displays.push_back (morph::Gdisplay (1700, 300, 100, 600, winTitle.c_str(), rhoInit, 0.0, 0.0, displays[0].win));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    winTitle = worldName + ": c[0] to c[4]";
    displays.push_back (morph::Gdisplay (1700, 300, 100, 900, winTitle.c_str(), rhoInit, 0.0, 0.0, displays[0].win));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    // SW - Contours
    winTitle = worldName + ": contours";
    displays.push_back (morph::Gdisplay (500, 500, 100, 900, winTitle.c_str(), rhoInit, 0.0, 0.0, displays[0].win));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    // Instantiate the model object
    RD_2D_Karb M;
    try {
        M.init (displays);
    } catch (const exception& e) {
        cerr << "Exception initialising RD_2D_Karb object: " << e.what() << endl;
    }

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
            // Save some frames ('c' variable only for now)
            // FIXME: Think about this.
            if (M.stepCount % 100 == 0) {
                M.saveC();
            }

            if (M.stepCount > 6000) {
                doing = false;
            }

        } catch (const exception& e) {
            cerr << "Caught exception calling M.plot(): " << e.what() << endl;
            doing = false;
        }
    }

    return 0;
};

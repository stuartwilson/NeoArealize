/*!
 * Implementation of HexGrid
 */
#include "HexGrid.h"
#include <cmath>
#include <iostream>
#include <sstream>

using std::ceil;
using std::abs;
using std::cout;
using std::endl;
using std::stringstream;

morph::HexGrid::HexGrid (float d, float width, float height, float fixz)
{
    this->hextohex = d;
    this->x_span = width;
    this->y_span = height;
    this->z = fixz;

    this->init();
}

unsigned int
morph::HexGrid::num (void) const
{
    return this->hexen.size();
}

string
morph::HexGrid::output (void) const
{
    stringstream ss;
    ss << "Hex grid with " << this->hexen.size() << " hexes." << endl;
    auto i = this->hexen.begin();
    float lasty = this->hexen.front().y;
    unsigned int rownum = 0;
    ss << endl << "Row " << rownum++ << ":" << endl;
    while (i != this->hexen.end()) {
        if (i->y > lasty) {
            ss << endl << "Row " << rownum++ << ":" << endl;
            lasty = i->y;
        }
        ss << i->output() << endl;
        ++i;
    }
    return ss.str();
}

void
morph::HexGrid::init (void)
{
    // First get the number of g indices that we need to go down from the centre of the HexGrid:
    float halfY = this->y_span/2.0f;
    float dv = (this->hextohex*morph::SQRT_OF_3_F)/2.0f;
    int distG = ceil(abs(halfY/dv));
    //cout << "distG: " << distG << " g-indices" << endl;

    // Now add up the distance we've gone in the horizontal direction.
    float x_travel = distG * this->hextohex / 2.0f;
    float x_rest = (this->x_span/2.0f) - x_travel;
    //cout << "x_travel: " << x_travel << " cart. units." << endl;
    //cout << "x_rest: " << x_rest  << " cart. units."<< endl;

    // Work out how many r-indices this is to the left and to the right for the first row.
    int distR_L = ceil(abs(x_rest/this->hextohex));
    int distR_R = ceil(abs( ((x_span/2.0f)+x_travel)/this->hextohex ));
    //cout << "distR_L: " << distR_L << endl;
    //cout << "distR_R: " << distR_R << endl;

    /*
     * First run through and create the grid.
     */
    cout << "Creating hex grid..." << endl;

    // The "vector iterator" - this is an identity iterator that is
    // added to each Hex in the grid. If hexes are removed from the
    // grid, then these may need to be updated in the remaining hexes.
    unsigned int vi = 0;
    unsigned int rowCount = 0;
    for (int gi = -distG; gi<=distG; ++gi) {
        // Add the hexes for the row.
        for (int ri = -distR_L; ri <= distR_R; ++ri) {
            //cout << "gi:" << gi << " ri:" << ri << " ";
            Hex h(vi++, this->hextohex, ri, gi);
            this->hexen.push_back (h);
        }
        if (rowCount++ % 2 == 0) {
            distR_L++;
            distR_R--;
        }
    }

    /*
     * Second run through; re-define the grid according to boundary.
     */
    cout << "Applying boundary..." << endl;

    /*
     * Final run through. Re-compute a suitable vector iterator for
     * the list of Hexes. Populate neighbour pointers for each Hex.
     */

}

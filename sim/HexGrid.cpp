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
morph::HexGrid::num (void)
{
    return this->hexen.size();
}

string
morph::HexGrid::output (void)
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
morph::HexGrid::initcart (void)
{
    // Number of hexes in the first row. The second row will have
    // rowN+2 hexes; the third row rowN hexes and so on, alternating.
    unsigned int rowN = ceil(this->x_span/this->hextohex);
    unsigned int colN = ceil(this->y_span/this->hextohex);

    /*
     * First run through and create the grid.
     */

    // How much we will move the centre of the hexes as we advance up
    // each row of hexes:
    float yhextoyhex = 2.0f*this->hextohex/morph::SQRT_OF_3_F;

    // The "vector iterator" - this is an identity iterator that is
    // added to each Hex in the grid. If hexes are removed from the
    // grid, then these may need to be updated in the remaining hexes.
    unsigned int vi = 0;

    for (unsigned int yi = 0; yi<colN; ++yi) {
        // Holds the number of hexes in the current row (this alternates).
        unsigned int _rowN = 0;
        // The running Cartesian x, y of the Hexes being added to the HexGrid.
        float xpos, ypos;
        if (yi%2 == 0) { // even row
            _rowN = rowN;
            // On even row, reset x position to 0
            xpos = 0.0f;
        } else { // odd row
            _rowN = rowN+1;
            // On even row, reset x position to a little left of 0.
            xpos = 0.0f-this->hextohex/2.0f;
        }
        // Add the hexes for the row.
        for (unsigned int xi = 0; xi < _rowN; ++xi) {
            Hex h(vi++, this->hextohex, xpos, ypos);
            this->hexen.push_back (h);
            xpos += this->hextohex;
        }
        ypos += yhextoyhex;
    }

    /*
     * Second run through; re-define the grid according to boundary.
     */
}

void
morph::HexGrid::init (void)
{
    // Number of hexes in the first row. The second row will have
    // rowN+2 hexes; the third row rowN hexes and so on, alternating.
    //int halfRowN = ceil(abs(this->x_span/2*this->hextohex));

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
}

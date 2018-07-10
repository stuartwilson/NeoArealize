/*!
 * Implementation of HexGrid
 */
#include "HexGrid.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <morph/BezCurvePath.h>
#include <morph/BezCoord.h>

using std::ceil;
using std::abs;
using std::cout;
using std::endl;
using std::stringstream;
using std::vector;
using std::tuple;
using std::make_tuple;
using std::runtime_error;

using morph::BezCurvePath;
using morph::BezCoord;
using morph::Hex;

morph::HexGrid::HexGrid (float d_, float x_span_, float y_span_, float z_)
{
    this->d = d_;
    this->x_span = x_span_;
    this->y_span = y_span_;
    this->z = z_;

    this->init();
}

void
morph::HexGrid::setBoundary (const BezCurvePath& p)
{
    this->boundary = p;

    /*
     * Second run through; re-define the grid according to boundary.
     */
    if (!this->boundary.isNull()) {
        cout << "Applying boundary..." << endl;

        // Compute the points on the boundary using the hex to hex
        // spacing as the step size.
        vector<BezCoord> bpoints = this->boundary.getPoints (this->d);

        auto bpi = bpoints.begin();
        list<Hex>::iterator lhi = this->hexen.begin();
        while (bpi != bpoints.end()) {
            lhi = this->setBoundary (*bpi++, lhi);
        }

        // Now something like:
        // this->discardOutside();
    }


    /*
     * Final run through. Re-compute a suitable vector iterator for
     * the list of Hexes. Populate neighbour pointers for each Hex.
     */
    this->setupNeighbours();
}

list<Hex>::iterator
morph::HexGrid::setBoundary (const BezCoord& point, list<Hex>::iterator startFrom)
{
    auto i = startFrom;
    while (i != this->hexen.end()) {
        ++i;
    }
    return startFrom;
}

void
morph::HexGrid::setupNeighbours (void)
{
    auto i = this->hexen.begin();

    // Matchers for neighbours on all sides.
#if 0
    tuple<int,int,int> neighE = make_tuple(0, 0, 0);
    tuple<int,int,int> neighNE = neighE;
    tuple<int,int,int> neighNW = neighE;
    tuple<int,int,int> neighW = neighE;
    tuple<int,int,int> neighSW = neighE;
    tuple<int,int,int> neighSE = neighE;
#endif

    vector<int> neighbours;
    neighbours.resize (18, 0);

    while (i != this->hexen.end()) {

        if (i->bi != 0) {
            throw runtime_error ("HexGrid::setupNeighbours assumes bi is not used and always set to 0!");
        }

        // NB: Assume bi always 0 and hence no degenerate solutions exist for neigh* tuples.
#if 0
        neighE  = make_tuple (i->ri+1, i->gi,   i->bi);
        neighNE = make_tuple (i->ri,   i->gi+1, i->bi);
        neighNW = make_tuple (i->ri-1, i->gi+1, i->bi);
        neighW  = make_tuple (i->ri-1, i->gi,   i->bi);
        neighSW = make_tuple (i->ri,   i->gi-1, i->bi);
        neighSE = make_tuple (i->ri+1, i->gi-1, i->bi);
#endif

        neighbours[0] = i->ri+1;
        neighbours[1] = i->gi;
        // neighbours[2] = bi
        neighbours[3] = i->ri;
        neighbours[4] = i->gi+1;
        // neighbours[5] = bi
        neighbours[6] = i->ri-1;
        neighbours[7] = i->gi+1;
        // neighbours[8] = bi
        neighbours[9] = i->ri-1;
        neighbours[10] = i->gi;
        // neighbours[11] = bi
        neighbours[12] = i->ri;
        neighbours[13] = i->gi-1;
        // neighbours[14] = bi
        neighbours[15] = i->ri+1;
        neighbours[16] = i->gi-1;
        // neighbours[17] = bi

        auto j = this->hexen.begin();
        while (j != this->hexen.end()) {
            if (j != i) {
                int neighbourNum = this->checkNeighbour (*j, neighbours);
                switch (neighbourNum) {
                case 1: // East
                    //cout << "East neighbour of (" << i->ri << "," << i->gi << "," << i->bi << ") is (";
                    //cout << j->ri << "," << j->gi << "," << j->bi << ")" << endl;
                    i->ne = &(*j);
                    break;
                case 2: // NE
                    i->nne = &(*j);
                    break;
                case 3: // NW
                    i->nnw = &(*j);
                    break;
                case 4: // West
                    i->nw = &(*j);
                    break;
                case 5: // SW
                    i->nsw = &(*j);
                    break;
                case 6: // SE
                    i->nse = &(*j);
                    break;
                default:
                    // j is not a neighbour of i
                    break;
                }
            }
            ++j;
        }

        ++i;
    }
}

// Opportunity to do this vectorised with SIMD/SSE or similar.
int
morph::HexGrid::checkNeighbour (const Hex& candidate, const vector<int>& neighbourRGB)
{
    int neighbourNum = 0;
    for (unsigned int i = 0; i<18; i+=3) {
        if (neighbourRGB[i] == candidate.ri
            && neighbourRGB[i+1] == candidate.gi
            /*&& neighbourRBG[i+2] == candidate.bi*/) {
            neighbourNum = (i+1)/3;
            break;
        }
    }
    return neighbourNum;
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
morph::HexGrid::initRect (void)
{
    // First get the number of g indices that we need to go down from the centre of the HexGrid:
    float halfY = this->y_span/2.0f;
    float dv = (this->d*morph::SQRT_OF_3_F)/2.0f;
    int distG = ceil(abs(halfY/dv));
    //cout << "distG: " << distG << " g-indices" << endl;

    // Now add up the distance we've gone in the horizontal direction.
    float x_travel = distG * this->d / 2.0f;
    float x_rest = (this->x_span/2.0f) - x_travel;
    //cout << "x_travel: " << x_travel << " cart. units." << endl;
    //cout << "x_rest: " << x_rest  << " cart. units."<< endl;

    // Work out how many r-indices this is to the left and to the right for the first row.
    int distR_L = ceil(abs(x_rest/this->d));
    int distR_R = ceil(abs( ((x_span/2.0f)+x_travel)/this->d ));
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
            Hex h(vi++, this->d, ri, gi);
            this->hexen.push_back (h);
        }
        if (rowCount++ % 2 == 0) {
            distR_L++;
            distR_R--;
        }
    }

    // Now setup neighbour information
    cout << "Set up neighbours..." << endl;
    this->setupNeighbours();


    cout << "init() done" << endl;
}

void
morph::HexGrid::init (void)
{
    // Use span_x to determine how many rings out to traverse.
    float halfX = this->x_span/2.0f;
    unsigned int maxRing = abs(ceil(halfX/this->d));

    cout << "Creating hexagonal hex grid..." << endl;

    // The "vector iterator" - this is an identity iterator that is
    // added to each Hex in the grid.
    unsigned int vi = 0;

    // Vectors of pointers to hexes in this->hexen. Used to keep a
    // track of nearest neighbours. I'm using vector, rather than a
    // list as this allows fast random access of elements and I'll not
    // be inserting or erasing in the middle of the arrays.
    vector<Hex*> prevRingEven;
    vector<Hex*> prevRingOdd;

    // Swap pointers between rings.
    vector<Hex*>* prevRing = &prevRingEven;
    vector<Hex*>* nextPrevRing = &prevRingOdd;

    // Direction iterators used in the loop.
    int ri = 0;
    int gi = 0;

    // Create central "ring" first (the single hex)
    this->hexen.emplace (vi++, this->d, ri, gi);

    // Put central ring in the prevRing vector:
    prevRing->push_back (&this->hexen.back());

    // Now build up the rings around it, setting neighbours as we
    // go. Each ring has 6 more hexes than the previous one (except
    // for ring 1, which has 6 instead of 1 in the centre).
    unsigned int numInRing = 6;

    // How many hops in the same direction before turning a corner?
    // Increases for each ring. Increases by 1 in each ring.
    unsigned int ringSideLen = 1;

    for (unsigned int ring = 0; ring <= maxRing; ++ring) {

        // Set start ri, gi. This moves up a hex and left a hex.
        --ri;
        ++gi;

        // Now walk around the ring:
        //
        // Walk in r direction
        for (unsigned int i = 0; i<ringSideLen; ++i) {
            this->hexen.emplace (vi++, this->d, ri++, gi);
            Hex& h = this->hexen.back();
            // Set my neighbour(s):
            h.nse = prevRing[0];
            // Set me as neighbour to those in prevRing:
            prevRing[0]->nnw = &h;

            // Put in me nextPrevRing:
            nextPrevRing->push_back (&h);
        }
        // Walk in -b direction
        for (unsigned int i = 0; i<ringSideLen; ++i) {
            this->hexen.emplace (vi++, this->d, ri++, gi--);

        }
        // Walk in -g direction
        for (unsigned int i = 0; i<ringSideLen; ++i) {
            this->hexen.emplace (vi++, this->d, ri, gi--);

        }
        // Walk in -r direction
        for (unsigned int i = 0; i<ringSideLen; ++i) {
            this->hexen.emplace (vi++, this->d, ri--, gi);

        }
        // Walk in b direction
        for (unsigned int i = 0; i<ringSideLen; ++i) {
            this->hexen.emplace (vi++, this->d, ri--, gi++);

        }
        // Should now be on the last hex.

        // Always 6 additional hexes in the next ring out
        numInRing += 6;
        // And ring side length goes up by 1
        ringSideLen++;

        // Swap prevRing and nextPrevRing.
        vector<Hex*>* tmp = prevRing;
        prevRing = nextPrevRing;
        nextPrevRing = tmp;
    }

    cout << "init() done" << endl;
}

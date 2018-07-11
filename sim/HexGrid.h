#ifndef _HEXGRID_H_
#define _HEXGRID_H_

#include "Hex.h"
#include <morph/BezCurvePath.h>

#include <list>
#include <string>

using std::list;
using std::string;
using morph::Hex;

namespace morph {

    class HexGrid
    {
    public:
        /*!
         * Construct the hex grid with a hex to hex distance of
         * d (centre to centre) and approximate spans of width
         * and height. Set z to fixz.
         */
        HexGrid (float d_, float x_span_, float y_span_, float z_ = 0.0f);

        /*!
         * Sets boundry to p, then runs the code to discard hexes
         * lying outside this boundary. Finishes up by calling
         * setupNeighbours.
         */
        void setBoundary (const BezCurvePath& p);

        /*!
         * Find the closest Hex in this->hexen to the coordinate
         * point, and set its onBoundary attribute to true.
         */
        list<Hex>::iterator setBoundary (const BezCoord& point, list<Hex>::iterator startFrom);

        /*!
         * Return the number of hexes in the grid.
         */
        unsigned int num (void) const;

        /*!
         * Output some text information about the hexgrid.
         */
        string output (void) const;

        /*!
         * The list of hexes.
         */
        list<Hex> hexen;

    private:
        /*!
         * Initialise a grid of hexes in a hex spiral, setting
         * neighbours as we go in, hopefully, a fairly effecient
         * algorithm.
         */
        void init (void);

        /*!
         * Check to see if the Hex candidate has the same ri, gi, bi
         * as the six candidates held in the vector of ints
         * neighbourRGB. The return code is 0: candidate does not
         * match any of neighbourRGB and is not a neighbour; 1:
         * candidate is an East neighbour; 2: candidate is a NE
         * neighbour; 3: candidate is a NW neighbour; 4: candidate is
         * a West neighbour; 5: candidate is a SW neighbour; 6:
         * candidate is a SE neighbour.
         */
        int checkNeighbour (const Hex& candidate, const vector<int>& neighbourRGB);

        /*!
         * Centre to centre hex distance.
         */
        float d = 1.0f;

        /*!
          * Make hexes in the horizonal direction over a span of
         * x_span.
         */
        float x_span = 10.0f;

        /*!
         * Make hexes in the vertical direction over a span of y_span.
         */
        float y_span = 10.0f;

        /*!
         * The z coordinate of this hex grid layer
         */
        float z;

        /*!
         * A boundary to apply to the initial, rectangular grid.
         */
        BezCurvePath boundary;
    };

} // namespace morph

#endif // _HEXGRID_H_

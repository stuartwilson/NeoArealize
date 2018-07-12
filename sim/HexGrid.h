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
         * Construct the hexagonal hex grid with a hex to hex distance
         * of d (centre to centre) and approximate diameter of
         * x_span_. Set z to z_.
         */
        HexGrid (float d_, float x_span_, float z_ = 0.0f);

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
         * Show the coordinates of the vertices of the overall hex
         * grid generated.
         */
        string extent (void) const;

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
         * Recursively mark hexes to be kept if they are inside the boundary.
         */
        void markHexesInside (list<Hex>::iterator hi);

        /*!
         * Discard hexes in this->hexen that are outside the boundary.
         */
        void discardOutside (void);

        /*!
         * Set true when a new boundary has been applied. This means
         * that the vertexE, vertexW, etc pointers are no longer
         * valid.
         */
        bool gridReducedToBoundary = false;

#ifdef DEPRECATED
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
#endif

        /*!
         * Centre to centre hex distance.
         */
        float d = 1.0f;

        /*!
         * Give the hexagonal hex grid a diameter of approximately
         * x_span.
         */
        float x_span = 10.0f;

        /*!
         * The z coordinate of this hex grid layer
         */
        float z;

        /*!
         * A boundary to apply to the initial, rectangular grid.
         */
        BezCurvePath boundary;

        /*!
         * Hex references to the hexes on the vertices of the
         * hexagonal grid. Configured during init().
         */
        //@{
        list<Hex>::iterator vertexE;
        list<Hex>::iterator vertexNE;
        list<Hex>::iterator vertexNW;
        list<Hex>::iterator vertexW;
        list<Hex>::iterator vertexSW;
        list<Hex>::iterator vertexSE;
        //@}
    };

} // namespace morph

#endif // _HEXGRID_H_

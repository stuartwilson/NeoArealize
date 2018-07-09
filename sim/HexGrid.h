#ifndef _HEXGRID_H_
#define _HEXGRID_H_

#include "Hex.h"

#include <list>
#include <string>

using std::list;
using std::string;

namespace morph {

    class HexGrid
    {
    public:
        /*!
         * Construct the hex grid with a hex to hex distance of
         * d (centre to centre) and approximate spans of width
         * and height. Set z to fixz.
         */
        HexGrid (float d, float width, float height, float fixz = 0.0f);

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
         * Initialise the Hex elements.
         */
        void init (void);

        /*!
         * Centre to centre hex distance.
         */
        float hextohex = 1.0f;

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
    };

} // namespace morph

#endif // _HEXGRID_H_

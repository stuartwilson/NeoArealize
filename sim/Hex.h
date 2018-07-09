#ifndef _HEX_H_
#define _HEX_H_

#include <tuple>
#include <string>

using std::tuple;
using std::get;
using std::string;
using std::to_string;

//#define DEBUG_DECONSTRUCT 1
#ifdef DEBUG_DECONSTRUCT
#include <iostream>
using std::cout;
using std::endl;
#endif

namespace morph {

    const float SQRT_OF_3_F = 1.73205081;
    /*!
     * Describes a regular hexagon arranged with vertices pointing
     * vertically and two flat sides perpendicular to the horizontal
     * axis:
     *
     *            *
     *         *     *
     *         *     *
     *            *
     *
     * The centre of the hex in a Cartesian right hand coordinate
     * system is represented with x, y and z:
     *
     *  y
     *  ^
     *  |
     *  |
     *  0-----> x     z out of screen/page
     *
     * Directions are "r" "g" and "b" and their negatives:
     *
     *         b  * g
     * -r <--  *     * ---> r
     *         *     *
     *         -g * -b
     *
     */
    class Hex
    {
    public:
        /*!
         * Constructor - will need to determine the best argument list for this.
         */
        Hex (const unsigned int& idx, const float& d_,
             const float& x_, const float& y_) {
            this->vi = idx;
            this->d = d_;
            this->x = x_;
            this->y = y_;
        }

        Hex (const unsigned int& idx, const float& d_,
             const int& r_, const int& g_) {
            this->vi = idx;
            this->d = d_;
            this->ri = r_;
            this->gi = g_;
            this->computeCartesian();
        }

        string output (void) {
            string s("Hex ");
            s += to_string(this->vi) + " x:";
            s += to_string(this->x) + " y:";
            s += to_string(this->y) + ". ";
            return s;
        }

        /*!
         * Convert ri, gi and bi indices into x and y coordinates
         * based on the hex-to-hex distance d.
         */
        void computeCartesian (void) {
            this->x = this->d*this->ri + (d/2.0f)*this->gi - (d/2.0f)*this->bi;
            float dv = (this->d*morph::SQRT_OF_3_F)/2.0f;
            this->y = dv*this->gi + dv*this->bi;
            //cout << "x:" << x << " y:" << y << endl;
        }

        /*!
         * Vector index. This is the index into those data vectors
         * which hold the relevant data pertaining to this hex. This
         * is a scheme which allows me to keep the data in separate
         * vectors and all the hex position information in this class.
         * What happens when I delete some hex elements?  Simple - I
         * can re-set the vi indices after creating a grid of Hex
         * elements and then pruning down.
         */
        unsigned int vi;

        /*!
         * Cartesian coordinates of the centre of the Hex.
         */
        //@{
        float x = 0.0f;
        float y = 0.0f;
        float z = 0.0f;
        //@}

        /*!
         * The centre-to-centre distance from one Hex to an
         * immediately adjacent Hex.
         */
        float d = 1.0f;

        /*!
         * The distance from the centre of the Hex to any of the
         * vertices.
         */
        float getRv (void) {
            float rv = this->d/morph::SQRT_OF_3_F;
            return rv;
        }

        /*!
         * The vertical distance between hex centres on adjacent rows.
         */
        float getDv (void) {
            float dv = (this->d*morph::SQRT_OF_3_F)/2.0f;
            return dv;
        }

        /*!
         * Indices in hex directions. These lie in the x-y
         * plane. They index in positive and negative directions,
         * starting from the Hex at (0,0,z)
         */
        //@{
        /*!
         * Index in r direction - positive "East", that is in the +x
         * direction.
         */
        int ri = 0;
        /*!
         * Index in g direction - positive "NorthEast". In a direction
         * 30 degrees East of North or 60 degrees North of East.
         */
        int gi = 0;
        /*!
         * Index in b direction - positive "NorthEast". In a direction
         * 30 degrees East of North or 60 degrees North of East.
         */
        int bi = 0;
        //@}

        /*!
         * Nearest neighbours
         */
        //@{
        /*!
         * Nearest neighbour to the East; in the plus r direction.
         */
        Hex* nr = (Hex*)0;
        /*!
         * Nearest neighbour to the NorthEast; in the plus g
         * direction.
         */
        Hex* ng = (Hex*)0;
        /*!
         * Nearest neighbour to the NorthWest; in the plus b
         * direction.
         */
        Hex* nb = (Hex*)0;
        /*!
         * Nearest neighbour to the West; in the minus r direction.
         */
        Hex* nmr = (Hex*)0;
        /*!
         * Nearest neighbour to the SoutWest; in the minus g
         * direction.
         */
        Hex* nmg = (Hex*)0;
        /*!
         * Nearest neighbour to the SouthEast; in the minus b
         * direction.
         */
        Hex* nmb = (Hex*)0;
        //@}
    };

} // namespace morph

#endif // _HEX_H_

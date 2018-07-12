#include "HexGrid.h"
#include <iostream>
#include <morph/ReadCurves.h>

using namespace morph;
using namespace std;

int main()
{
    int rtn = 0;
    try {
        ReadCurves r("../trial.svg");

        HexGrid hg(0.05, 3);
        hg.setBoundary (r.getCorticalPath());

        cout << hg.extent() << endl;
        //cout << hg.output() << endl;

    } catch (const exception& e) {
        cerr << "Caught exception reading trial.svg: " << e.what() << endl;
        rtn = -1;
    }
    return rtn;
}

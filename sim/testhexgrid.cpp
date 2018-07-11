#include "HexGrid.h"
#include <iostream>

using namespace morph;
using namespace std;

int main()
{
    HexGrid hg(2.309401, 7, 6);
    cout << hg.output() << endl;
    cout << "Set up " << hg.num() << " hexes in a grid." << endl;
    return 0;
}

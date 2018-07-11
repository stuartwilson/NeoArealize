#include "HexGrid.h"
#include <iostream>

using namespace morph;
using namespace std;

int main()
{
    HexGrid hg(1, 6);
    cout << hg.output() << endl;
    cout << "Set up " << hg.num() << " hexes in a grid." << endl;
    return 0;
}

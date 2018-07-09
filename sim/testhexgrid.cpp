#include "HexGrid.h"
#include <iostream>

using namespace morph;
using namespace std;

int main()
{
    HexGrid hg(2, 12, 6);
    cout << hg.output() << endl;
    return 0;
}

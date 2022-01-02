#include <iostream>
#include <parasolid.h>

#include <vector>
#include <map>
#include <part.h>


int main(int argc, char** argv) {
    PartOptions options;
    auto part = Part(TEST_PART, options);

    return 0;
}

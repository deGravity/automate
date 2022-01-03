#include <iostream>
#include <parasolid.h>

#include <vector>
#include <map>
#include <part.h>


int main(int argc, char** argv) {
    PartOptions options;
    options.onshape_style = false;
    options.default_mcfs_only_face_axes = false;

    auto part = Part(TEST_PART, options);

    return 0;
}

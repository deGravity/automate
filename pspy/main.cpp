#include <iostream>
#include <parasolid.h>

#include <vector>
#include <map>
#include <part.h>


int main(int argc, char** argv) {
    PartOptions options;
    options.onshape_style = false;
    options.default_mcfs_only_face_axes = false;
    options.num_uv_samples = 100;
    options.collect_inferences = false;
    options.default_mcfs = false;

    auto part = Part(TEST_PART, options);

    return 0;
}

#include <iostream>
#include <parasolid.h>

#include <vector>
#include <map>
#include <part.h>


int main(int argc, char** argv) {
    PartOptions options;
    options.onshape_style = false;
    options.default_mcfs_only_face_axes = false;

    auto failed_part0 = Part("C:/Users/ben/data/models/dfd74ebbece1644bca68cd7b/8e741c3ba52c3840d00ccdf0/83238564274046cede121261/kfydfb2i/jjiui.x_t", options);
    auto failed_part1 = Part("C:/Users/ben/data/models/0e2510289c094498d83e1212/93bb0b6ef5c3a379ed4dfb28/aa1012fdd1fd69d52bb11167/default/jjdei.x_t", options);
    auto failed_part2 = Part("C:/Users/ben/data/models/fc98db0d70734fb15dd6f9f2/bb20bbbc6177b278e2e775d4/040b8c7f84fd101a6809cd66/default/jndeiqq.x_t", options);
    auto failed_part3 = Part("C:/Users/ben/data/models/edb3a5a93bc79766e77d40e8/dc4fefcbd8ce119875755cc9/1d76fd01693d9ee9cbbdd0ad/default/jjdei.x_t", options);
    auto failed_part4 = Part("C:/Users/ben/data/models/57385bc2e4b06c68b35e411a/eba206fd89154785816b4302/b2f2847ed483aca02d514222/default/jndheoa.x_t", options);
    auto failed_part5 = Part("C:/Users/ben/data/models/5a965d854700c9832920381a/3b272aa7f9ba91cea7373cef/4e294aae95d5f3d62c521550/default/jjgui.x_t", options);
    auto failed_part6 = Part("C:/Users/ben/data/models/bafddbcdbd60c5bddd92eec1/03346a8e7880805f2f0f2bca/e06904b567fd265896916128/default/jndgusi.x_t", options);
    auto failed_part7 = Part("C:/Users/ben/data/models/589dfa6b4f6a920fd166cc74/be912a8b37b2f5f32c91880c/070eca70fd1babc3f1613289/default/jndeyqy.x_t", options);
    auto failed_part8 = Part("C:/Users/ben/data/models/821a27e0f6a7a25a467ea008/b507d8bc4c2004d61412943a/da39610cddf74769ed2abafa/default/kjeeera.x_t", options);
    auto failed_part9 = Part("C:/Users/ben/data/models/56eccd07e4b0673b5de76f3e/a2b345ad9eee6d52868bee3b/ad1ef5e091da0cab79cf05d0/default/jjdei.x_t", options);


    //auto part = Part(TEST_PART, options);

    return 0;
}

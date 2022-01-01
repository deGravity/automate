#ifndef VERTEX_H_INCLUDED
#define VERTEX_H_INCLUDED 1

#include <Eigen/Core>
#include <vector>
#include "types.h"

struct Vertex {
    Vertex(int id);

    int _id;
    Eigen::Vector3d position;

    std::vector<Inference> get_inferences();

};

#endif // !VERTEX_H_INCLUDED

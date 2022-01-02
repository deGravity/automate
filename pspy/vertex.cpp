#include "vertex.h"
#include "vertex.h"

#include <parasolid.h>

Vertex::Vertex(int id) {
    _id = id;

    PK_ERROR_t err = PK_ERROR_no_errors;

    PK_POINT_t point;
    err = PK_VERTEX_ask_point(_id, &point);
    assert(err == PK_ERROR_no_errors); // PK_VERTEX_ask_point
    PK_POINT_sf_t point_sf;
    PK_POINT_ask(point, &point_sf);
    position = Eigen::Vector3d(point_sf.position.coord[0], point_sf.position.coord[1], point_sf.position.coord[2]);
}

std::vector<Inference> Vertex::get_inferences()
{
    std::vector<Inference> inferences;

    Inference inf;
    inf.inference_type = InferenceType::POINT;
    inf.origin = position;
    inf.z_axis = Eigen::Vector3d(0.0, 0.0, 1.0);

    inferences.push_back(inf);

    return inferences;
}

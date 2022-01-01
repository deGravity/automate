#ifndef LOOP_H_INCLUDED
#define LOOP_H_INCLUDED 1

#include <parasolid.h>
#include "types.h"
#include <assert.h>
#include <Eigen/Core>
#include <vector>

struct Loop {
    Loop(int id);

    std::vector<Inference> get_inferences();

    int _id;
    LoopType _type;

    bool _is_circle;

    double length;
    Eigen::Vector3d center_of_gravity;
    Eigen::MatrixXd moment_of_inertia;

    Eigen::Vector3d na_bb_center;
    Eigen::Vector3d na_bb_x;
    Eigen::Vector3d na_bb_z;
    Eigen::MatrixXd na_bounding_box;
};


#endif // LOOP_H_INCLUDED

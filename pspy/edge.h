#ifndef EDGE_H_INCLUDED
#define EDGE_H_INCLUDED 1

#include "types.h"
#include <Eigen/Core>
#include <vector>

struct Edge {
    Edge(int id);

    void init_line();

    void init_circle();

    void init_ellipse();

    void init_bb();

    void init_nabb();

    std::vector<Inference> get_inferences();

    void add_inferences_circle_or_ellipse(std::vector<Inference>& inferences);
    void add_inferences_line(std::vector<Inference>& inferences);
    void add_inferences_other(std::vector<Inference>& inferences);

    int _id;
    int _curve;

    CurveFunction function;
    std::vector<double> parameters;

    double t_start;
    double t_end;

    Eigen::Vector3d start;
    Eigen::Vector3d end;
    bool has_curve;
    bool _is_reversed; // Is edge opposite curve? Already baked-in to start and
                       // end, but needed for tangent computation.

    bool is_periodic;
    Eigen::Vector3d mid_point;

    double length;
    Eigen::Vector3d center_of_gravity;
    Eigen::MatrixXd moment_of_inertia;

    Eigen::MatrixXd bounding_box;
    Eigen::Vector3d na_bb_center;
    Eigen::Vector3d na_bb_x;
    Eigen::Vector3d na_bb_z;
    Eigen::MatrixXd na_bounding_box;

    void sample_points(
        const int num_points,
        const bool sample_tangents,
        std::vector<Eigen::VectorXd>& samples,
        Eigen::Vector2d& t_range);
};


#endif // !EDGE_H_INCLUDED

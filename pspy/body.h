#ifndef BODY_H_INCLUDED
#define BODY_H_INCLUDED 1

#include <parasolid.h>
#include "face.h"
#include "loop.h"
#include "edge.h"
#include "topology.h"


#include <map>
#include <Eigen/Core>
#include <iostream>
#include "types.h"

class Body {
public:
    Body(int id);
    ~Body();

    BREPTopology GetTopology();

    MassProperties GetMassProperties(double accuracy = MASS_ACC);

    Eigen::MatrixXd GetBoundingBox();

    int Transform(const Eigen::MatrixXd& xfrm);

    void Tesselate(
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        Eigen::VectorXi& FtoT,
        Eigen::MatrixXi& EtoT,
        Eigen::VectorXi& VtoT);

    void debug();

private:
    int _id;
    bool _valid; // If we need to re-compute due to transforms
};

// Helper Functions
bool is_body(int id);
std::vector<Body> read_xt(std::string path);

#endif // !BODY_H_INCLUDED

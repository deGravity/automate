#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED 1

#include <Eigen/Core>
#include <assert.h>
#include <parasolid.h>

enum class TopologyType {
    FACE,
    EDGE,
    VERTEX,
    LOOP
};

enum class SurfaceFunction {
    PLANE,
    CYLINDER,
    CONE,
    SPHERE,
    TORUS,
    SPUN,
    BSURF,
    OFFSET,
    SWEPT,
    BLENDSF,
    MESH,
    FSURF,
    NONE // Faces sometimes have no surface
};

enum class CurveFunction {
    LINE,
    CIRCLE,
    ELLIPSE,
    BCURVE,
    ICURVE,
    FCURVE,
    SPCURVE,
    TRCURVE,
    CPCURVE,
    PLINE,
    NONE // Edges sometimes have no curve
};

enum class LoopType {
    VERTEX,
    WIRE,
    OUTER,
    INNER,
    WINDING,
    INNER_SING,
    LIKELY_OUTER,
    LIKELY_INNER,
    UNCLEAR,
    ERROR
};

enum class InferenceType {
    CENTER, // CIRCLE, ELLIPSE, SPHERE
    CENTROID, // PLANE - ignores inner loops (COG centroid)
    MID_POINT, // EDGES that aren't circles or ellipses
    POINT, // VERTEX, tip of CONE
    TOP_AXIS_POINT, // CYLINDER, CONE, TORUS, SPUN
    BOTTOM_AXIS_POINT, // ...
    MID_AXIS_POINT, // ...
    LOOP_CENTER // First edge of planar loop, also not CIRCLE
};

const double XFRM_TOL = 0.999;
const double MASS_ACC = 0.999;

struct MassProperties {
    MassProperties(int* ids, double accuracy = MASS_ACC, int num_ids = 1) {
        PK_ERROR_t err = PK_ERROR_no_errors;
        PK_TOPOL_eval_mass_props_o_t mass_props_options;
        PK_TOPOL_eval_mass_props_o_m(mass_props_options);

        m_of_i.resize(3, 3);

        err = PK_TOPOL_eval_mass_props(num_ids, ids, accuracy,
            &mass_props_options,
            &amount,
            &mass,
            c_of_g.data(),
            m_of_i.data(),
            &periphery);
        assert(err == PK_ERROR_no_errors || err == PK_ERROR_missing_geom); // PK_TOPOL_eval_mass_props
        if (err == PK_ERROR_missing_geom) {
            amount = 0;
            mass = 0;
            c_of_g.setZero();
            m_of_i = Eigen::MatrixXd::Zero(3, 3);
            periphery = 0;
        }
    }
    double amount;
    double mass;
    Eigen::Vector3d c_of_g;
    Eigen::MatrixXd m_of_i;
    double periphery;
};

struct Inference {
    Eigen::Vector3d z_axis;
    Eigen::Vector3d origin;
    InferenceType inference_type;

    // For interfacing with the Onshape CAD System
    // Onshape has special case rules for some coordinate
    // inferences that can exclude an inference for redundancy
    // (e.g. all circular inner loops), or flip the axis
    // (e.g. circular or elliptical curves oriented opposite
    // of a planar face containing them)
    bool onshape_inference = true; // Onshape allows this inference
    bool flipped_in_onshape = false; // Onshape would flip our inference
};

#endif // !TYPES_H_INCLUDED
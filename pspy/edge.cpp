#include "edge.h"
#include "edge.h"
#include "edge.h"
#include "edge.h"
#include "edge.h"
#include "edge.h"

#include <parasolid.h>
#include <assert.h>
#include <float.h>

Edge::Edge(int id) {
    _id = id;

    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_CURVE_t curve;
    PK_CLASS_t curve_class;
    PK_VECTOR_t ends[2];
    PK_INTERVAL_t t_int;
    PK_LOGICAL_t sense;
    err = PK_EDGE_ask_geometry(_id, PK_LOGICAL_true /*yes interval*/, &curve, &curve_class, ends, &t_int, &sense);
    assert(err == PK_ERROR_no_errors); // PK_EDGE_ask_geometry

    _curve = curve;
    _is_reversed = false;

    if (curve != PK_ENTITY_null) {
        has_curve = true;
        if (sense) { // interval from start to end
            start = Eigen::Vector3d(ends[0].coord[0], ends[0].coord[1], ends[0].coord[2]);
            end = Eigen::Vector3d(ends[1].coord[0], ends[1].coord[1], ends[1].coord[2]);
            t_start = t_int.value[0];
            t_end = t_int.value[1];
        }
        else { // interval from end to start
            end = Eigen::Vector3d(ends[0].coord[0], ends[0].coord[1], ends[0].coord[2]);
            start = Eigen::Vector3d(ends[1].coord[0], ends[1].coord[1], ends[1].coord[2]);
            t_end = t_int.value[0];
            t_start = t_int.value[1];
            _is_reversed = true; // Needed for tangent evaluation
        }
        PK_PARAM_sf_t param;
        err = PK_CURVE_ask_param(curve, &param);
        assert(err == PK_ERROR_no_errors); // PK_CURVE_ask_param

        is_periodic = (param.periodic == PK_PARAM_periodic_no_c) ? false : true;

        PK_VECTOR_t p;
        PK_CURVE_eval(curve, (t_start + t_end) / 2, 0, &p);
        assert(err == PK_ERROR_no_errors);
        mid_point = Eigen::Vector3d(p.coord[0], p.coord[1], p.coord[2]);

        auto m = MassProperties(&_id);
        length = m.amount;
        center_of_gravity = m.c_of_g;
        moment_of_inertia = m.m_of_i;

    }
    else {
        has_curve = false;
        start = Eigen::Vector3d(0, 0, 0);
        end = Eigen::Vector3d(0, 0, 0);
        t_start = 0;
        t_end = 0;
        is_periodic = false;
        mid_point = Eigen::Vector3d(0, 0, 0);

        length = 0;
        center_of_gravity.setZero();
        moment_of_inertia = Eigen::MatrixXd::Zero(3, 3);
    }

    init_bb();

    init_nabb();

    if (!has_curve) {
        function = CurveFunction::NONE;
        return;
    }

    switch (curve_class) {
    case PK_CLASS_line:
        init_line();
        break;
    case PK_CLASS_circle:
        init_circle();
        break;
    case PK_CLASS_ellipse:
        init_ellipse();
        break;
    case PK_CLASS_bcurve:
        function = CurveFunction::BCURVE;
        break;
    case PK_CLASS_icurve:
        function = CurveFunction::ICURVE;
        break;
    case PK_CLASS_fcurve:
        function = CurveFunction::FCURVE;
        break;
    case PK_CLASS_spcurve:
        function = CurveFunction::SPCURVE;
        break;
    case PK_CLASS_trcurve:
        function = CurveFunction::TRCURVE;
        break;
    case PK_CLASS_cpcurve:
        function = CurveFunction::CPCURVE;
        break;
    case PK_CLASS_pline:
        function = CurveFunction::PLINE;
        break;
    }

}

void Edge::init_line() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_LINE_sf_t line_sf;
    err = PK_LINE_ask(_curve, &line_sf);
    assert(err == PK_ERROR_no_errors); // PK_LINE_ask
    function = CurveFunction::LINE;
    parameters.push_back(line_sf.basis_set.location.coord[0]);
    parameters.push_back(line_sf.basis_set.location.coord[1]);
    parameters.push_back(line_sf.basis_set.location.coord[2]);
    parameters.push_back(line_sf.basis_set.axis.coord[0]);
    parameters.push_back(line_sf.basis_set.axis.coord[1]);
    parameters.push_back(line_sf.basis_set.axis.coord[2]);
}

void Edge::init_circle() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_CIRCLE_sf_t circle_sf;
    err = PK_CIRCLE_ask(_curve, &circle_sf);
    assert(err == PK_ERROR_no_errors); // PK_CIRCLE_ask
    function = CurveFunction::CIRCLE;
    parameters.push_back(circle_sf.basis_set.location.coord[0]);
    parameters.push_back(circle_sf.basis_set.location.coord[1]);
    parameters.push_back(circle_sf.basis_set.location.coord[2]);
    parameters.push_back(circle_sf.basis_set.axis.coord[0]);
    parameters.push_back(circle_sf.basis_set.axis.coord[1]);
    parameters.push_back(circle_sf.basis_set.axis.coord[2]);
    parameters.push_back(circle_sf.radius);
}

void Edge::init_ellipse() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_ELLIPSE_sf_t ellipse_sf;
    err = PK_ELLIPSE_ask(_curve, &ellipse_sf);
    assert(err == PK_ERROR_no_errors); // PK_ELLIPSE_ask
    function = CurveFunction::LINE;
    parameters.push_back(ellipse_sf.basis_set.location.coord[0]);
    parameters.push_back(ellipse_sf.basis_set.location.coord[1]);
    parameters.push_back(ellipse_sf.basis_set.location.coord[2]);
    parameters.push_back(ellipse_sf.basis_set.axis.coord[0]);
    parameters.push_back(ellipse_sf.basis_set.axis.coord[1]);
    parameters.push_back(ellipse_sf.basis_set.axis.coord[2]);
    parameters.push_back(ellipse_sf.R1);
    parameters.push_back(ellipse_sf.R2);
}

void Edge::init_bb() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    // Get Bounding Box
    PK_BOX_t box;
    err = PK_TOPOL_find_box(_id, &box);
    assert(err == PK_ERROR_no_errors || err == PK_ERROR_missing_geom);
    if (err == PK_ERROR_missing_geom) {
        bounding_box = Eigen::MatrixXd::Zero(2, 3);
    }
    else {
        bounding_box.resize(2, 3);
        bounding_box <<
            box.coord[0], box.coord[1], box.coord[2],
            box.coord[3], box.coord[3], box.coord[5];
    }
}

void Edge::init_nabb() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    // Get Non-Aligned Bounding Box
    // Axes are ordered so X is largest, Y is second largest, Z is smallest

    if (!has_curve) {
        na_bb_center.setZero();
        na_bb_x.setZero();
        na_bb_z.setZero();
        na_bounding_box = Eigen::MatrixXd::Zero(2, 3);
        return;
    }

    PK_NABOX_sf_t nabox;
    PK_TOPOL_find_nabox_o_t nabox_options;
    PK_TOPOL_find_nabox_o_m(nabox_options);
    err = PK_TOPOL_find_nabox(1, &_id, NULL, &nabox_options, &nabox);
    assert(err == PK_ERROR_no_errors); // PK_TOPOL_find_nabox

    na_bb_center <<
        nabox.basis_set.location.coord[0],
        nabox.basis_set.location.coord[1],
        nabox.basis_set.location.coord[2];
    na_bb_x <<
        nabox.basis_set.ref_direction.coord[0],
        nabox.basis_set.ref_direction.coord[1],
        nabox.basis_set.ref_direction.coord[2];
    na_bb_z <<
        nabox.basis_set.axis.coord[0],
        nabox.basis_set.axis.coord[1],
        nabox.basis_set.axis.coord[2];
    na_bounding_box.resize(2, 3);
    na_bounding_box <<
        nabox.box.coord[0], nabox.box.coord[1], nabox.box.coord[2],
        nabox.box.coord[3], nabox.box.coord[4], nabox.box.coord[5];
}

std::vector<Inference> Edge::get_inferences()
{
    std::vector<Inference> inferences;

    switch (function) {
    case CurveFunction::CIRCLE:
        add_inferences_circle_or_ellipse(inferences);
        break;
    case CurveFunction::ELLIPSE:
        add_inferences_circle_or_ellipse(inferences);
        break;
    case CurveFunction::LINE:
        add_inferences_line(inferences);
        break;
    default:
        add_inferences_other(inferences);
        break;
    }

    return inferences;
}

void Edge::add_inferences_circle_or_ellipse(std::vector<Inference>& inferences)
{
    Inference inf;
    inf.inference_type = InferenceType::CENTER;
    inf.origin = Eigen::Vector3d(parameters[0], parameters[1], parameters[2]);
    inf.z_axis = Eigen::Vector3d(parameters[3], parameters[4], parameters[5]);

    // Onshape Special Case: If a circle or ellipse is on a planar face,
    // orient with the face normal instead of the cirlce or ellipse axis
    PK_ERROR_t err = PK_ERROR_no_errors;
    int num_faces;
    PK_FACE_t* faces;

    PK_EDGE_ask_faces(_id, &num_faces, &faces);
    assert(err == PK_ERROR_no_errors); // PK_EDGE_ask_faces

    auto is_parallel = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
        double dot_product = a.dot(b);
        double diff = (1 - (dot_product * dot_product));
        double abs_diff = diff * diff;
        return (abs_diff < FLT_EPSILON);
    };

    for (int i = 0; i < num_faces; ++i) {
        PK_SURF_t surface;
        PK_LOGICAL_t oriented;
        PK_FACE_ask_oriented_surf(faces[i], &surface, &oriented);
        assert(err == PK_ERROR_no_errors); // PK_FACE_ask_oriented_surf
        PK_CLASS_t face_class;
        PK_ENTITY_ask_class(surface, &face_class);
        assert(err == PK_ERROR_no_errors); // PK_ENTITY_ask_class
        if (face_class == PK_CLASS_plane) {
            PK_PLANE_sf_t plane_sf;
            PK_PLANE_ask(surface, &plane_sf);
            assert(err == PK_ERROR_no_errors);
            Eigen::Vector3d plane_axis(
                plane_sf.basis_set.axis.coord[0],
                plane_sf.basis_set.axis.coord[1],
                plane_sf.basis_set.axis.coord[2]);
            if (!oriented) {
                plane_axis = -plane_axis;
            }

            // If parallel then the circle lies in the face
            if (is_parallel(plane_axis, inf.z_axis)) {
                // If parallel and not equal, then they are opposite
                // and we should flip the axis if using with Onshape
                if ((plane_axis.normalized() - 
                    inf.z_axis.normalized()).norm() > FLT_EPSILON) {
                    inf.flipped_in_onshape = true;
                }
                break; // Stop after we find the containing plane, if it exists
            }
            
        }
    }

    PK_MEMORY_free(faces); // Clean up PK_EDGE_ask_faces
    inferences.push_back(inf);
}

void Edge::add_inferences_line(std::vector<Inference>& inferences)
{
    Inference inf;
    inf.inference_type = InferenceType::MID_POINT;
    inf.origin = mid_point;
    inf.z_axis = Eigen::Vector3d(parameters[3], parameters[4], parameters[5]);

    inferences.push_back(inf);
}

void Edge::add_inferences_other(std::vector<Inference>& inferences)
{
    Inference inf;
    inf.inference_type = InferenceType::MID_POINT;
    inf.origin = mid_point;
    inf.z_axis = Eigen::Vector3d(0.0, 0.0, 1.0);

    inferences.push_back(inf);
}

void Edge::sample_points(const int num_points, const bool sample_tangents, std::vector<Eigen::VectorXd>& samples, Eigen::Vector2d& t_range)
{
    // t_start and t_end already account for edge "sense"
    t_range(0) = t_start;
    t_range(1) = t_end;

    if (!has_curve) { // No curve, return all 0s
        samples.push_back(Eigen::VectorXd::Zero(num_points)); // x
        samples.push_back(Eigen::VectorXd::Zero(num_points)); // y
        samples.push_back(Eigen::VectorXd::Zero(num_points)); // z
        if (sample_tangents) {
            samples.push_back(Eigen::VectorXd::Zero(num_points)); // t_x
            samples.push_back(Eigen::VectorXd::Zero(num_points)); // t_y
            samples.push_back(Eigen::VectorXd::Zero(num_points)); // t_z
            samples.push_back(Eigen::VectorXd::Zero(num_points)); // curvature
        }
    }
    else {

        Eigen::ArrayXd ts = Eigen::ArrayXd::LinSpaced(num_points, t_start, t_end);

        Eigen::VectorXd x, y, z, t_x, t_y, t_z, c;
        x.resize(num_points);
        y.resize(num_points);
        z.resize(num_points);
        if (sample_tangents) {
            t_x.resize(num_points);
            t_y.resize(num_points);
            t_z.resize(num_points);
            c.resize(num_points);
        }
        for (int i = 0; i < num_points; ++i) {
            double t = ts(i);
            PK_ERROR_t err = PK_ERROR_no_errors;
            PK_VECTOR_t point;
            err = PK_CURVE_eval(_curve, t, 0 /*# derivatives*/, &point);
            assert(err == PK_ERROR_no_errors); // PK_CURVE_eval
            x(i) = point.coord[0];
            y(i) = point.coord[1];
            z(i) = point.coord[2];

            if (sample_tangents) {
                PK_VECTOR1_t tangent, principal_normal, binormal;
                double curvature;
                err = PK_CURVE_eval_curvature(
                    _curve,
                    t,
                    &tangent,
                    &principal_normal,
                    &binormal,
                    &curvature
                );
                assert(err == PK_ERROR_no_errors ||
                        err == PK_ERROR_at_terminator ||
                        err == PK_ERROR_bad_parameter ||
                        err == PK_ERROR_eval_failure); // PK_CURVE_eval_curvature

                if (err == PK_ERROR_at_terminator ||
                    err == PK_ERROR_bad_parameter ||
                    err == PK_ERROR_eval_failure) {
                    t_x(i) = 0;
                    t_y(i) = 0;
                    t_z(i) = 0;
                    c(i) = 0;
                }
                else {
                    t_x(i) = tangent.coord[0];
                    t_y(i) = tangent.coord[1];
                    t_z(i) = tangent.coord[2];
                    c(i) = curvature;
                }
            }
        }

        samples.push_back(x);
        samples.push_back(y);
        samples.push_back(z);

        if (sample_tangents) {
            if (_is_reversed) { // Flip tangents if edge is reversed
                t_x = -t_x;
                t_y = -t_y;
                t_z = -t_z;
            }
            samples.push_back(t_x);
            samples.push_back(t_y);
            samples.push_back(t_z);
            samples.push_back(c);
        }
    }
}

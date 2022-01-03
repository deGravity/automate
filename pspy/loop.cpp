#include "loop.h"
#include "loop.h"

Loop::Loop(int id)
{
    _id = id;

    PK_ERROR_t err = PK_ERROR_no_errors;

    PK_LOOP_type_t loop_type;
    err = PK_LOOP_ask_type(_id, &loop_type);
    assert(err == PK_ERROR_no_errors); // PK_LOOP_ask_type
    // LOOP_type enum is range 5410 (vertex_c) - 5419
    // this subtraction puts it in range 0-9 instead
    _type = static_cast<LoopType>(loop_type - PK_LOOP_type_vertex_c);

    PK_EDGE_t* edges;
    int num_edges;
    err = PK_LOOP_ask_edges(_id, &num_edges, &edges);
    assert(err == PK_ERROR_no_errors);

    // A loop can be degenerate (e.g. the tip of a cone)
    // so we can't always get mass properties
    if (num_edges > 0) {

        // Check if this loop is just a single circle
        // this is needed in inference calculation
        if (num_edges == 1) {
            PK_EDGE_t _edge_id = edges[0];
            PK_LOGICAL_t want_interval = PK_LOGICAL_false;
            PK_CURVE_t curve;
            PK_CLASS_t curve_class;
            PK_VECTOR_t ends[2];
            PK_INTERVAL_t t_int;
            PK_LOGICAL_t sense;
            err = PK_EDGE_ask_geometry(
                _edge_id, 
                want_interval, 
                &curve, 
                &curve_class, 
                ends, 
                &t_int, 
                &sense);
            if (curve_class == PK_CLASS_circle) {
                _is_circle = true;
            }
        }
        else {
            _is_circle = false;
        }

        auto m = MassProperties((int*)edges, MASS_ACC, num_edges);
        length = m.amount;
        center_of_gravity = m.c_of_g;
        moment_of_inertia = m.m_of_i;

        // Get Non-Aligned Bounding Box
        // Axes are ordered so X is largest, Y is second largest, Z is smallest
        PK_NABOX_sf_t nabox;
        PK_TOPOL_find_nabox_o_t nabox_options;
        PK_TOPOL_find_nabox_o_m(nabox_options);
        err = PK_TOPOL_find_nabox(num_edges, edges, NULL, &nabox_options, &nabox);
        assert(err == PK_ERROR_no_errors || err == PK_ERROR_missing_geom); // PK_TOPOL_find_nabox

        if (err == PK_ERROR_missing_geom) { // At least one curve is missing geometry
            na_bb_center.setZero();
            na_bb_x.setZero();
            na_bb_z.setZero();
            na_bounding_box = Eigen::MatrixXd::Zero(2, 3);
        }
        else {
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
    }
    else {
        length = 0;
        center_of_gravity = Eigen::Vector3d(0, 0, 0);
        moment_of_inertia = Eigen::MatrixXd::Zero(3, 3);

        // We could techincally get this for vertices if they exist
        na_bb_center = Eigen::Vector3d(0, 0, 0);
        na_bb_x = Eigen::Vector3d(0, 0, 0);
        na_bb_z = Eigen::Vector3d(0, 0, 0);
        na_bounding_box = Eigen::MatrixXd::Zero(2, 3);
    }

    PK_MEMORY_free(edges);
}

std::vector<Inference> Loop::get_inferences()
{
    std::vector<Inference> inferences;

    // Only make inferences for inner loops
    if (_type != LoopType::OUTER && _type != LoopType::LIKELY_OUTER) {

        // Only make inferences for loops in planes
        PK_ERROR_t err = PK_ERROR_no_errors;
        PK_FACE_t face;
        err = PK_LOOP_ask_face(_id, &face);
        assert(err == PK_ERROR_no_errors); // PK_LOOP_ask_face
        PK_SURF_t surf;
        PK_LOGICAL_t oriented;
        err = PK_FACE_ask_oriented_surf(face, &surf, &oriented);
        assert(err == PK_ERROR_no_errors); // PK_FACE_ask_surf
        PK_CLASS_t face_class;
        if (surf == PK_ENTITY_null) {
            return inferences;
        }
        err = PK_ENTITY_ask_class(surf, &face_class);
        assert(err == PK_ERROR_no_errors);
        if (face_class == PK_CLASS_plane) {

            // Get the plane's axis to use
            PK_PLANE_sf_t plane_sf;
            PK_PLANE_ask(surf, &plane_sf);
            assert(err == PK_ERROR_no_errors);
            Eigen::Vector3d plane_axis(
                plane_sf.basis_set.axis.coord[0],
                plane_sf.basis_set.axis.coord[1],
                plane_sf.basis_set.axis.coord[2]);
            if (!oriented) {
                plane_axis = -plane_axis;
            }

            Inference inf;
            inf.inference_type = InferenceType::LOOP_CENTER;
            inf.origin = center_of_gravity;
            inf.z_axis = plane_axis;

            // Onshape Special Case - Onshape does not allow
            // loop center inferences for circular loops
            if (_is_circle) {
                inf.onshape_inference = false;
            }


            inferences.push_back(inf);
        }
    }

    return inferences;
}

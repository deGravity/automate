#include "face.h"
#include "face.h"
#include "face.h"
#include "face.h"
#include "face.h"
#include "face.h"
#include <parasolid.h>
#include <assert.h>

Face::Face(int id) {
    _id = id;

    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_SURF_t surface;
    PK_LOGICAL_t oriented;

    err = PK_FACE_ask_oriented_surf(_id, &surface, &oriented);
    assert(err == PK_ERROR_no_errors); // PK_FACE_ask_oriented_surf

    _surf = surface;
    orientation = oriented;

    init_parametric_function();

    init_bb();

    init_nabb();

    init_mass_props();

    // TODO - Inferences
}

void Face::init_parametric_function() {
    PK_ERROR_t err = PK_ERROR_no_errors;

    PK_CLASS_t face_class;
    err = PK_ENTITY_ask_class(_surf, &face_class);
    assert(err == PK_ERROR_no_errors); //PK_ENTITY_ask_class

    switch (face_class) {
    case PK_CLASS_plane:
        init_plane();
        break;
    case PK_CLASS_cyl:
        init_cyl();
        break;
    case PK_CLASS_cone:
        init_cone();
        break;
    case PK_CLASS_sphere:
        init_sphere();
        break;
    case PK_CLASS_torus:
        init_torus();
        break;
    case PK_CLASS_spun:
        init_spun();
        break;
    case PK_CLASS_bsurf:
        function = SurfaceFunction::BSURF;
        break;
    case PK_CLASS_offset:
        function = SurfaceFunction::OFFSET;
        break;
    case PK_CLASS_swept:
        function = SurfaceFunction::SWEPT;
        break;
    case PK_CLASS_blendsf:
        function = SurfaceFunction::BLENDSF;
        break;
    case PK_CLASS_mesh:
        function = SurfaceFunction::MESH;
        break;
    case PK_CLASS_fsurf:
        function = SurfaceFunction::FSURF;
        break;
    }
}

void Face::init_bb() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    // Get Bounding Box
    PK_BOX_t box;
    err = PK_TOPOL_find_box(_id, &box);
    assert(err == PK_ERROR_no_errors);
    bounding_box.resize(2, 3);
    bounding_box <<
        box.coord[0], box.coord[1], box.coord[2],
        box.coord[3], box.coord[3], box.coord[5];
}

void Face::init_nabb() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    // Get Non-Aligned Bounding Box
    // Axes are ordered so X is largest, Y is second largest, Z is smallest
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

void Face::init_mass_props() {
    auto m = MassProperties(&_id);
    surface_area = m.amount;
    center_of_gravity = m.c_of_g;
    moment_of_inertia = m.m_of_i;
    circumference = m.periphery;
}

void Face::init_plane() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_PLANE_sf_t plane_sf;
    PK_PLANE_ask(_surf, &plane_sf);
    assert(err == PK_ERROR_no_errors); // PK_PLANE_ask
    function = SurfaceFunction::PLANE;
    parameters.push_back(plane_sf.basis_set.location.coord[0]);
    parameters.push_back(plane_sf.basis_set.location.coord[1]);
    parameters.push_back(plane_sf.basis_set.location.coord[2]);
    parameters.push_back(plane_sf.basis_set.axis.coord[0]);
    parameters.push_back(plane_sf.basis_set.axis.coord[1]);
    parameters.push_back(plane_sf.basis_set.axis.coord[2]);
}

void Face::init_cyl() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_CYL_sf_t cylinder;
    err = PK_CYL_ask(_surf, &cylinder);
    assert(err == PK_ERROR_no_errors); // PK_CYL_ask
    function = SurfaceFunction::CYLINDER;
    parameters.push_back(cylinder.basis_set.location.coord[0]);
    parameters.push_back(cylinder.basis_set.location.coord[1]);
    parameters.push_back(cylinder.basis_set.location.coord[2]);
    parameters.push_back(cylinder.basis_set.axis.coord[0]);
    parameters.push_back(cylinder.basis_set.axis.coord[1]);
    parameters.push_back(cylinder.basis_set.axis.coord[2]);
    parameters.push_back(cylinder.radius);
}

void Face::init_cone() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_CONE_sf_t cone;
    err = PK_CONE_ask(_surf, &cone);
    assert(err == PK_ERROR_no_errors); // PK_CONE_ask
    function = SurfaceFunction::CONE;
    parameters.push_back(cone.basis_set.location.coord[0]); // 0 - origin_x
    parameters.push_back(cone.basis_set.location.coord[1]); // 1 - origin_y
    parameters.push_back(cone.basis_set.location.coord[2]); // 2 - origin_z
    parameters.push_back(cone.basis_set.axis.coord[0]);     // 3 - axis_x
    parameters.push_back(cone.basis_set.axis.coord[1]);     // 4 - axis_y
    parameters.push_back(cone.basis_set.axis.coord[2]);     // 5 - axis_z
    parameters.push_back(cone.radius);                      // 6 - radius
    parameters.push_back(cone.semi_angle);                  // 7 - semi-angle
}

void Face::init_sphere() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_SPHERE_sf_t sphere;
    PK_SPHERE_ask(_surf, &sphere);
    assert(err == PK_ERROR_no_errors); // PK_SPHERE_ask
    function = SurfaceFunction::SPHERE;
    parameters.push_back(sphere.basis_set.location.coord[0]);
    parameters.push_back(sphere.basis_set.location.coord[1]);
    parameters.push_back(sphere.basis_set.location.coord[2]);
    parameters.push_back(sphere.basis_set.axis.coord[0]);
    parameters.push_back(sphere.basis_set.axis.coord[1]);
    parameters.push_back(sphere.basis_set.axis.coord[2]);
    parameters.push_back(sphere.radius);
}

void Face::init_torus() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_TORUS_sf_t torus;
    err = PK_TORUS_ask(_surf, &torus);
    assert(err == PK_ERROR_no_errors); // PK_CONE_ask
    function = SurfaceFunction::TORUS;
    parameters.push_back(torus.basis_set.location.coord[0]);
    parameters.push_back(torus.basis_set.location.coord[1]);
    parameters.push_back(torus.basis_set.location.coord[2]);
    parameters.push_back(torus.basis_set.axis.coord[0]);
    parameters.push_back(torus.basis_set.axis.coord[1]);
    parameters.push_back(torus.basis_set.axis.coord[2]);
    parameters.push_back(torus.major_radius);
    parameters.push_back(torus.minor_radius);
}

void Face::init_spun() {
    PK_ERROR_t err = PK_ERROR_none;
    PK_SPUN_sf_t spun;
    err = PK_SPUN_ask(_surf, &spun);
    assert(err == PK_ERROR_no_errors); // PK_SPUN_ask
    function = SurfaceFunction::SPUN;
    parameters.push_back(spun.axis.location.coord[0]);
    parameters.push_back(spun.axis.location.coord[1]);
    parameters.push_back(spun.axis.location.coord[2]);
    parameters.push_back(spun.axis.axis.coord[0]);
    parameters.push_back(spun.axis.axis.coord[1]);
    parameters.push_back(spun.axis.axis.coord[2]);
}
std::vector<Inference> Face::get_inferences()
{
    std::vector<Inference> inferences;
    switch (function) {
    case SurfaceFunction::PLANE:
        add_inferences_plane(inferences);
        break;
    case SurfaceFunction::CYLINDER:
        add_inferences_axial(inferences);
        break;
    case SurfaceFunction::CONE:
        add_inferences_cone(inferences);
        add_inferences_axial(inferences);
        break;
    case SurfaceFunction::SPHERE:
        add_inferences_sphere(inferences);
        break;
    case SurfaceFunction::TORUS:
        add_inferences_axial(inferences);
        break;
    case SurfaceFunction::SPUN:
        add_inferences_axial(inferences);
        break;
    default:
        break;
    }
    return inferences;
}
void Face::add_inferences_plane(std::vector<Inference>& inferences)
{
    // Centroid of plane, oriented with the _face_ normal
    // (so flipped if the plane geometry is also)
    Inference inf;
    inf.inference_type = InferenceType::CENTROID;
    inf.origin = Eigen::Vector3d(
        center_of_gravity[0],
        center_of_gravity[1],
        center_of_gravity[2]
    );
    inf.z_axis = Eigen::Vector3d(
        parameters[3],
        parameters[4],
        parameters[5]
    );
    if (!orientation) {
        inf.z_axis = -inf.z_axis;
    }
    inferences.push_back(inf);
}
void Face::add_inferences_cone(std::vector<Inference>& inferences)
{
    // Tip of cone, oriented along cone axis
    double radius = parameters[6];
    double semi_angle = parameters[7];
    Inference tip_inf;
    tip_inf.inference_type = InferenceType::POINT;

    double length = radius / tan(semi_angle);
    for (int i = 0; i < 3; ++i) {
        // parmaters[i] is i-th origin coord
        // parameters[i+3] is i-th axis coord
        tip_inf.origin(i) =
            parameters[i] - parameters[i + 3] * length;
        tip_inf.z_axis(i) = parameters[i + 3];
    }
    inferences.push_back(tip_inf);
}
void Face::add_inferences_sphere(std::vector<Inference>& inferences)
{
    // Center of sphere, oriented with z-axis
    Inference inf;
    inf.inference_type = InferenceType::CENTER;
    inf.origin = Eigen::Vector3d(parameters[0], parameters[1], parameters[2]);
    inf.z_axis = Eigen::Vector3d(0.0, 0.0, 1.0);
    inferences.push_back(inf);
}
void Face::add_inferences_axial(std::vector<Inference>& inferences)
{
    // Top, Mid-Point, and Bottom of face as projected onto central
    // axis, oriented parallel to the central axis

    // Find the limits of the face parameterization. Used to determine
    // where the extreme and mid-points are
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_UVBOX_t uv_box;
    err = PK_FACE_find_uvbox(_id, &uv_box);
    assert(err == PK_ERROR_no_errors); // PK_FACE_find_uvbox

    // Extract the UV box corners 
    PK_UV_t min_uv, max_uv, center_uv;
    min_uv.param[0] = uv_box.param[0];
    min_uv.param[1] = uv_box.param[1];
    max_uv.param[0] = uv_box.param[2];
    max_uv.param[1] = uv_box.param[3];
    center_uv.param[0] = (min_uv.param[0] + max_uv.param[0]) / 2;
    center_uv.param[1] = (min_uv.param[1] + max_uv.param[1]) / 2;

    // Evaluate the surface function at the extreme and center
    // coordinates for projection onto central axis
    PK_VECTOR_t min_surf_coord, max_surf_coord, center_surf_coord;
    err = PK_SURF_eval(_surf, min_uv, 0, 0, PK_LOGICAL_false, &min_surf_coord);
    assert(err == PK_ERROR_no_errors); // PK_SURF_eval - min coord
    err = PK_SURF_eval(_surf, max_uv, 0, 0, PK_LOGICAL_false, &max_surf_coord);
    assert(err == PK_ERROR_no_errors); // PK_SURF_eval - max coord
    err = PK_SURF_eval(_surf, center_uv, 0, 0, PK_LOGICAL_false, &center_surf_coord);
    assert(err == PK_ERROR_no_errors); // PK_SURF_eval - center coord

    Eigen::Vector3d surf_min = Eigen::Map<Eigen::Vector3d>(min_surf_coord.coord);
    Eigen::Vector3d surf_max = Eigen::Map<Eigen::Vector3d>(max_surf_coord.coord);
    Eigen::Vector3d surf_center = Eigen::Map<Eigen::Vector3d>(center_surf_coord.coord);

    auto project_to_line = [](
        const Eigen::Vector3d& point, 
        const Eigen::Vector3d& origin, 
        const Eigen::Vector3d& axis) {
            Eigen::Vector3d v = point - origin;
            double length = v.dot(axis) / axis.dot(axis);
            return origin + axis * length;
    };

    // Origin and axis x,y,z components are the first 6 parameters
    Eigen::Vector3d origin = Eigen::Vector3d(parameters[0], parameters[1], parameters[2]);
    Eigen::Vector3d axis = Eigen::Vector3d(parameters[3], parameters[4], parameters[5]);

    Eigen::Vector3d min_axis_point = project_to_line(surf_min, origin, axis);
    Eigen::Vector3d max_axis_point = project_to_line(surf_max, origin, axis);
    Eigen::Vector3d center_axis_point = project_to_line(surf_center, origin, axis);

    Inference bottom_axis, top_axis, mid_axis;
    
    bottom_axis.inference_type = InferenceType::BOTTOM_AXIS_POINT;
    bottom_axis.origin = min_axis_point;
    bottom_axis.z_axis = axis;

    top_axis.inference_type = InferenceType::TOP_AXIS_POINT;
    top_axis.origin = max_axis_point;
    top_axis.z_axis = axis;

    mid_axis.inference_type = InferenceType::MID_AXIS_POINT;
    mid_axis.origin = center_axis_point;
    mid_axis.z_axis = axis;

    inferences.push_back(bottom_axis);
    inferences.push_back(top_axis);
    inferences.push_back(mid_axis);
}

void Face::sample_points(
    const int num_points,
    const bool sample_normals,
    std::vector<Eigen::MatrixXd>& samples,
    Eigen::MatrixXd& uv_range) {

    // Find UV bounding box
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_UVBOX_t uv_box;
    err = PK_FACE_find_uvbox(_id, &uv_box);
    assert(err == PK_ERROR_no_errors);
    uv_range.resize(2, 2);
    uv_range <<
        uv_box.param[0], uv_box.param[1],
        uv_box.param[2], uv_box.param[3];

    // Return Empty Samples if there is no surface
    if (_surf == PK_ENTITY_null) {
        samples.push_back(Eigen::MatrixXd::Zero(num_points, num_points)); // x
        samples.push_back(Eigen::MatrixXd::Zero(num_points, num_points)); // y
        samples.push_back(Eigen::MatrixXd::Zero(num_points, num_points)); // z
        samples.push_back(Eigen::MatrixXd::Zero(num_points, num_points)); // mask
        if (sample_normals) {
            samples.push_back(Eigen::MatrixXd::Zero(num_points, num_points)); // n_x
            samples.push_back(Eigen::MatrixXd::Zero(num_points, num_points)); // n_y
            samples.push_back(Eigen::MatrixXd::Zero(num_points, num_points)); // n_y
            samples.push_back(Eigen::MatrixXd::Zero(num_points, num_points)); // pc_1
            samples.push_back(Eigen::MatrixXd::Zero(num_points, num_points)); // pc_2
        }
        return;
    }

    // Generated Equally Spaced Sample Points
    Eigen::ArrayXd u_samples = Eigen::ArrayXd::LinSpaced(
        num_points, uv_box.param[0], uv_box.param[2]);
    Eigen::ArrayXd v_samples = Eigen::ArrayXd::LinSpaced(
        num_points, uv_box.param[1], uv_box.param[3]);

    // Evaluate surface on this grid
    PK_VECTOR_t* points = new PK_VECTOR_t[num_points * num_points];
    err = PK_SURF_eval_grid(
        _surf, // surface 
        num_points, // # u parameters
        u_samples.data(), // u parameters
        num_points, // # v parameters
        v_samples.data(), // v parameters
        0, // # u derivatives (must be 0)
        0, // # v derivatives (must be 0)
        PK_LOGICAL_false, // triangular derivative array required
        points
    );
    assert(err == PK_ERROR_no_errors); // PK_SURF_eval_grid

    for (int i = 0; i < 3; ++i) { // x,y,z
        Eigen::MatrixXd s_i(num_points, num_points);
        int k = 0;
        for (int v = 0; v < num_points; ++v) {
            for (int u = 0; u < num_points; ++u) {
                s_i(u, v) = points[k].coord[i];
                ++k;
            }
        }
        samples.push_back(s_i);
    }

    // Evaluate Normals and Curvature if requested
    if (sample_normals) {
        Eigen::MatrixXd n_x, n_y, n_z;
        n_x.resize(num_points, num_points);
        n_y.resize(num_points, num_points);
        n_z.resize(num_points, num_points);
        Eigen::MatrixXd pc_grid_1, pc_grid_2;
        pc_grid_1.resize(num_points, num_points);
        pc_grid_2.resize(num_points, num_points);
        for (int i = 0; i < num_points; ++i) {
            for (int j = 0; j < num_points; ++j) {
                double u = u_samples[i];
                double v = v_samples[j];
                PK_UV_t uv;
                uv.param[0] = u;
                uv.param[1] = v;
                PK_VECTOR_t normal, pd_1, pd_2;
                double pc_1, pc_2;
                err = PK_SURF_eval_curvature(
                    _surf,
                    uv,
                    &normal,
                    &pd_1,
                    &pd_2,
                    &pc_1,
                    &pc_2
                );
                assert(err == PK_ERROR_no_errors); // PK_SURF_eval_curvature
                n_x(i, j) = normal.coord[0];
                n_y(i, j) = normal.coord[1];
                n_z(i, j) = normal.coord[2];
                pc_grid_1(i, j) = pc_1;
                pc_grid_2(i, j) = pc_2;
            }
        }
        if (!orientation) { // If the surface is reversed, flip normals and pcs
            n_x = -1 * n_x;
            n_y = -1 * n_y;
            n_z = -1 * n_z;
            pc_grid_1 = -1 * pc_grid_1;
            pc_grid_2 = -1 * pc_grid_2;
        }
        samples.push_back(n_x);
        samples.push_back(n_y);
        samples.push_back(n_z);
        samples.push_back(pc_grid_1);
        samples.push_back(pc_grid_2);
    }

    // Check the grid points to see if they are in the face or not
    PK_TOPOL_t* point_topos = new PK_TOPOL_t[num_points * num_points];
    PK_FACE_contains_vectors_o_t contains_vectors_opt;
    PK_FACE_contains_vectors_o_m(contains_vectors_opt);
    contains_vectors_opt.is_on_surf = PK_LOGICAL_true;
    contains_vectors_opt.n_vectors = 100;
    contains_vectors_opt.vectors = points;
    err = PK_FACE_contains_vectors(_id, &contains_vectors_opt, point_topos);
    assert(err == PK_ERROR_no_errors); // PK_FACE_contains_vectors

    Eigen::MatrixXd mask(num_points, num_points);
    int k = 0;
    for (int v = 0; v < num_points; ++v) {
        for (int u = 0; u < num_points; ++u) {
            mask(u, v) = (point_topos[k] == NULL) ? 0.0 : 1.0;
            ++k;
        }
    }
    samples.push_back(mask);

    delete[] points;
    delete[] point_topos;
}

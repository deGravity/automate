#include "body.h"

#include <set>
#include <map>
#include <vector>

bool is_body(int id) {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_CLASS_t entity_class;
    err = PK_ENTITY_ask_class(id, &entity_class);
    return (err == 0) && (entity_class == PK_CLASS_body);
}

std::vector<Body> read_xt(std::string path) {
    ensure_parasolid_session();
    PK_PART_receive_o_t receive_opts;
    PK_PART_receive_o_m(receive_opts);
    receive_opts.transmit_format = PK_transmit_format_text_c;
    int n_parts = 0;
    PK_PART_t* parts = NULL;
    PK_ERROR_t err = PK_ERROR_no_errors;
    err = PK_PART_receive(path.c_str(), &receive_opts, &n_parts, &parts);
    std::vector<Body> parts_vec;
    if (err == 0) {
        for (int i = 0; i < n_parts; ++i) {
            if (is_body(parts[i])) {
                parts_vec.emplace_back(parts[i]);
            }
            else {
                int num_deleted;
                PK_ENTITY_delete_attribs(parts[i], PK_ENTITY_null, &num_deleted);
                PK_ENTITY_delete(1, &parts[i]);
            }
        }
    }
    PK_MEMORY_free(parts); // Do I need to do this, or is this causing segfaults?
    return parts_vec;
}

Body::Body(int id) {
    _id = id;
    _valid = true;
}

Body::~Body() {
    int num_deleted;
    PK_ENTITY_delete_attribs(_id, PK_ENTITY_null, &num_deleted);
    PK_ENTITY_delete(1, &_id);
}

BREPTopology Body::GetTopology() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_BODY_ask_topology_o_t ask_topology_options;
    PK_BODY_ask_topology_o_m(ask_topology_options);
    ask_topology_options.want_fins = PK_LOGICAL_false;
    PK_TOPOL_t* topols;
    PK_CLASS_t* classes;
    int n_topols;
    int n_relations;
    int* parents;
    int* children;
    PK_TOPOL_sense_t* senses;
    err = PK_BODY_ask_topology(
        _id,
        &ask_topology_options,
        &n_topols,
        &topols,
        &classes,
        &n_relations,
        &parents,
        &children,
        &senses);
    assert(err == PK_ERROR_no_errors); // PK_BODY_ask_topology
    BREPTopology topology;

    std::map<int, int> cat_idx;

    for (int i = 0; i < n_topols; ++i) {
        topology.pk_to_class[topols[i]] = classes[i];
        switch (classes[i]) {
        case PK_CLASS_face:
            cat_idx[i] = topology.faces.size();
            topology.faces.emplace_back(topols[i]);
            break;
        case PK_CLASS_loop:
            cat_idx[i] = topology.loops.size();
            topology.loops.emplace_back(topols[i]);
            break;
        case PK_CLASS_edge:
            cat_idx[i] = topology.edges.size();
            topology.edges.emplace_back(topols[i]);
            break;
        case PK_CLASS_vertex:
            cat_idx[i] = topology.vertices.size();
            topology.vertices.emplace_back(topols[i]);
            break;
        default:
            break;
        }
        // Update map from pk index to index with each entity type
        topology.pk_to_idx[topols[i]] = cat_idx[i];
    }

    // Make loopups for relation types instead of branching
    const int FACE_LOOP = PK_CLASS_face + 2 * PK_CLASS_loop;
    const int LOOP_EDGE = PK_CLASS_loop + 2 * PK_CLASS_edge;
    const int EDGE_VERTEX = PK_CLASS_edge + 2 * PK_CLASS_vertex;
    const int LOOP_VERTEX = PK_CLASS_loop + 2 * PK_CLASS_vertex; // degenerate loops

    // Maps for finding face-face edges
    std::map<int, int> loop_to_face;
    std::map<int, std::vector<int> > edge_to_loop;

    std::map<int, std::vector<int> > face_to_loops;
    std::map<int, std::vector<int> > loop_to_edges;
    std::map<int, std::vector<int> > loop_to_vertex;
    std::map<int, std::vector<int> > edge_to_vertices;
    for (int i = 0; i < n_relations; ++i) {
        // Index with their type vectors of parent and child
        int parent = cat_idx[parents[i]];
        int child = cat_idx[children[i]];
        int parent_class = classes[parents[i]];
        int child_class = classes[children[i]];
        int relation_type = parent_class + 2 * child_class;

        switch (relation_type) {
        case FACE_LOOP:
            loop_to_face[child] = parent;
            topology.face_to_loop.emplace_back(parent, child, senses[i]);
            face_to_loops[parent].push_back(child);
            break;
        case LOOP_EDGE:
            edge_to_loop[child].push_back(parent);
            topology.loop_to_edge.emplace_back(parent, child, senses[i]);
            loop_to_edges[parent].push_back(child);
            break;
        case EDGE_VERTEX:
            topology.edge_to_vertex.emplace_back(parent, child, senses[i]);
            edge_to_vertices[parent].push_back(child);
            break;
        case LOOP_VERTEX:
            topology.loop_to_vertex.emplace_back(parent, child, senses[i]);
            loop_to_vertex[parent].push_back(child);
            break;
        default:
            // Ignore relation types we don't use
            break;
        }
    }

    // Also find Face-Face Edges
    for (auto edgeloop : edge_to_loop) {
        int edge = edgeloop.first;
        auto loops = edgeloop.second;
        assert(loops.size() == 2);
        int face1 = loop_to_face[loops[0]];
        int face2 = loop_to_face[loops[1]];
        topology.face_to_face.emplace_back(face1, face2, edge);
    }

    // Construct Adjacency List Maps

    // Compute Lower Face Adjacencies
    std::map<int, std::vector<int> > face_to_edges;
    std::map<int, std::vector<int> > face_to_vertices;
    for (int face = 0; face < topology.faces.size(); ++face) {
        std::set<int> edge_neighbors;
        std::set<int> vertex_neighbors;
        for (int loop : face_to_loops[face]) {
            for (int edge : loop_to_edges[loop]) {
                edge_neighbors.insert(edge);
                for (int vertex : edge_to_vertices[edge]) {
                    vertex_neighbors.insert(vertex);
                }
            }
        }
        for (int edge : edge_neighbors) {
            face_to_edges[face].push_back(edge);
        }
        for (int vertex : vertex_neighbors) {
            face_to_vertices[face].push_back(vertex);
        }
    }

    // Compute Lower Loop Adjacencies
    std::map<int, std::vector<int> > loop_to_vertices;
    for (int loop = 0; loop < topology.loops.size(); ++loop) {
        std::set<int> vertex_neighbors;
        for (int edge : loop_to_edges[loop]) {
            for (int vertex : edge_to_vertices[edge]) {
                vertex_neighbors.insert(vertex);
            }
        }
        for (int vertex : vertex_neighbors) {
            loop_to_vertices[loop].push_back(vertex);
        }
    }

    // Assign to structure
    // TODO - don't use the temporary variables
    // TODO - int -> size_t
    topology.face_loop = face_to_loops;
    topology.face_edge = face_to_edges;
    topology.face_vertex = face_to_vertices;
    topology.loop_edge = edge_to_loop;
    topology.loop_vertex = loop_to_vertices;
    topology.edge_vertex = edge_to_vertices;

    // Clean Up PK_BODY_ask_topology
    PK_MEMORY_free(topols);
    PK_MEMORY_free(classes);
    PK_MEMORY_free(parents);
    PK_MEMORY_free(children);
    PK_MEMORY_free(senses);

    return topology;
}

MassProperties Body::GetMassProperties(double accuracy) {
    return MassProperties(&_id, accuracy);
}

Eigen::MatrixXd Body::GetBoundingBox() {
    PK_ERROR_code_t err = PK_ERROR_no_errors;
    PK_BOX_t box;
    err = PK_TOPOL_find_box(_id, &box);
    assert(err == PK_ERROR_no_errors); // PK_TOPOL_find_box
    Eigen::MatrixXd corners(2, 3);
    corners <<
        box.coord[0], box.coord[1], box.coord[2],
        box.coord[3], box.coord[3], box.coord[5];
    return corners;
}

int Body::Transform(const Eigen::MatrixXd& xfrm) {
    // Apply Transform
    PK_BODY_transform_o_t transform_options;
    PK_BODY_transform_o_m(transform_options);
    PK_TRANSF_t transformation;
    PK_TRANSF_sf_t transform_mat;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            transform_mat.matrix[i][j] = xfrm(i, j);
        }
    }
    PK_TOPOL_track_r_t tracking;
    PK_TOPOL_local_r_t results;
    PK_TRANSF_create(&transform_mat, &transformation);
    PK_ERROR_t err = PK_ERROR_no_errors;
    err = PK_BODY_transform_2(_id, transformation, XFRM_TOL, &transform_options, &tracking, &results);

    if (!(err == PK_ERROR_no_errors &&
        (results.status == PK_local_status_ok_c ||
            results.status == PK_local_status_nocheck_c)))
    {
        err = 1;
    }
    // Free Memory
    PK_TOPOL_track_r_f(&tracking);
    PK_TOPOL_local_r_f(&results);

    return err;
}

void Body::Tesselate(
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::VectorXi& FtoT,
    Eigen::MatrixXi& EtoT,
    Eigen::VectorXi& VtoT) {
    PK_ERROR_t err = PK_ERROR_no_errors;

    // Setup faceting call options
    PK_TOPOL_facet_2_o_t facet_options;
    PK_TOPOL_facet_2_o_m(facet_options);
    facet_options.choice.data_curv_idx = PK_LOGICAL_true;
    facet_options.choice.fin_data = PK_LOGICAL_true;
    facet_options.choice.data_param_idx = PK_LOGICAL_true;
    facet_options.choice.facet_face = PK_LOGICAL_true;
    facet_options.choice.facet_fin = PK_LOGICAL_true;
    facet_options.choice.fin_edge = PK_LOGICAL_true;
    facet_options.choice.fin_topol = PK_LOGICAL_true;
    facet_options.choice.data_point_idx = PK_LOGICAL_true;
    facet_options.choice.point_vec = PK_LOGICAL_true;
    facet_options.choice.point_topol = PK_LOGICAL_true;
    facet_options.choice.fin_edge = PK_LOGICAL_true;

    // Facet the body
    PK_TOPOL_facet_2_r_t facets;
    err = PK_TOPOL_facet_2(1, &_id, NULL, &facet_options, &facets);
    // TODO - Do we need more detailed checking (like for transforms)
    assert(err == PK_ERROR_no_errors); // PK_TOPOL_facet_2


    // Figure out which table is which, then point these names at them
    // for easier access
    PK_TOPOL_fctab_point_vec_t* point_vec = NULL;
    PK_TOPOL_fctab_facet_fin_t* facet_fin = NULL;
    PK_TOPOL_fctab_fin_data_t* fin_data = NULL;
    PK_TOPOL_fctab_data_point_t* data_point_idx = NULL;
    PK_TOPOL_fctab_facet_face_t* facet_face = NULL;
    PK_TOPOL_fctab_fin_topol_t* fin_topol = NULL;
    PK_TOPOL_fctab_point_topol_t* point_topol = NULL;
    PK_TOPOL_fctab_fin_edge_t* fin_edge = NULL;

    // Extract the tables we care about
    for (int i = 0; i < facets.number_of_tables; ++i) {
        switch (facets.tables[i].fctab) {
        case PK_TOPOL_fctab_point_vec_c:
            point_vec = facets.tables[i].table.point_vec;
            break;
        case PK_TOPOL_fctab_facet_fin_c:
            facet_fin = facets.tables[i].table.facet_fin;
            break;
        case PK_TOPOL_fctab_fin_data_c:
            fin_data = facets.tables[i].table.fin_data;
            break;
        case PK_TOPOL_fctab_data_point_c:
            data_point_idx = facets.tables[i].table.data_point_idx;
            break;
        case PK_TOPOL_fctab_facet_face_c:
            facet_face = facets.tables[i].table.facet_face;
            break;
        case PK_TOPOL_fctab_fin_topol_c:
            fin_topol = facets.tables[i].table.fin_topol;
            break;
        case PK_TOPOL_fctab_point_topol_c:
            point_topol = facets.tables[i].table.point_topol;
            break;
        case PK_TOPOL_fctab_fin_edge_c:
            fin_edge = facets.tables[i].table.fin_edge;
        }
    }


    // TODO - Previously we guarded this with facets.number_of_tables > 0
    // will error-catching deal with this case?

    // Populate Mesh Vertices
    V.resize(point_vec->length, 3);
    for (int i = 0; i < point_vec->length; ++i) {
        auto vec = point_vec->vec[i];
        V(i, 0) = vec.coord[0];
        V(i, 1) = vec.coord[1];
        V(i, 2) = vec.coord[2];
    }

    // Populate Mesh Faces
    F.resize(facet_fin->length / 3, 3);
    Eigen::VectorXi facet_offsets = Eigen::VectorXi::Zero(F.rows());
    Eigen::MatrixXi fin_indices(facet_fin->length, 2);
    for (int i = 0; i < facet_fin->length; ++i) {
        int facet_idx = facet_fin->data[i].facet;
        int fin_idx = facet_fin->data[i].fin;
        int data_idx = fin_data->data[fin_idx];
        int point_idx = data_point_idx->point[data_idx];
        F(facet_idx, facet_offsets[facet_idx]) = point_idx;
        fin_indices(i, 0) = facet_idx;
        fin_indices(i, 1) = facet_offsets[facet_idx];
        facet_offsets[facet_idx] += 1;
    }

    // PK_ENTITY_NULL == 0 - so use this as a sentinel topo reference

    // Populate Mesh to Topology Face References
    FtoT.resize(facet_face->length);
    for (int i = 0; i < facet_face->length; ++i) {
        FtoT[i] = facet_face->face[i];
    }
    // Populate Mesh to Topology Edge References
    EtoT = Eigen::MatrixXi::Zero(F.rows(), 3);
    for (int i = 0; i < fin_edge->length; ++i) {
        int fin = fin_edge->data[i].fin;
        PK_EDGE_t edge = fin_edge->data[i].edge;
        EtoT(fin_indices(fin, 0), fin_indices(fin, 1)) = edge;
    }
    // Populate Mesh to Topology Vertex References
    VtoT = Eigen::VectorXi::Zero(V.rows());
    int num_assigned = 0;
    for (int i = 0; i < point_topol->length; ++i) {
        VtoT(point_topol->data[i].point) = point_topol->data[i].topol;
    }

    PK_TOPOL_facet_2_r_f(&facets);
}

void Body::debug() {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_CLASS_t entity_class;
    err = PK_ENTITY_ask_class(_id, &entity_class);
    std::cout << entity_class << std::endl;
    auto topo = GetTopology();
    auto mass = GetMassProperties();

    Eigen::MatrixXd xfrm(4, 4);
    xfrm <<
        1, 0, 0, 0.0,
        0, 1, 0, 0.0,
        0, 0, 1, 0.0,
        0, 0, 0, 0.5;
    err = Transform(xfrm);
    std::cout << "xfrm error = " << err << std::endl;
    auto topo_xfrmed = GetTopology();
    auto mass_xfrmed = GetMassProperties();


    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXi FtoT;
    Eigen::MatrixXi EtoT;
    Eigen::VectorXi VtoT;
    Tesselate(V, F, FtoT, EtoT, VtoT);

    std::cout << "Num Vertices: " << topo.vertices.size() << std::endl;
    std::cout << "Num face-face relations: " << topo.face_to_face.size() << std::endl;
    std::cout << "Num Tess Verts: " << V.rows() << std::endl;
    std::cout << "Num Tess Faces: " << F.rows() << std::endl;
    std::cout << "Mass Amount: " << mass.amount << std::endl;
    std::cout << "Mass Mass: " << mass.mass << std::endl;
    std::cout << "Mass c of G:" << std::endl << mass.c_of_g << std::endl;
    std::cout << "Mass m of i:" << std::endl << mass.m_of_i << std::endl;
    std::cout << "Num Tess Faces: " << F.rows() << std::endl;

}
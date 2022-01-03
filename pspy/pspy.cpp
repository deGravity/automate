#include <torch/extension.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "eclass.h"
#include "part.h"

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
	// part.h

	py::class_<PartOptions>(m, "PartOptions")
		.def(py::init<>())
		.def_readwrite("just_bb", &PartOptions::just_bb)
		.def_readwrite("normalize", &PartOptions::normalize)
		.def_readwrite("transform", &PartOptions::transform)
		.def_readwrite("transform_matrix", &PartOptions::transform_matrix)
		.def_readwrite("num_uv_samples", &PartOptions::num_uv_samples)
		.def_readwrite("sample_normals", &PartOptions::num_uv_samples)
		.def_readwrite("sample_tangents", &PartOptions::sample_tangents)
		.def_readwrite("tesselate", &PartOptions::tesselate)
		.def_readwrite("default_mcfs", &PartOptions::default_mcfs)
		.def_readwrite("default_mcfs_only_face_axes", &PartOptions::default_mcfs_only_face_axes)
		.def_readwrite("onshape_style", &PartOptions::onshape_style)
		.def_readwrite("collect_inferences", &PartOptions::collect_inferences);

	py::class_<Part>(m, "Part")
		.def(py::init<const std::string&>())
		.def(py::init<const std::string&, PartOptions>())
		.def_readonly("mesh", &Part::mesh)
		.def_readonly("mesh_topology", &Part::mesh_topology)
		.def_readonly("brep", &Part::brep)
		.def_readonly("samples", &Part::samples)
		.def_readonly("summary", &Part::summary)
		.def_readonly("inferences", &Part::inferences)
		.def_readonly("default_mcfs", &Part::default_mcfs)
		.def_readonly("is_valid", &Part::_is_valid);

	py::class_<Mesh>(m, "Mesh")
		.def_readonly("V", &Mesh::V)
		.def_readonly("F", &Mesh::F);

	py::class_<MeshTopology>(m, "MeshTopology")
		.def_readonly("face_to_topology", &MeshTopology::face_to_topology)
		.def_readonly("edge_to_topology", &MeshTopology::edge_to_topology)
		.def_readonly("point_to_topology", &MeshTopology::point_to_topology);

	py::class_<InferenceReference>(m, "InferenceReference")
		.def_readonly("reference_index", &InferenceReference::reference_index)
		.def_readonly("reference_type", &InferenceReference::reference_type)
		.def_readonly("inference_type", &InferenceReference::inference_type);

	py::class_<PartInference>(m, "PartInference")
		.def_readonly("origin", &PartInference::origin)
		.def_readonly("axis", &PartInference::axis)
		.def_readonly("onshape_inference", &PartInference::onshape_inference)
		.def_readonly("flipped_in_onshape", &PartInference::flipped_in_onshape)
		.def_readonly("reference", &PartInference::reference);

	py::class_<PartFace>(m, "PartFace")
		.def_readonly("index", &PartFace::index)
		.def_readonly("function", &PartFace::function)
		.def_readonly("parameters", &PartFace::parameters)
		.def_readonly("orientation", &PartFace::orientation)
		.def_readonly("bounding_box", &PartFace::bounding_box)
		.def_readonly("na_bounding_box", &PartFace::na_bounding_box)
		.def_readonly("surface_area", &PartFace::surface_area)
		.def_readonly("circumference", &PartFace::circumference)
		.def_readonly("center_of_gravity", &PartFace::center_of_gravity)
		.def_readonly("moment_of_inertia", &PartFace::moment_of_inertia)
		.def_readonly("loop_neighbors", &PartFace::loop_neighbors)
		.def_readonly("edge_neighbors", &PartFace::edge_neighbors)
		.def_readonly("vertex_neighbors", &PartFace::vertex_neighbors)
		.def_readonly("inferences", &PartFace::inferences);

	py::class_<PartLoop>(m, "PartLoop")
		.def_readonly("index", &PartLoop::index)
		.def_readonly("type", &PartLoop::type)
		.def_readonly("length", &PartLoop::length)
		.def_readonly("center_of_gravity", &PartLoop::center_of_gravity)
		.def_readonly("moment_of_inertia", &PartLoop::moment_of_inertia)
		.def_readonly("na_bounding_box", &PartLoop::na_bounding_box)
		.def_readonly("edge_neighbors", &PartLoop::edge_neighbors)
		.def_readonly("vertex_neighbors", &PartLoop::vertex_neighbors)
		.def_readonly("inferences", &PartLoop::inferences);

	py::class_<PartEdge>(m, "PartEdge")
		.def_readonly("index", &PartEdge::index)
		.def_readonly("function", &PartEdge::function)
		.def_readonly("parameters", &PartEdge::parameters)
		.def_readonly("t_range", &PartEdge::t_range)
		.def_readonly("start", &PartEdge::start)
		.def_readonly("end", &PartEdge::end)
		.def_readonly("is_periodic", &PartEdge::is_periodic)
		.def_readonly("mid_point", &PartEdge::mid_point)
		.def_readonly("length", &PartEdge::length)
		.def_readonly("center_of_gravity", &PartEdge::center_of_gravity)
		.def_readonly("moment_of_inertia", &PartEdge::moment_of_inertia)
		.def_readonly("bounding_box", &PartEdge::bounding_box)
		.def_readonly("na_bounding_box", &PartEdge::na_bounding_box)
		.def_readonly("vertex_neighbors", &PartEdge::vertex_neighbors)
		.def_readonly("inferences", &PartEdge::inferences);

	py::class_<PartVertex>(m, "PartVertex")
		.def_readonly("index", &PartVertex::index)
		.def_readonly("position", &PartVertex::position)
		.def_readonly("inferences", &PartVertex::inferences);

	py::class_<PartTopologyNodes>(m, "PartTopologyNodes")
		.def_readonly("faces", &PartTopologyNodes::faces)
		.def_readonly("loops", &PartTopologyNodes::loops)
		.def_readonly("edges", &PartTopologyNodes::edges)
		.def_readonly("vertices", &PartTopologyNodes::vertices);
	
	py::class_<PartTopologyRelations>(m, "PartTopologyRelations")
		.def_readonly("face_to_loop", &PartTopologyRelations::face_to_loop)
		.def_readonly("loop_to_edge", &PartTopologyRelations::loop_to_edge)
		.def_readonly("edge_to_vertex", &PartTopologyRelations::edge_to_vertex)
		.def_readonly("face_to_face", &PartTopologyRelations::face_to_face);

	py::class_<PartTopology>(m, "PartTopology")
		.def_readonly("nodes", &PartTopology::nodes)
		.def_readonly("relations", &PartTopology::relations);

	py::class_<PartSamples>(m, "PartSamples")
		.def_readonly("face_samples", &PartSamples::face_samples)
		.def_readonly("edge_samples", &PartSamples::edge_samples);

	py::class_<PartSummary>(m, "PartSummary")
		.def_readonly("bounding_box", &PartSummary::bounding_box)
		.def_readonly("volume", &PartSummary::volume)
		.def_readonly("mass", &PartSummary::mass)
		.def_readonly("center_of_gravity", &PartSummary::center_of_gravity)
		.def_readonly("moment_of_inertia", &PartSummary::moment_of_inertia)
		.def_readonly("surface_area", &PartSummary::surface_area)
		.def_readonly("topo_type_counts", &PartSummary::topo_type_counts)
		.def_readonly("surface_type_counts", &PartSummary::surface_type_counts)
		.def_readonly("curve_type_counts", &PartSummary::curve_type_counts)
		.def_readonly("loop_type_counts", &PartSummary::loop_type_counts)
		.def_readonly("fingerprint", &PartSummary::fingerprint);

	py::class_<PartUniqueInferences>(m, "PartUniqueInferences")
		.def_readonly("origins", &PartUniqueInferences::origins)
		.def_readonly("axes", &PartUniqueInferences::axes)
		.def_readonly("frames", &PartUniqueInferences::frames)
		.def_readonly("origin_references", &PartUniqueInferences::origin_references)
		.def_readonly("axes_references", &PartUniqueInferences::axes_references)
		.def_readonly("frame_references", &PartUniqueInferences::frame_references);

	py::class_<MCFReference>(m, "MCFReference")
		.def_readonly("origin_ref", &MCFReference::origin_ref)
		.def_readonly("axis_ref", &MCFReference::axis_ref);

	py::class_<MCF>(m, "MCF")
		.def_readonly("origin", &MCF::origin)
		.def_readonly("axis", &MCF::axis)
		.def_readonly("ref", &MCF::ref);

	// eclass.h
	m.def("find_equivalence_classes", &find_equivalence_classes);

	
	// body.h
	py::class_<Body>(m, "Body")
		.def("GetTopology", &Body::GetTopology)
		.def("GetMassProperties", &Body::GetMassProperties)
		.def("GetBoundingBox", &Body::GetBoundingBox)
		.def("Transform", &Body::Transform)
		.def("Tesselate", &Body::Tesselate)
		.def("__repr__",
			[](const Body& p) {
				return "<Parasolid Body>";
			});

	m.def("read_xt", &read_xt);

	

	// types.h

	py::enum_<TopologyType>(m, "TopologyType")
		.value("FACE", TopologyType::FACE)
		.value("EDGE", TopologyType::EDGE)
		.value("VERTEX", TopologyType::VERTEX)
		.value("LOOP", TopologyType::LOOP);

	py::enum_<SurfaceFunction>(m, "SurfaceFunction")
		.value("PLANE", SurfaceFunction::PLANE)
		.value("CYLINDER", SurfaceFunction::CYLINDER)
		.value("CONE", SurfaceFunction::CONE)
		.value("SPHERE", SurfaceFunction::SPHERE)
		.value("TORUS", SurfaceFunction::TORUS)
		.value("SPUN", SurfaceFunction::SPUN)
		.value("BSURF", SurfaceFunction::BSURF)
		.value("OFFSET", SurfaceFunction::OFFSET)
		.value("SWEPT", SurfaceFunction::SWEPT)
		.value("BLENDSF", SurfaceFunction::BLENDSF)
		.value("MESH", SurfaceFunction::MESH)
		.value("FSURF", SurfaceFunction::FSURF)
		.value("NONE", SurfaceFunction::NONE);

	py::enum_<CurveFunction>(m, "CurveFunction")
		.value("LINE", CurveFunction::LINE)
		.value("CIRCLE", CurveFunction::CIRCLE)
		.value("ELLIPSE", CurveFunction::ELLIPSE)
		.value("BCURVE", CurveFunction::BCURVE)
		.value("ICURVE", CurveFunction::ICURVE)
		.value("FCURVE", CurveFunction::FCURVE)
		.value("SPCURVE", CurveFunction::SPCURVE)
		.value("TRCURVE", CurveFunction::TRCURVE)
		.value("CPCURVE", CurveFunction::CPCURVE)
		.value("PLINE", CurveFunction::PLINE)
		.value("NONE", CurveFunction::NONE);

	py::enum_<LoopType>(m, "LoopType")
		.value("VERTEX", LoopType::VERTEX)
		.value("WIRE", LoopType::WIRE)
		.value("OUTER", LoopType::OUTER)
		.value("INNER", LoopType::INNER)
		.value("WINDING", LoopType::WINDING)
		.value("INNER_SING", LoopType::INNER_SING)
		.value("LIKELY_OUTER", LoopType::LIKELY_OUTER)
		.value("LIKELY_INNER", LoopType::LIKELY_INNER)
		.value("UNCLEAR", LoopType::UNCLEAR)
		.value("ERROR", LoopType::ERROR);

	py::enum_<InferenceType>(m, "InferenceType")
		.value("CENTER", InferenceType::CENTER)
		.value("CENTROID", InferenceType::CENTROID)
		.value("MID_POINT", InferenceType::MID_POINT)
		.value("POINT", InferenceType::POINT)
		.value("TOP_AXIS_POINT", InferenceType::TOP_AXIS_POINT)
		.value("BOTTOM_AXIS_POINT", InferenceType::BOTTOM_AXIS_POINT)
		.value("MID_AXIS_POINT", InferenceType::MID_AXIS_POINT)
		.value("LOOP_CENTER", InferenceType::LOOP_CENTER);

	py::class_<MassProperties>(m, "MassProperties")
		.def_readonly("amount", &MassProperties::amount)
		.def_readonly("mass", &MassProperties::mass)
		.def_readonly("c_of_g", &MassProperties::c_of_g)
		.def_readonly("m_of_i", &MassProperties::m_of_i)
		.def_readonly("periphery", &MassProperties::periphery)
		.def("__repr__",
			[](const MassProperties& m) {
				return "MassProperties(amount=" + std::to_string(m.amount) + ")";
			});

	py::class_<Inference>(m, "Inference")
		.def_readonly("z_axis", &Inference::z_axis)
		.def_readonly("origin", &Inference::origin)
		.def_readonly("inference_type", &Inference::inference_type)
		.def_readonly("onshape_inference", &Inference::onshape_inference)
		.def_readonly("flipped_in_onshape", &Inference::flipped_in_onshape)
		.def("__repr__",
			[](const Inference& i) {
				return "<Inference>";
				/*
				return "Inference(origin=" + 
					"(" + std::to_string(i.origin(0)) + "," + std::to_string(i.origin(1)) + "," + std::to_string(i.origin(2)) + ")" +
					" , z_axis=" + 
					"(" + std::to_string(i.origin(0)) + "," + std::to_string(i.origin(1)) + "," + std::to_string(i.origin(2)) + ")" +
					")";
				*/
			});

	

	// topology.h
	py::enum_<TopoRelationSense>(m, "TopoRelationSense")
		.value("None", TopoRelationSense::None)
		.value("Negative", TopoRelationSense::Negative)
		.value("Positive", TopoRelationSense::Positive);

	py::class_<TopoRelation>(m, "TopoRelation")
		.def_readonly("_parent", &TopoRelation::_parent)
		.def_readonly("_child", &TopoRelation::_child)
		.def_readonly("_sense", &TopoRelation::_sense)
		.def("__repr__",
			[](const TopoRelation& tr) {
				return "TopoRelation(" + 
					std::to_string(tr._parent) + "," + 
					std::to_string(tr._child) + ")";
			});

	py::class_<BREPTopology>(m, "BREPTopology")
		.def_readonly("faces", &BREPTopology::faces)
		.def_readonly("loops", &BREPTopology::loops)
		.def_readonly("edges", &BREPTopology::edges)
		.def_readonly("vertices", &BREPTopology::vertices)
		.def_readonly("face_to_loop", &BREPTopology::face_to_loop)
		.def_readonly("loop_to_edge", &BREPTopology::loop_to_edge)
		.def_readonly("edge_to_vertex", &BREPTopology::edge_to_vertex)
		.def_readonly("loop_to_vertex", &BREPTopology::loop_to_vertex)
		.def_readonly("face_to_face", &BREPTopology::face_to_face)
		.def_readonly("face_loop", &BREPTopology::face_loop)
		.def_readonly("face_edge", &BREPTopology::face_edge)
		.def_readonly("face_vertex", &BREPTopology::face_vertex)
		.def_readonly("loop_edge", &BREPTopology::loop_edge)
		.def_readonly("loop_vertex", &BREPTopology::loop_vertex)
		.def_readonly("edge_vertex", &BREPTopology::edge_vertex)
		.def_readonly("pk_to_idx", &BREPTopology::pk_to_idx)
		.def_readonly("pk_to_class", &BREPTopology::pk_to_class)
		.def("__repr__",
			[](const BREPTopology& t) {
				return "<BREPTopology(" +
					std::to_string(t.faces.size()) + " faces, " +
					std::to_string(t.loops.size()) + " loops, " +
					std::to_string(t.edges.size()) + " edges, " +
					std::to_string(t.vertices.size()) + " vertices>";
			});


	// face.h
	py::class_<Face>(m, "Face")
		.def("get_inferences", &Face::get_inferences)
		.def("sample_points", &Face::sample_points)
		.def_readonly("function", &Face::function)
		.def_readonly("parameters", &Face::parameters)
		.def_readonly("orientation", &Face::orientation)
		.def_readonly("bounding_box", &Face::bounding_box)
		.def_readonly("na_bb_center", &Face::na_bb_center)
		.def_readonly("na_bb_x", &Face::na_bb_x)
		.def_readonly("na_bb_z", &Face::na_bb_z)
		.def_readonly("na_bounding_box", &Face::na_bounding_box)
		.def_readonly("surface_area", &Face::surface_area)
		.def_readonly("circumference", &Face::circumference)
		.def_readonly("center_of_gravity", &Face::center_of_gravity)
		.def_readonly("moment_of_inertia", &Face::moment_of_inertia)
		.def("__repr__",
			[](const Face& f) {
				std::string message = "";
				switch (f.function) {
				case SurfaceFunction::PLANE:
					message += "PLANE(";
					break;
				case SurfaceFunction::CYLINDER:
					message += "CYLINDER(";
					break;
				case SurfaceFunction::CONE:
					message += "CONE(";
					break;
				case SurfaceFunction::SPHERE:
					message += "SPHERE(";
					break;
				case SurfaceFunction::TORUS:
					message += "TORUS(";
					break;
				case SurfaceFunction::SPUN:
					message += "SPUN(";
					break;
				case SurfaceFunction::BSURF:
					message += "BSURF(";
					break;
				case SurfaceFunction::OFFSET:
					message += "OFFSET(";
					break;
				case SurfaceFunction::SWEPT:
					message += "SWEPT(";
					break;
				case SurfaceFunction::BLENDSF:
					message += "BLENDSF(";
					break;
				case SurfaceFunction::MESH:
					message += "MESH(";
					break;
				case SurfaceFunction::FSURF:
					message += "FSURF";
					break;
				}
				for (int i = 0; i < f.parameters.size(); ++i) {
					if (i > 0) {
						message += ",";
					}
					message += std::to_string(f.parameters[i]);
				}
				message += ")";
				return message;
			});


	// loop.h
		py::class_<Loop>(m, "Loop")
			.def("get_inferences", &Loop::get_inferences)
			.def_readonly("_type", &Loop::_type)
			.def_readonly("_is_circle", &Loop::_is_circle)
			.def_readonly("length", &Loop::length)
			.def_readonly("center_of_gravity", &Loop::center_of_gravity)
			.def_readonly("moment_of_inertia", &Loop::moment_of_inertia)
			.def_readonly("na_bb_center", &Loop::na_bb_center)
			.def_readonly("na_bb_x", &Loop::na_bb_x)
			.def_readonly("na_bb_z", &Loop::na_bb_z)
			.def_readonly("na_bounding_box", &Loop::na_bounding_box)
			.def("__repr__",
				[](const Loop& l) {
					std::string message = "";
					switch (l._type) {
					case LoopType::OUTER:
						message = "<Outer Loop>";
						break;
					case LoopType::LIKELY_OUTER:
						message = "<Likely Outer Loop>";
						break;
					case LoopType::INNER:
						message = "<Inner Loop>";
						break;
					case LoopType::LIKELY_INNER:
						message = "<Likely Inner Loop>";
						break;
					case LoopType::INNER_SING:
						message = "<Inner Singular Loop>";
						break;
					case LoopType::VERTEX:
						message = "<Vertex Loop>";
						break;
					case LoopType::UNCLEAR:
						message = "<Unclear Loop>";
						break;
					case LoopType::WINDING:
						message = "<Winding Loop>";
						break;
					case LoopType::WIRE:
						message = "<Wire Loop>";
						break;
					case LoopType::ERROR:
						message = "<Error Loop>";
						break;
					}
					return message;
				});

		// edge.h

		py::class_<Edge>(m, "Edge")
			.def("get_inferences", &Edge::get_inferences)
			.def("sample_points", &Edge::sample_points)
			.def_readonly("function", &Edge::function)
			.def_readonly("parameters", &Edge::parameters)
			.def_readonly("t_start", &Edge::t_start)
			.def_readonly("t_end", &Edge::t_end)
			.def_readonly("start", &Edge::start)
			.def_readonly("end", &Edge::end)
			.def_readonly("has_curve", &Edge::has_curve)
			.def_readonly("_is_reversed", &Edge::_is_reversed)
			.def_readonly("is_periodic", &Edge::is_periodic)
			.def_readonly("mid_point", &Edge::mid_point)
			.def_readonly("length", &Edge::length)
			.def_readonly("center_of_gravity", &Edge::center_of_gravity)
			.def_readonly("moment_of_inertia", &Edge::moment_of_inertia)
			.def_readonly("bounding_box", &Edge::bounding_box)
			.def_readonly("na_bb_center", &Edge::na_bb_center)
			.def_readonly("na_bb_x", &Edge::na_bb_x)
			.def_readonly("na_bb_z", &Edge::na_bb_z)
			.def_readonly("nn_bounding_box", &Edge::na_bounding_box)
			.def("__repr__",
				[](const Edge& e) {
					std::string message = "";

					switch (e.function) {
					case CurveFunction::LINE:
						message += "LINE(";
						break;
					case CurveFunction::CIRCLE:
						message += "CIRCLE(";
						break;
					case CurveFunction::ELLIPSE:
						message += "ELLIPSE(";
						break;
					case CurveFunction::BCURVE:
						message += "BCURVE(";
						break;
					case CurveFunction::ICURVE:
						message += "ICURVE(";
						break;
					case CurveFunction::FCURVE:
						message += "FCURVE(";
						break;
					case CurveFunction::SPCURVE:
						message += "SPCURVE(";
						break;
					case CurveFunction::TRCURVE:
						message += "TRCURVE(";
						break;
					case CurveFunction::CPCURVE:
						message += "CPCURVE(";
						break;
					case CurveFunction::PLINE:
						message += "PLINE(";
						break;
					}

					for (int i = 0; i < e.parameters.size(); ++i) {
						if (i > 0) {
							message += ",";
						}
						message += std::to_string(e.parameters[i]);
					}
					message += ")";
					return message;
				});
		

		// vertex.h

		py::class_<Vertex>(m, "Vertex")
			.def("get_inferences", &Vertex::get_inferences)
			.def_readonly("position", &Vertex::position)
			.def("__repr__",
				[](const Vertex& v) {
					return "Vertex(" + std::to_string(v.position(0)) + "," + std::to_string(v.position(1)) + "," + std::to_string(v.position(2)) + ")";
				});

}
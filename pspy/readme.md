# pspy

## A Python wrapper for loading Parasolid parts

This directory contains pspy, a Python wrapper that can
load parts in the Parasolid transmit (.x_t) format.

To run this, you must have a copy of Parasolid! Look at the
`readme.md` files in `./parasolid` and `./parasolid/lib` for
information on how to add your Parasolid distribution files
to this repository in order to build pspy.

Once you have all of the necessary files in place, pspy can
be installed with either

`python setup.py install`

or

`pip install .`

in this directory. Getting the compilation and linking working
across multiple platforms has been tricky - if you find that 
one of those commands does not work on your system, the other
may.

## Using pspy

The pspy module gives you one main entry point for loading a 
Parasolid text transmit file, the `pspart.Part` class. Its
constructor has one required argument of the path to the part file,
or its raw text

```
from pspy import Part

part_from_file = Part('/path/to/my/model.x_t')
part_from_text = Part('RAW TEXT OF .X_T FILE)

```

### Part features

Loaded parts have many feature types that can be accessed hierarchically.
An important property to check before using a part is `part.is_valid`. If
this is `False`, there was an error loading the part and it may contain
incorrect data.

 - Part
  - mesh
    - V
	- F
  - mesh_topology
    - face_to_topology
	- edge_to_topology
	- point_to_topology
  - brep
    - nodes
	  - faces (list)
		  - index
		  - function
		  - parameters
		  - orientation
		  - bounding_box
		  - na_bounding_box
		  - surface_area
		  - circumference
		  - center_of_gravity
		  - moment_of_inertia
		  - loop_neighbors
		  - edge_neighbors
		  - vertex_neighbors
		  - inferences
	  - loops (list)
		  - index
		  - type
		  - length
		  - center_of_gravity
		  - moment_of_inertia
		  - na_bounding_box
		  - edge_neighbors
		  - vertex_neighbors
		  - inferences
	  - edges (list)
	      - index
		  - function
		  - parameters
		  - t_range
		  - start
		  - end
		  - is_periodic
		  - mid_point
		  - length
		  - center_of_gravity
		  - moment_of_inertia
		  - na_bounding_box
		  - vertex_neighbors
		  - inferences
	  - vertices (list)
	      - index
		  - position
		  - inferences
	- relations
	  - face_to_loop
	  - loop_to_edge
	  - edge_to_vertex
	  - face_to_face
  - samples
    - face_samples
	- edge_samples
  - summary
    - bounding_box
	- volume
	- mass
	- center_of_gravity
	- moment_of_inertia
	- surface_area
	- topo_type_counts
	- surface_type_counts
	- loop_type_counts
	- fingerprint
  - inferences
    - origins
	- axes
	- frames
	- origin_references
	- axes_references
	- frame_references
  - default_mcfs
  - is_valid


Most of these features are fairly self-explanatory, but a few need clarification.
The `na_bounding_box` features are non-axis-aligned bounding boxes. These are
stored as (5,3) arrays, where the first 3 rows define the non-axis-aligned coordinate
system (center, x-direction, z-direction), and the last 2 rows are corners of a
bounding box within that coordinate system.

The `face_samples` and `edge_samples` are evenly spaced samples of a surface
or curve. In `face_samples` they are a list of 2d arrays [x, y, z, (n_x, n_y, n_z), mask], where all points are evaluated on a regular uv-grid. The normal directions are optional (see below). The `mask` grid indicates whether a point is within the clipping curve of a face. In `edge_samples`, a similar scheme is used but the arrays are 1-dimensional, and the order is [x,y,z, (t_x,t_y,t_z,pc_1,pc_2)]. In this case the tangent and principal curvatures are optional, and there is no mask since the curve is only evaluated within its clipped domain.

### Turning on and off features

The part class can also take an optional PartOptions
object which can turn on or off certain feature extraction steps
(such as mesh tesselation, or grid sampling of surfaces and curves).


```
from pspy import Part, PartOptions

opts = PartOptions()
opts.tesselate = false
opts.collect_inferences = false

part_with_options = Part(path, opts)
```

The options, and their default values, are as follows:

 - `just_bb = False` only compute the overall part bounding box
 - `normalize = False` normalize the part to a unit cube at the origin. Overrides `transform`
 - `transform = False` apply a transformation matrix to the part before loading.
 - `transform_matrix` transformation matrix to be applied if `transform == True`. No default value, must be set. 
 An invalid (non-rigid) transform will result in an invalid part.
 - `num_uv_samples = 10` number of grid samples per dimensios for surfaces and curves. Set to 0 to stop sampling.
 - `sample_normals = True` also sample grid normals for surfaces
 - `sample_tangents = True` also sample grid tangents and principle curvatures for curves
 - `tesselate = True` generate a mesh tesselation of the part and associations between mesh elements and topological entities
 - `default_mcfs = True` compute all potential MCFs as described in the Automate paper
 - `default_mcfs_only_face_axes = True` restrict potential MCFs to those oriented off of faces
 - `onshape_style = True` flip and exclude coordinate inferences in potential MCFs to match those from the Onshape CAD system
 - `collect_inferences` deduplication and collate coordinate inferences by origin, axis, and frame for easier searching
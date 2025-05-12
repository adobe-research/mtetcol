# MTetCol: A Mini Data Structure for Simplicial Columns

MTetCol is a C++ library that provides an efficient and minimal data structure for working with
simplicial columns, particularly useful for swept volume analysis and contour extraction.

## Features

- **Efficient Storage**: Compact representation of simplicial columns
- **Contouring**: Extract contours based on time derivative of swept volume function
- **Flexible Function Support**: 
  - Handle both cyclic and non-cyclic swept volume functions
  - Support for arbitrary time samples and function values
- **Triangulation**: Built-in support for cycle triangulation
- **Isocontour Extraction**: Extract zero-crossing isocontours from triangulated contours
- **File I/O**: Support for saving contours in MSH and OBJ formats

## Requirements

- C++20 or later
- CMake 3.28 or later
- A modern C++ compiler

## Installation

```bash
# Clone the repository
git clone https://github.com/qnzhou/mtetcol.git
cd mtetcol

# Create build directory
mkdir build && cd build

# Configure and build
cmake ..
make -j

# Install (optional)
make install
```

## Usage

### Basic Example

```c++
#include <mtetcol/simplicial_column.h>

// Initialize with tet mesh data
mtetcol::SimplicialColumn<4> columns;

// Set up the mesh
// tet_vertices: Flat array of tet vertices (3D coordinates)
// tet_mesh: Flat array of tet indices
columns.set_vertices(tet_vertices);
columns.set_simplices(tet_mesh);

// Set time samples and function values
// get_time_samples: Function returning time samples for a vertex
// get_function_values: Function returning function values for a vertex
columns.set_time_samples(get_time_samples, get_function_values);

// Extract contour at specific value
Scalar value = 0.0;
bool cyclic = false;
auto contour = columns.extract_contour(value, cyclic);

// The output contour is a 3D polyhedral mesh embedded in a 4D space-time domain.
// A cycle or polyhedron in the contour is considered "regular" if its spatial projection forms a
// non-degenerate simplex.
// 
// To check if a cycle/polyhedron is regular:
if (contour.is_cycle_regular(cycle_index)) { ... }
if (contour.is_polyhedron_regular(polyhedron_index)) { ... }

// The cycle/polyhedron regularity flag is computed only when extracting from simplicial columns.
// Subsequent operations, such as `triangulate_cycles` and `isocontour`, simply propagate these flags.

// Triangulate cycles if needed
auto triangulated_contours = columns.triangulate_cycles(contour);

// Extract isocontour
// Let `function_values` be a vector of function values at the contour vertices
auto isocontour = triangulated_contours.isocontour(function_values);
```

To save extracted contours:

```c++
#include <mtetcol/io.h>

mtetcol::io::save_contour("contour.msh", contour);  // MSH format
mtetcol::io::save_contour("contour.obj", contour);  // OBJ format
```


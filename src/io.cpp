#include <mtetcol/io.h>
#include <mtetcol/logger.h>

#include <mshio/mshio.h>

#include <fstream>
#include <stdexcept>

namespace mtetcol {

namespace {

mshio::MshSpec generate_spec(const Contour<4>& contour)
{
    size_t num_vertices = contour.get_num_vertices();
    size_t num_cycles = contour.get_num_cycles();

    mshio::MshSpec spec;
    spec.mesh_format.file_type = 1; // binary

    // Initialize nodes
    auto& nodes = spec.nodes;
    nodes.num_entity_blocks = 1;
    nodes.num_nodes = num_vertices;
    nodes.min_node_tag = 1;
    nodes.max_node_tag = nodes.num_nodes;
    nodes.entity_blocks.resize(1);

    auto& node_block = nodes.entity_blocks[0];
    node_block.entity_dim = 2;
    node_block.entity_tag = 1;
    node_block.parametric = 0;
    node_block.num_nodes_in_block = nodes.num_nodes;
    node_block.tags.reserve(nodes.num_nodes);
    node_block.data.reserve(nodes.num_nodes * 3);

    for (size_t i = 0; i < num_vertices; i++) {
        auto pos = contour.get_vertex(i);
        node_block.tags.push_back(i + 1);
        node_block.data.push_back(pos[0]);
        node_block.data.push_back(pos[1]);
        node_block.data.push_back(pos[2]);
    }

    // Initialize elements
    auto& elements = spec.elements;
    elements.num_entity_blocks = 1;
    elements.num_elements = num_cycles;
    elements.min_element_tag = 1;
    elements.max_element_tag = elements.num_elements;
    elements.entity_blocks.resize(1);

    auto& element_block = elements.entity_blocks[0];
    element_block.entity_dim = 2;
    element_block.entity_tag = 1;
    element_block.element_type = 2;
    element_block.num_elements_in_block = elements.num_elements;
    element_block.data.reserve(elements.num_elements * 4);

    for (size_t i = 0; i < num_cycles; i++) {
        auto cycle = contour.get_cycle(i);
        assert(cycle.size() == 3);
        element_block.data.push_back(i + 1);
        for (auto si : cycle) {
            Index seg_id = index(si);
            bool seg_ori = orientation(si);

            auto seg = contour.get_segment(seg_id);

            element_block.data.push_back((seg_ori ? seg[0] : seg[1]) + 1);
        }
    }

    // Initiailize time attribute
    mshio::Data node_data;
    node_data.header.string_tags.push_back("time");
    node_data.header.int_tags.push_back(0);
    node_data.header.int_tags.push_back(1);
    node_data.header.int_tags.push_back(num_vertices);

    auto& entries = node_data.entries;
    entries.resize(num_vertices);
    for (size_t i = 0; i < num_vertices; i++) {
        auto pos = contour.get_vertex(i);
        auto& entry = entries[i];
        entry.tag = i + 1;
        entry.data = {pos[3]};
    }

    spec.node_data.push_back(std::move(node_data));

    return spec;
}

mshio::MshSpec generate_spec(const Contour<3>& contour)
{
    size_t num_vertices = contour.get_num_vertices();
    size_t num_segments = contour.get_num_segments();

    mshio::MshSpec spec;
    spec.mesh_format.file_type = 1; // binary

    // Initialize nodes
    auto& nodes = spec.nodes;
    nodes.num_entity_blocks = 1;
    nodes.num_nodes = num_vertices;
    nodes.min_node_tag = 1;
    nodes.max_node_tag = nodes.num_nodes;
    nodes.entity_blocks.resize(1);

    auto& node_block = nodes.entity_blocks[0];
    node_block.entity_dim = 1;
    node_block.entity_tag = 1;
    node_block.parametric = 0;
    node_block.num_nodes_in_block = nodes.num_nodes;
    node_block.tags.reserve(nodes.num_nodes);
    node_block.data.reserve(nodes.num_nodes * 3);

    for (size_t i = 0; i < num_vertices; i++) {
        auto pos = contour.get_vertex(i);
        node_block.tags.push_back(i + 1);
        node_block.data.push_back(pos[0]);
        node_block.data.push_back(pos[1]);
        node_block.data.push_back(0);
    }

    // Initialize elements
    auto& elements = spec.elements;
    elements.num_entity_blocks = 1;
    elements.num_elements = num_segments;
    elements.min_element_tag = 1;
    elements.max_element_tag = elements.num_elements;
    elements.entity_blocks.resize(1);

    auto& element_block = elements.entity_blocks[0];
    element_block.entity_dim = 1;
    element_block.entity_tag = 1;
    element_block.element_type = 1;
    element_block.num_elements_in_block = elements.num_elements;
    element_block.data.reserve(elements.num_elements * 4);

    for (size_t i = 0; i < num_segments; i++) {
        auto seg = contour.get_segment(i);
        element_block.data.push_back(i + 1);
        element_block.data.push_back(seg[0] + 1);
        element_block.data.push_back(seg[1] + 1);
    }

    // Initiailize time attribute
    mshio::Data node_data;
    node_data.header.string_tags.push_back("time");
    node_data.header.int_tags.push_back(0);
    node_data.header.int_tags.push_back(1);
    node_data.header.int_tags.push_back(num_vertices);

    auto& entries = node_data.entries;
    entries.resize(num_vertices);
    for (size_t i = 0; i < num_vertices; i++) {
        auto pos = contour.get_vertex(i);
        auto& entry = entries[i];
        entry.tag = i + 1;
        entry.data = {pos[2]};
    }

    spec.node_data.push_back(std::move(node_data));

    return spec;
}


void add_scalar_field(mshio::MshSpec& spec, std::string name, std::span<Scalar> field)
{
    mshio::Data node_data;
    node_data.header.string_tags.push_back(name);
    node_data.header.int_tags.push_back(0);
    node_data.header.int_tags.push_back(1);
    node_data.header.int_tags.push_back(field.size());

    auto& entries = node_data.entries;
    entries.resize(field.size());
    for (size_t i = 0; i < field.size(); i++) {
        auto& entry = entries[i];
        entry.tag = i + 1;
        entry.data = {field[i]};
    }

    spec.node_data.push_back(std::move(node_data));
}

} // namespace

void save_contour(
    std::string_view filename,
    const Contour<3>& contour,
    std::span<Scalar> function_values)
{
    logger().info("Saving contour to {}", filename);
    if (filename.ends_with(".msh")) {
        auto spec = generate_spec(contour);
        if (!function_values.empty()) {
            assert(function_values.size() == contour.get_num_vertices());
            add_scalar_field(spec, "values", function_values);
        }
        mshio::save_msh(std::string(filename), spec);
    } else {
        std::ofstream fout(filename.data());

        size_t num_vertices = contour.get_num_vertices();
        size_t num_segments = contour.get_num_segments();

        for (size_t i = 0; i < num_vertices; i++) {
            auto pos = contour.get_vertex(i);
            fout << "v " << pos[0] << " " << pos[1] << " 0" << std::endl;
        }

        for (size_t i = 0; i < num_segments; i++) {
            auto seg = contour.get_segment(i);
            fout << "l " << seg[0] + 1 << " " << seg[1] + 1 << std::endl;
        }
    }
}

void save_contour(
    std::string_view filename,
    const Contour<4>& contour,
    std::span<Scalar> function_values)
{
    logger().info("Saving contour to {}", filename);
    if (filename.ends_with(".msh")) {
        auto spec = generate_spec(contour);
        if (!function_values.empty()) {
            assert(function_values.size() == contour.get_num_vertices());
            add_scalar_field(spec, "values", function_values);
        }
        mshio::save_msh(std::string(filename), spec);
    } else {
        std::ofstream fout(filename.data());

        size_t num_vertices = contour.get_num_vertices();
        size_t num_cycles = contour.get_num_cycles();

        for (size_t i = 0; i < num_vertices; i++) {
            auto pos = contour.get_vertex(i);
            fout << "v " << pos[0] << " " << pos[1] << " " << pos[2] << " " << std::endl;
        }

        for (size_t i = 0; i < num_cycles; i++) {
            auto cycle = contour.get_cycle(i);
            fout << "f ";
            for (auto si : cycle) {
                Index seg_id = index(si);
                bool seg_ori = orientation(si);

                auto seg = contour.get_segment(seg_id);

                fout << (seg_ori ? seg[0] : seg[1]) + 1 << " ";
            }
            fout << std::endl;
        }
    }
}

} // namespace mtetcol

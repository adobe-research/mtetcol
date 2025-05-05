#include <mtetcol/io.h>
#include <mtetcol/logger.h>

#include <fstream>
#include <stdexcept>

namespace mtetcol {

void save_contour(std::string_view filename, const Contour<3>& contour)
{
    assert(false);
    throw std::runtime_error("Saving 3D contour is not implemented yet");
    // Implement saving logic for 3D contour
}

void save_contour(std::string_view filename, const Contour<4>& contour)
{
    logger().info("Saving contour to {}", filename);
    std::ofstream fout(filename.data());

    size_t num_vertices = contour.get_num_vertices();
    size_t num_cycles = contour.get_num_cycles();

    for (size_t i=0; i<num_vertices; i++) {
        auto pos = contour.get_vertex(i);
        fout << "v " << pos[0] << " " << pos[1] << " " << pos[2] << " " << std::endl;
    }

    for (size_t i=0; i<num_cycles; i++) {
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

} // namespace mtetcol

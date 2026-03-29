#ifndef SOFTBODY_SOLVER_H
#define SOFTBODY_SOLVER_H

#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/classes/mesh_instance3d.hpp>
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/variant/packed_vector3_array.hpp>
#include <godot_cpp/variant/packed_int32_array.hpp>
#include <godot_cpp/variant/packed_float32_array.hpp>

namespace godot {

class SoftBodySolver : public Node3D {
    GDCLASS(SoftBodySolver, Node3D)

private:
    // Physics data
    PackedVector3Array pos;
    PackedVector3Array prev_pos;
    PackedVector3Array velocity;
    PackedFloat32Array inv_mass;
    
    // Mesh data
    PackedInt32Array edge_ids;
    PackedInt32Array tet_ids;
    PackedFloat32Array edge_lengths;
    PackedFloat32Array tet_volumes;
    
    // Parameters
    double edge_compliance = 0.0;
    double volume_compliance = 0.0;
    
    // Helper arrays for optimization
    Vector3 *pos_ptr;
    Vector3 *prev_pos_ptr;
    Vector3 *velocity_ptr;
    float *inv_mass_ptr;
    int32_t *edge_ids_ptr;
    int32_t *tet_ids_ptr;
    float *edge_lengths_ptr;
    float *tet_volumes_ptr;
    
    int num_particles;
    int num_edges;
    int num_tets;
    
    // Volume gradient order (cache as static array)
    static const int vol_id_order[4][3];

protected:
    static void _bind_methods();

public:
    SoftBodySolver();
    ~SoftBodySolver();
    
    // Core solver methods
    void solve_edges(double compliance, double dt);
    void solve_volumes(double compliance, double dt);
    double get_tet_volume(int tet_index);
    
    // Initialization
    void compute_edge_rest_lengths();
    void compute_tet_rest_volumes();
    void pre_solve(double dt, Vector3 force);
    void post_solve(double dt);
    void solve(double dt);
    
    // Getters/Setters
    void set_pos(const PackedVector3Array &p_pos);
    PackedVector3Array get_pos() const;
    
    void set_edge_ids(const PackedInt32Array &p_ids);
    void set_tet_ids(const PackedInt32Array &p_ids);
    
    void set_edge_compliance(double p_compliance);
    void set_volume_compliance(double p_compliance);
    
    // Optimization: process multiple tetrahedra with SIMD
    void solve_volumes_simd(double compliance, double dt);
};

} // namespace godot

#endif
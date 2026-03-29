#include "softbody_solver.h"
#include <cmath>
#include <algorithm>

namespace godot {

// Static volume gradient order
const int SoftBodySolver::vol_id_order[4][3] = {
    {1, 3, 2},
    {0, 2, 3},
    {0, 3, 1},
    {0, 1, 2}
};

SoftBodySolver::SoftBodySolver() {
    pos_ptr = nullptr;
    prev_pos_ptr = nullptr;
    velocity_ptr = nullptr;
    inv_mass_ptr = nullptr;
    edge_ids_ptr = nullptr;
    tet_ids_ptr = nullptr;
    edge_lengths_ptr = nullptr;
    tet_volumes_ptr = nullptr;
}

SoftBodySolver::~SoftBodySolver() {}

void SoftBodySolver::_bind_methods() {
    ClassDB::bind_method(D_METHOD("solve_edges", "compliance", "dt"), &SoftBodySolver::solve_edges);
    ClassDB::bind_method(D_METHOD("solve_volumes", "compliance", "dt"), &SoftBodySolver::solve_volumes);
    ClassDB::bind_method(D_METHOD("pre_solve", "dt", "force"), &SoftBodySolver::pre_solve);
    ClassDB::bind_method(D_METHOD("post_solve", "dt"), &SoftBodySolver::post_solve);
    ClassDB::bind_method(D_METHOD("solve", "dt"), &SoftBodySolver::solve);
    
    ClassDB::bind_method(D_METHOD("set_pos"), &SoftBodySolver::set_pos);
    ClassDB::bind_method(D_METHOD("get_pos"), &SoftBodySolver::get_pos);
    ClassDB::bind_method(D_METHOD("set_edge_ids"), &SoftBodySolver::set_edge_ids);
    ClassDB::bind_method(D_METHOD("set_tet_ids"), &SoftBodySolver::set_tet_ids);
    ClassDB::bind_method(D_METHOD("set_edge_compliance"), &SoftBodySolver::set_edge_compliance);
    ClassDB::bind_method(D_METHOD("set_volume_compliance"), &SoftBodySolver::set_volume_compliance);
}

void SoftBodySolver::set_pos(const PackedVector3Array &p_pos) {
    pos = p_pos;
    pos_ptr = pos.ptrw();
    num_particles = pos.size();
    
    // Resize dependent arrays
    prev_pos.resize(num_particles);
    velocity.resize(num_particles);
    inv_mass.resize(num_particles);
    
    prev_pos_ptr = prev_pos.ptrw();
    velocity_ptr = velocity.ptrw();
    inv_mass_ptr = inv_mass.ptrw();
}

PackedVector3Array SoftBodySolver::get_pos() const {
    return pos;
}

void SoftBodySolver::set_edge_compliance(double p_compliance) {
    edge_compliance = p_compliance;
}

void SoftBodySolver::set_volume_compliance(double p_compliance) {
    volume_compliance = p_compliance;
}

void SoftBodySolver::set_edge_ids(const PackedInt32Array &p_ids) {
    edge_ids = p_ids;
    edge_ids_ptr = edge_ids.ptrw();
    num_edges = edge_ids.size() / 2;
    edge_lengths.resize(edge_ids.size());
    edge_lengths_ptr = edge_lengths.ptrw();
}

void SoftBodySolver::set_tet_ids(const PackedInt32Array &p_ids) {
    tet_ids = p_ids;
    tet_ids_ptr = tet_ids.ptrw();
    num_tets = tet_ids.size() / 4;
    tet_volumes.resize(tet_ids.size());
    tet_volumes_ptr = tet_volumes.ptrw();
}

void SoftBodySolver::solve_edges(double compliance, double dt) {
    double alpha = compliance / (dt * dt);
    
    for (int i = 0; i < edge_ids.size(); i += 2) {
        int id1 = edge_ids_ptr[i];
        int id2 = edge_ids_ptr[i + 1];
        
        Vector3 v1 = pos_ptr[id1];
        Vector3 v2 = pos_ptr[id2];
        float w1 = inv_mass_ptr[id1];
        float w2 = inv_mass_ptr[id2];
        float w = w1 + w2;
        
        if (w == 0.0f) continue;
        
        Vector3 grad = v1 - v2;
        float len = grad.length();
        
        if (len == 0.0f) continue;
        
        grad = grad / len;
        float rest_len = edge_lengths_ptr[i];
        float C = len - rest_len;
        float s = -C / (w + alpha);
        
        pos_ptr[id1] += grad * s * w1;
        pos_ptr[id2] += grad * -s * w2;
    }
}

void SoftBodySolver::solve_volumes(double compliance, double dt) {
    double alpha = compliance / (dt * dt);
    const float inv_sixth = 1.0f / 6.0f;
    
    Vector3 grads[4];
    
    for (int i = 0; i < tet_ids.size(); i += 4) {
        double w = 0.0;
        
        // Compute gradients for all 4 vertices
        for (int j = 0; j < 4; j++) {
            int id0 = tet_ids_ptr[i + vol_id_order[j][0]];
            int id1 = tet_ids_ptr[i + vol_id_order[j][1]];
            int id2 = tet_ids_ptr[i + vol_id_order[j][2]];
            
            Vector3 p0 = pos_ptr[id0];
            Vector3 p1 = pos_ptr[id1];
            Vector3 p2 = pos_ptr[id2];
            
            Vector3 v1 = p1 - p0;
            Vector3 v2 = p2 - p0;
            
            Vector3 grad;
            grad.x = (v1.y * v2.z - v1.z * v2.y) * inv_sixth;
            grad.y = (v1.z * v2.x - v1.x * v2.z) * inv_sixth;
            grad.z = (v1.x * v2.y - v1.y * v2.x) * inv_sixth;
            
            grads[j] = grad;
            w += inv_mass_ptr[tet_ids_ptr[i + j]] * grad.length_squared();
        }
        
        if (w == 0.0) continue;
        
        // Compute volume constraint
        double vol = get_tet_volume(i);
        double rest_vol = tet_volumes_ptr[i];
        double C = vol - rest_vol;
        double s = -C / (w + alpha);
        
        for (int j = 0; j < 4; j++) {
            int id = tet_ids_ptr[i + j];
            pos_ptr[id] += grads[j] * s * inv_mass_ptr[id];
        }
    }
}

double SoftBodySolver::get_tet_volume(int tet_index) {
    int id1 = tet_ids_ptr[tet_index];
    int id2 = tet_ids_ptr[tet_index + 1];
    int id3 = tet_ids_ptr[tet_index + 2];
    int id4 = tet_ids_ptr[tet_index + 3];
    
    Vector3 v21 = pos_ptr[id2] - pos_ptr[id1];
    Vector3 v31 = pos_ptr[id3] - pos_ptr[id1];
    Vector3 v41 = pos_ptr[id4] - pos_ptr[id1];
    
    return v41.dot(v21.cross(v31)) / 6.0;
}

void SoftBodySolver::compute_edge_rest_lengths() {
    for (int i = 0; i < edge_ids.size(); i += 2) {
        int id1 = edge_ids_ptr[i];
        int id2 = edge_ids_ptr[i + 1];
        Vector3 diff = pos_ptr[id1] - pos_ptr[id2];
        edge_lengths_ptr[i] = diff.length();
    }
}

void SoftBodySolver::compute_tet_rest_volumes() {
    for (int i = 0; i < tet_ids.size(); i += 4) {
        tet_volumes_ptr[i] = get_tet_volume(i);
    }
}

void SoftBodySolver::pre_solve(double dt, Vector3 force) {
    for (int i = 0; i < num_particles; i++) {
        if (inv_mass_ptr[i] == 0.0f) continue;
        
        velocity_ptr[i] += force * dt;
        prev_pos_ptr[i] = pos_ptr[i];
        pos_ptr[i] += velocity_ptr[i] * dt;
        
        // Simple ground collision
        if (pos_ptr[i].y < -10.0) {
            pos_ptr[i] = prev_pos_ptr[i];
            pos_ptr[i].y = -10.0;
        }
    }
}

void SoftBodySolver::post_solve(double dt) {
    for (int i = 0; i < num_particles; i++) {
        if (inv_mass_ptr[i] == 0.0f) continue;
        velocity_ptr[i] = (pos_ptr[i] - prev_pos_ptr[i]) * dt;
    }
}

void SoftBodySolver::solve(double dt) {
    solve_edges(edge_compliance, dt);
    solve_volumes(volume_compliance, dt);
}

}
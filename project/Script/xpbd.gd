extends Node

var solver: SoftBodySolver
var tet_mesh: TetGenMesh
var tet_array_constructor: TetArrayConstructor
var mesh_instance: MeshInstance3D

func _ready():
	tet_mesh = TetGenMesh.new()
	tet_mesh.load_from_base_name("res://Mesh/LowPoly/Suzanne")
	
	tet_array_constructor = TetArrayConstructor.new()
	var edge_ids: PackedInt32Array = tet_array_constructor.construct_edge(tet_mesh.edges)
	var tet_ids: PackedInt32Array = tet_array_constructor.construct_tet(tet_mesh.tetrahedra)
	
	solver = SoftBodySolver.new()
	add_child(solver)
	
	solver.set_pos(tet_mesh.vertices)
	solver.set_edge_ids(edge_ids)
	solver.set_tet_ids(tet_ids)
	
	var inv_mass: PackedFloat32Array = PackedFloat32Array()
	inv_mass.resize(solver.get_pos().size())
	inv_mass.fill(0.5)
	solver.set_inv_mass(inv_mass)
	
	var surface_mesh = tet_mesh.create_surface_mesh()
	mesh_instance = MeshInstance3D.new()
	mesh_instance.mesh = surface_mesh
	add_child(mesh_instance)
	
	solver.compute_edge_rest_lengths()
	solver.compute_tet_rest_volumes()
	solver.set_edge_compliance(0.0)
	solver.set_volume_compliance(0.0)

func _process(delta: float):
	var force: Vector3 = Vector3(0, -5, 0)
	var dt: float = 0.05
	
	var substeps: float = 10
	var sub_dt: float = dt / substeps
	
	for step in range(substeps):
		solver.pre_solve(sub_dt, force)
		solver.solve(sub_dt)
		solver.post_solve(sub_dt)
	
	var updated_mesh: Mesh = tet_mesh.update_mesh(solver.get_pos())
	mesh_instance.mesh = updated_mesh

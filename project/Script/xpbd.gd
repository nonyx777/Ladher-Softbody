extends Node

var bunny_solver: SoftBodySolver
var bunny_tet_mesh: TetGenMesh
var bunny_tet_array_constructor: TetArrayConstructor
var bunny_mesh_instance: MeshInstance3D

var suzanne_solver: SoftBodySolver
var suzanne_tet_mesh: TetGenMesh
var suzanne_tet_array_constructor: TetArrayConstructor
var suzanne_mesh_instance: MeshInstance3D

var push_force: Vector3 = Vector3(0, 0, 0)

func bunny_init() -> void:
	bunny_tet_mesh = TetGenMesh.new()
	bunny_tet_mesh.load_from_base_name("res://Mesh/StanfordBunny/StanfordBunny")
	
	bunny_tet_array_constructor = TetArrayConstructor.new()
	var bunny_edge_ids: PackedInt32Array = bunny_tet_array_constructor.construct_edge(bunny_tet_mesh.edges)
	var bunny_tet_ids: PackedInt32Array = bunny_tet_array_constructor.construct_tet(bunny_tet_mesh.tetrahedra)
	
	bunny_solver = SoftBodySolver.new()
	add_child(bunny_solver)
	
	bunny_solver.set_pos(bunny_tet_mesh.vertices)
	bunny_solver.set_edge_ids(bunny_edge_ids)
	bunny_solver.set_tet_ids(bunny_tet_ids)
	
	var bunny_inv_mass: PackedFloat32Array = PackedFloat32Array()
	bunny_inv_mass.resize(bunny_solver.get_pos().size())
	bunny_inv_mass.fill(0.5)
	bunny_solver.set_inv_mass(bunny_inv_mass)
	
	var bunny_surface_mesh = bunny_tet_mesh.create_surface_mesh()
	bunny_mesh_instance = MeshInstance3D.new()
	bunny_mesh_instance.mesh = bunny_surface_mesh
	add_child(bunny_mesh_instance)
	
	bunny_solver.compute_edge_rest_lengths()
	bunny_solver.compute_tet_rest_volumes()
	bunny_solver.set_edge_compliance(0.0)
	bunny_solver.set_volume_compliance(0.0)

func suzanne_init() -> void:
	suzanne_tet_mesh = TetGenMesh.new()
	suzanne_tet_mesh.load_from_base_name("res://Mesh/Suzanne/LowPoly/Suzanne")
	
	suzanne_tet_array_constructor = TetArrayConstructor.new()
	var suzanne_edge_ids: PackedInt32Array = suzanne_tet_array_constructor.construct_edge(suzanne_tet_mesh.edges)
	var suzanne_tet_ids: PackedInt32Array = suzanne_tet_array_constructor.construct_tet(suzanne_tet_mesh.tetrahedra)
	
	suzanne_solver = SoftBodySolver.new()
	add_child(suzanne_solver)
	
	suzanne_solver.set_pos(suzanne_tet_mesh.vertices)
	suzanne_solver.set_edge_ids(suzanne_edge_ids)
	suzanne_solver.set_tet_ids(suzanne_tet_ids)
	
	var suzanne_inv_mass: PackedFloat32Array = PackedFloat32Array()
	suzanne_inv_mass.resize(suzanne_solver.get_pos().size())
	suzanne_inv_mass.fill(0.5)
	suzanne_solver.set_inv_mass(suzanne_inv_mass)
	
	var suzanne_surface_mesh = suzanne_tet_mesh.create_surface_mesh()
	suzanne_mesh_instance = MeshInstance3D.new()
	suzanne_mesh_instance.mesh = suzanne_surface_mesh
	add_child(suzanne_mesh_instance)
	
	suzanne_solver.compute_edge_rest_lengths()
	suzanne_solver.compute_tet_rest_volumes()
	suzanne_solver.set_edge_compliance(0.0)
	suzanne_solver.set_volume_compliance(0.0)

func _ready():
	bunny_init()
	suzanne_init()

func _process(delta: float):
	push_force *= 0.0
	if Input.is_key_pressed(KEY_S):
		push_force += Vector3(2, 2, 0)
	if Input.is_key_pressed(KEY_W):
		push_force += Vector3(-2, 2, 0)
	if Input.is_key_pressed(KEY_D):
		push_force += Vector3(0, 2, -2)
	if Input.is_key_pressed(KEY_A):
		push_force += Vector3(0, 2, 2)
	if Input.is_key_pressed(KEY_SPACE):
		push_force += Vector3(0, 10, 0)
	
	var force: Vector3 = Vector3(0, -5, 0) + push_force
	var dt: float = 0.04
	
	var substeps: float = 10
	var sub_dt: float = dt / substeps
	
	for step in range(substeps):
		bunny_solver.pre_solve(sub_dt, force)
		bunny_solver.solve(sub_dt)
		bunny_solver.post_solve(sub_dt)
		
		suzanne_solver.pre_solve(sub_dt, force)
		suzanne_solver.solve(sub_dt)
		suzanne_solver.post_solve(sub_dt)
	
	var bunny_updated_mesh: Mesh = bunny_tet_mesh.update_mesh(bunny_solver.get_pos())
	bunny_mesh_instance.mesh = bunny_updated_mesh
	
	var suzanne_updated_mesh: Mesh = suzanne_tet_mesh.update_mesh(suzanne_solver.get_pos())
	suzanne_mesh_instance.mesh = suzanne_updated_mesh

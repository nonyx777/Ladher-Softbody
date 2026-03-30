class_name  TetGenMesh

var vertices: PackedVector3Array = PackedVector3Array()
var tetrahedra: Array = []
var surface_faces: Array = []
var edges: Array = []
var face_markers: Array = []  # Store boundary markers for faces

func load_from_base_name(base_path: String):
	load_nodes(base_path + ".node")
	load_ele(base_path + ".ele")
	load_face(base_path + ".face")
	load_edge(base_path + ".edge")

func load_nodes(filepath: String):
	vertices.clear()
	var file = FileAccess.open(filepath, FileAccess.READ)
	if not file:
		print("Error: Cannot open ", filepath)
		return
	
	# Read first line
	var line = file.get_line()
	print("Node file header: ", line)
	
	var parts = line.strip_edges().split(" ", false)
	if parts.size() < 2:
		print("Invalid .node file format")
		return
	
	var num_vertices = int(parts[0])
	var dimension = int(parts[1]) if parts.size() > 1 else 3
	var num_attributes = int(parts[2]) if parts.size() > 2 else 0
	var has_markers = int(parts[3]) if parts.size() > 3 else 0
	
	print("Loading ", num_vertices, " vertices (dim:", dimension, ", attr:", num_attributes, ", markers:", has_markers, ")")
	
	# Read each vertex
	for i in range(num_vertices):
		line = file.get_line()
		while line.strip_edges() == "":
			line = file.get_line()
		
		parts = line.strip_edges().split(" ", false)
		if parts.size() < 1 + dimension:
			print("Error parsing vertex at line ", i+2)
			continue
		
		# Index is first token, skip it
		var x = float(parts[1])
		var y = float(parts[2]) if dimension >= 2 else 0.0
		var z = float(parts[3]) if dimension >= 3 else 0.0
		
		vertices.append(Vector3(x, y, z))
	
	file.close()
	#print("Loaded ", vertices.size(), " vertices from ", filepath)
	print("Vertices: ", vertices.slice(0, 10))

func load_ele(filepath: String):
	tetrahedra.clear()
	var file = FileAccess.open(filepath, FileAccess.READ)
	if not file:
		print("Error: Cannot open ", filepath)
		return
	
	var line = file.get_line()
	print("Ele file header: ", line)
	
	var parts = line.strip_edges().split(" ", false)
	if parts.size() < 2:
		print("Invalid .ele file format")
		return
	
	var num_tets = int(parts[0])
	var nodes_per_tet = int(parts[1])
	var num_attributes = int(parts[2]) if parts.size() > 2 else 0
	
	print("Loading ", num_tets, " tetrahedra (nodes:", nodes_per_tet, ", attr:", num_attributes, ")")
	
	for i in range(num_tets):
		line = file.get_line()
		while line.strip_edges() == "":
			line = file.get_line()
		
		parts = line.strip_edges().split(" ", false)
		if parts.size() < 1 + nodes_per_tet:
			print("Error parsing tet at line ", i+2)
			continue
		
		# Convert from 1-based to 0-based indices
		var v0 = int(parts[1]) - 1
		var v1 = int(parts[2]) - 1
		var v2 = int(parts[3]) - 1
		var v3 = int(parts[4]) - 1
		
		tetrahedra.append([v0, v1, v2, v3])
	
	file.close()
	#print("Loaded ", tetrahedra.size(), " tetrahedra from ", filepath)
	print("Tetrahedra: ", tetrahedra.slice(0, 10))

func load_face(filepath: String):
	surface_faces.clear()
	face_markers.clear()
	var file = FileAccess.open(filepath, FileAccess.READ)
	if not file:
		print("Error: Cannot open ", filepath)
		return
	
	var line = file.get_line()
	print("Face file header: ", line)
	
	var parts = line.strip_edges().split(" ", false)
	if parts.size() < 1:
		print("Invalid .face file format")
		return
	
	var num_faces = int(parts[0])
	var has_markers = int(parts[1]) if parts.size() > 1 else 0
	
	print("Loading ", num_faces, " faces (markers:", has_markers, ")")
	
	for i in range(num_faces):
		line = file.get_line()
		while line.strip_edges() == "":
			line = file.get_line()
		
		parts = line.strip_edges().split(" ", false)
		if parts.size() < 4:
			print("Error parsing face at line ", i+2)
			continue
		
		# Convert from 1-based to 0-based indices
		var v0 = int(parts[1]) - 1
		var v1 = int(parts[2]) - 1
		var v2 = int(parts[3]) - 1
		
		surface_faces.append([v0, v1, v2])
		
		# Store marker if present
		if has_markers > 0 and parts.size() > 4:
			face_markers.append(int(parts[4]))
		else:
			face_markers.append(0)
	
	file.close()
	#print("Loaded ", surface_faces.size(), " surface faces from ", filepath)
	print("Surfaces Faces: ", surface_faces.slice(0, 10))

func load_edge(filepath: String):
	edges.clear()
	var file = FileAccess.open(filepath, FileAccess.READ)
	if not file:
		print("Error: Cannot open ", filepath)
		return
	
	var line = file.get_line()
	print("Edge file header: ", line)
	
	var parts = line.strip_edges().split(" ", false)
	if parts.size() < 1:
		print("Invalid .edge file format")
		return
	
	var num_edges = int(parts[0])
	var has_markers = int(parts[1]) if parts.size() > 1 else 0
	
	print("Loading ", num_edges, " edges (markers:", has_markers, ")")
	
	for i in range(num_edges):
		line = file.get_line()
		while line.strip_edges() == "":
			line = file.get_line()
		
		parts = line.strip_edges().split(" ", false)
		if parts.size() < 3:
			print("Error parsing edge at line ", i+2)
			continue
		
		# Convert from 1-based to 0-based indices
		var v0 = int(parts[1]) - 1
		var v1 = int(parts[2]) - 1
		
		edges.append([v0, v1])
	
	file.close()
	#print("Loaded ", edges.size(), " edges from ", filepath)
	print("Edges: ", edges.slice(0, 10))

func create_surface_mesh() -> ArrayMesh:
	if surface_faces.is_empty():
		print("No surface faces to create mesh")
		return null
	
	var arrays = []
	arrays.resize(Mesh.ARRAY_MAX)
	
	# Vertices
	arrays[Mesh.ARRAY_VERTEX] = vertices
	
	# Indices (triangles)
	var indices = PackedInt32Array()
	for face in surface_faces:
		indices.append(face[0])
		indices.append(face[1])
		indices.append(face[2])
	arrays[Mesh.ARRAY_INDEX] = indices
	
	# Calculate normals
	var normals = PackedVector3Array()
	normals.resize(vertices.size())
	for i in range(normals.size()):
		normals[i] = Vector3.ZERO
	
	# First pass: accumulate normals
	for i in range(0, indices.size(), 3):
		var i0 = indices[i]
		var i1 = indices[i+1]
		var i2 = indices[i+2]
		
		var v0 = vertices[i0]
		var v1 = vertices[i1]
		var v2 = vertices[i2]
		
		var normal = (v1 - v0).cross(v2 - v0)
		
		if normal.length() < 0.0001:
			continue
			
		normal = -normal.normalized()
		
		normals[i0] += normal
		normals[i1] += normal
		normals[i2] += normal
	
	# Normalize normals
	for i in range(normals.size()):
		normals[i] = normals[i].normalized()
	
	arrays[Mesh.ARRAY_NORMAL] = normals
	
	# Create mesh
	var mesh = ArrayMesh.new()
	mesh.add_surface_from_arrays(Mesh.PRIMITIVE_TRIANGLES, arrays)
	
	print("Created surface mesh with ", vertices.size(), " vertices and ", surface_faces.size(), " triangles")
	print("Mesh bounds: ", mesh.get_aabb())
	
	return mesh

func update_mesh(vert: PackedVector3Array) -> ArrayMesh:
	var arrays = []
	arrays.resize(Mesh.ARRAY_MAX)
	
	# Vertices
	arrays[Mesh.ARRAY_VERTEX] = vert
	
	# Indices (triangles)
	var indices = PackedInt32Array()
	for face in surface_faces:
		indices.append(face[0])
		indices.append(face[1])
		indices.append(face[2])
	arrays[Mesh.ARRAY_INDEX] = indices
	
	# Calculate normals
	var normals = PackedVector3Array()
	normals.resize(vertices.size())
	for i in range(normals.size()):
		normals[i] = Vector3.ZERO
	
	# First pass: accumulate normals
	for i in range(0, indices.size(), 3):
		var i0 = indices[i]
		var i1 = indices[i+1]
		var i2 = indices[i+2]
		
		var v0 = vertices[i0]
		var v1 = vertices[i1]
		var v2 = vertices[i2]
		
		var normal = (v1 - v0).cross(v2 - v0)
		
		if normal.length() < 0.0001:
			continue
			
		normal = -normal.normalized()
		
		normals[i0] += normal
		normals[i1] += normal
		normals[i2] += normal
	
	# Normalize normals
	for i in range(normals.size()):
		normals[i] = normals[i].normalized()
	
	arrays[Mesh.ARRAY_NORMAL] = normals
	
	# Create mesh
	var mesh = ArrayMesh.new()
	mesh.add_surface_from_arrays(Mesh.PRIMITIVE_TRIANGLES, arrays)
	
	return mesh

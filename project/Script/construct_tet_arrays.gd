class_name TetArrayConstructor

func construct_compact_edge(edges: Array, edge_stride: PackedInt32Array, edge_vertex_indices: PackedInt32Array) -> PackedInt32Array:
	# Build adjacency list
	var adjacency = {}
	
	# Add all edges to adjacency list
	for edge in edges:
		var v0 = edge[0]
		var v1 = edge[1]
		
		# Add bidirectional connections
		if not adjacency.has(v0):
			adjacency[v0] = []
		if not adjacency.has(v1):
			adjacency[v1] = []
		
		adjacency[v0].append(v1)
		adjacency[v1].append(v0)
	
	# Sort each vertex's neighbors for consistency
	for v in adjacency.keys():
		adjacency[v].sort()
	
	# Create compact array
	var compact = PackedInt32Array()
	
	# Process vertices in order (you can change the order if needed)
	var vertices = adjacency.keys()
	vertices.sort()
	
	
	edge_stride.resize(vertices.size())
	var main: int = 0
	for v in vertices:
		var neighbors: int = 0
		compact.append(v)
			
		# Add all its neighbors
		for neighbor in adjacency[v]:
			compact.append(neighbor)
			neighbors += 1
		edge_stride[main] = neighbors
		main += 1
	
	edge_vertex_indices.resize(edge_stride.size())
	for i in range(edge_stride.size()):
		if i == 0:
			edge_vertex_indices[0] = 0
		for j in range(0, i):
			edge_vertex_indices[i] += edge_stride[j] + 1
	
	return compact

func construct_compact_tet(tets: Array, tet_stride: PackedInt32Array, tet_vertex_indices: PackedInt32Array) -> PackedInt32Array:
	# Build adjacency list
	var adjacency = {}
	
	# Add all tets to adjacency list
	for tet in tets:
		var t0 = tet[0]
		var t1 = tet[1]
		var t2 = tet[2]
		var t3 = tet[3]
		
		# Add bidirectional connections
		if not adjacency.has(t0):
			adjacency[t0] = []
		if not adjacency.has(t1):
			adjacency[t1] = []
		if not adjacency.has(t2):
			adjacency[t2] = []
		if not adjacency.has(t3):
			adjacency[t3] = []
		
		adjacency[t0].append(t1)
		adjacency[t0].append(t2)
		adjacency[t0].append(t3)
		
		adjacency[t1].append(t0)
		adjacency[t1].append(t2)
		adjacency[t1].append(t3)
		
		adjacency[t2].append(t0)
		adjacency[t2].append(t1)
		adjacency[t2].append(t3)
		
		adjacency[t3].append(t0)
		adjacency[t3].append(t1)
		adjacency[t3].append(t2)
	
	# Create compact array
	var compact = PackedInt32Array()
	
	var vertices = adjacency.keys()
	
	
	tet_stride.resize(vertices.size())
	var main: int = 0
	for v in vertices:
		var neighbors: int = 0
		compact.append(v)
			
		# Add all its neighbors
		for neighbor in adjacency[v]:
			compact.append(neighbor)
			neighbors += 1
		tet_stride[main] = neighbors
		main += 1
	
	tet_vertex_indices.resize(tet_stride.size())
	for i in range(tet_stride.size()):
		if i == 0:
			tet_vertex_indices[0] = 0
		for j in range(0, i):
			tet_vertex_indices[i] += tet_stride[j] + 1
	
	return compact
	
func construct_edge(edges: Array) -> PackedInt32Array:
	var edge_array: PackedInt32Array = PackedInt32Array()
	# Add all edges to adjacency list
	for edge in edges:
		var v0 = edge[0]
		var v1 = edge[1]
		edge_array.append(v0)
		edge_array.append(v1)
	
	return edge_array


func construct_tet(tets: Array) -> PackedInt32Array:
	# Add all tets to adjacency list
	var tet_array = PackedInt32Array()
	for tet in tets:
		var t0 = tet[0]
		var t1 = tet[1]
		var t2 = tet[2]
		var t3 = tet[3]
		
		tet_array.append(t0)
		tet_array.append(t1)
		tet_array.append(t2)
		tet_array.append(t3)
	return tet_array

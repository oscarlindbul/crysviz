// Shared indexed-mesh types for isosurface smoothing (included by
// marching_cubes.cpp and the iso_*.hpp smoothing modules).
//
// marching_cubes.cpp emits a triangle soup; when a mesh smoothing method is
// active it welds that soup into an IsoMesh (every MC vertex sits on exactly
// one grid edge, so the edge id is the weld key), hands it to one of the
// smoothers/subdividers, and de-indexes the result again for three.js.
//
// Positions are in the normalised grid box used by MarchingCubes
// (grid point i -> i/(nx-1)), so one voxel is (dx, dy, dz) there.
#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>

struct IsoMesh {
	std::vector<float> pos;      // 3 floats per vertex
	std::vector<float> val;      // 1 float per vertex when has_val (colour-by values), else empty
	std::vector<uint32_t> tri;   // 3 indices per triangle, consistent winding
	bool has_val = false;

	uint32_t vertex_count() const { return (uint32_t)(pos.size() / 3); }
	uint32_t triangle_count() const { return (uint32_t)(tri.size() / 3); }
};

inline uint64_t iso_edge_key(uint32_t a, uint32_t b) {
	if (a > b) std::swap(a, b);
	return ((uint64_t)a << 32) | (uint64_t)b;
}
inline uint32_t iso_edge_a(uint64_t key) { return (uint32_t)(key >> 32); }
inline uint32_t iso_edge_b(uint64_t key) { return (uint32_t)(key & 0xffffffffu); }

// Vertex adjacency and boundary information of an IsoMesh.
//
// A vertex is "boundary" when it touches an edge used by exactly one triangle
// (the open rims where the surface meets the cell faces or the NCI mask), or
// an edge used by more than two (non-manifold). Smoothers keep boundary
// vertices FIXED: the periodic images of the surface are the same mesh
// translated by a lattice vector, so moving a rim vertex would open a seam
// between neighbouring cells.
struct MeshTopology {
	uint32_t nv = 0;
	std::vector<uint32_t> ring_start;  // nv+1 offsets into ring (CSR)
	std::vector<uint32_t> ring;        // unique neighbour vertex ids, per vertex
	std::vector<uint8_t> boundary;     // 1 = boundary / non-manifold vertex (keep fixed)
	std::vector<uint64_t> edges;       // unique undirected edges, sorted (iso_edge_key)
	std::vector<uint8_t> edge_faces;   // triangles using each edge (saturates at 255)

	uint32_t valence(uint32_t v) const { return ring_start[v + 1] - ring_start[v]; }

	// Index of edge (a,b) in `edges`, or -1 when absent.
	int64_t find_edge(uint32_t a, uint32_t b) const {
		const uint64_t key = iso_edge_key(a, b);
		auto it = std::lower_bound(edges.begin(), edges.end(), key);
		if (it == edges.end() || *it != key) return -1;
		return (int64_t)(it - edges.begin());
	}
	bool is_boundary_edge(size_t e) const { return edge_faces[e] != 2; }
};

inline MeshTopology build_topology(const IsoMesh& m) {
	MeshTopology t;
	const uint32_t nv = m.vertex_count();
	const size_t nt = m.triangle_count();
	t.nv = nv;

	std::vector<uint64_t> all;
	all.reserve(nt * 3);
	for (size_t f = 0; f < nt; f++) {
		const uint32_t a = m.tri[3*f], b = m.tri[3*f+1], c = m.tri[3*f+2];
		all.push_back(iso_edge_key(a, b));
		all.push_back(iso_edge_key(b, c));
		all.push_back(iso_edge_key(c, a));
	}
	std::sort(all.begin(), all.end());

	t.edges.reserve(all.size() / 2 + 1);
	t.edge_faces.reserve(all.size() / 2 + 1);
	for (size_t i = 0; i < all.size();) {
		size_t j = i + 1;
		while (j < all.size() && all[j] == all[i]) j++;
		t.edges.push_back(all[i]);
		t.edge_faces.push_back((uint8_t)std::min<size_t>(j - i, 255));
		i = j;
	}

	t.boundary.assign(nv, 0);
	std::vector<uint32_t> deg(nv, 0);
	for (size_t e = 0; e < t.edges.size(); e++) {
		const uint32_t a = iso_edge_a(t.edges[e]), b = iso_edge_b(t.edges[e]);
		deg[a]++; deg[b]++;
		if (t.edge_faces[e] != 2) { t.boundary[a] = 1; t.boundary[b] = 1; }
	}
	t.ring_start.assign(nv + 1, 0);
	for (uint32_t v = 0; v < nv; v++) t.ring_start[v + 1] = t.ring_start[v] + deg[v];
	t.ring.resize(t.ring_start[nv]);
	std::vector<uint32_t> fill(t.ring_start.begin(), t.ring_start.end() - 1);
	for (size_t e = 0; e < t.edges.size(); e++) {
		const uint32_t a = iso_edge_a(t.edges[e]), b = iso_edge_b(t.edges[e]);
		t.ring[fill[a]++] = b;
		t.ring[fill[b]++] = a;
	}
	// Isolated vertices (valence 0) are treated as boundary so nothing moves them.
	for (uint32_t v = 0; v < nv; v++) if (deg[v] == 0) t.boundary[v] = 1;
	return t;
}

// Subdivision surfaces for the welded isosurface mesh (see iso_mesh.hpp).
//
//   iso_subdivide_loop          Loop (1987), every triangle -> 4 per level.
//   iso_subdivide_catmull_clark Catmull-Clark (1978) on the triangle mesh,
//                               every triangle -> 6 per level.
//
// Both keep boundary / non-manifold vertices (MeshTopology::boundary) EXACTLY
// fixed and put new points on boundary edges at the exact edge midpoint, so
// the open rims of neighbouring periodic images keep matching bit for bit
// (the rim polyline is only refined, never moved). The optional per-vertex
// scalar (m.val) is interpolated with the same weights as the positions.
//
// Size cap: a level is only applied while its output stays at or below
// ISO_SUBDIVIDE_MAX_TRIANGLES triangles (50M, i.e. 150M uint32 indices, well
// below 2^31); further levels are silently skipped.
#pragma once

#include "iso_mesh.hpp"

static const size_t ISO_SUBDIVIDE_MAX_TRIANGLES = 50000000;

// Edge ids (index into t.edges) of triangle f's edges (a,b), (b,c), (c,a).
inline std::vector<uint32_t> iso_triangle_edges(const IsoMesh& m, const MeshTopology& t) {
	const size_t nt = m.triangle_count();
	std::vector<uint32_t> te(nt * 3);
	for (size_t f = 0; f < nt; f++) {
		const uint32_t a = m.tri[3*f], b = m.tri[3*f+1], c = m.tri[3*f+2];
		te[3*f]   = (uint32_t)t.find_edge(a, b);
		te[3*f+1] = (uint32_t)t.find_edge(b, c);
		te[3*f+2] = (uint32_t)t.find_edge(c, a);
	}
	return te;
}

// Loop (1987) subdivision: every triangle -> 4, `levels` times.
inline void iso_subdivide_loop(IsoMesh& m, int levels) {
	const bool hv = m.has_val && m.val.size() == m.pos.size() / 3;
	for (int level = 0; level < levels; level++) {
		const uint32_t nv = m.vertex_count();
		const size_t nt = m.triangle_count();
		if (nt == 0 || nv == 0) return;
		if (nt * 4 > ISO_SUBDIVIDE_MAX_TRIANGLES) return;

		const MeshTopology t = build_topology(m);
		const size_t ne = t.edges.size();
		const std::vector<uint32_t> te = iso_triangle_edges(m, t);
		const size_t nv2 = (size_t)nv + ne;

		const std::vector<float>& P = m.pos;
		const std::vector<float>& V = m.val;
		std::vector<float> P2(nv2 * 3, 0.0f);
		std::vector<float> V2(hv ? nv2 : 0, 0.0f);

		// Old (even) vertices.
		for (uint32_t v = 0; v < nv; v++) {
			if (t.boundary[v]) {
				P2[3*v] = P[3*v]; P2[3*v+1] = P[3*v+1]; P2[3*v+2] = P[3*v+2];
				if (hv) V2[v] = V[v];
				continue;
			}
			const uint32_t n = t.valence(v);
			const double beta = n == 3 ? 3.0 / 16.0 : 3.0 / (8.0 * n);
			double sx = 0, sy = 0, sz = 0, sv = 0;
			for (uint32_t k = t.ring_start[v]; k < t.ring_start[v + 1]; k++) {
				const uint32_t u = t.ring[k];
				sx += P[3*u]; sy += P[3*u+1]; sz += P[3*u+2];
				if (hv) sv += V[u];
			}
			const double w = 1.0 - n * beta;
			P2[3*v]   = (float)(w * P[3*v]   + beta * sx);
			P2[3*v+1] = (float)(w * P[3*v+1] + beta * sy);
			P2[3*v+2] = (float)(w * P[3*v+2] + beta * sz);
			if (hv) V2[v] = (float)(w * V[v] + beta * sv);
		}

		// New (odd) edge vertices. An interior edge receives 3/16 (a+b) + 1/8 c
		// from each of its two triangles -> 3/8 (a+b) + 1/8 (c+d).
		for (size_t f = 0; f < nt; f++) {
			for (int k = 0; k < 3; k++) {
				const uint32_t a = m.tri[3*f + k], b = m.tri[3*f + (k+1)%3], c = m.tri[3*f + (k+2)%3];
				const uint32_t e = te[3*f + k];
				const size_t o = (size_t)nv + e;
				if (t.is_boundary_edge(e)) {
					for (int d = 0; d < 3; d++) P2[3*o+d] = 0.5f * (P[3*a+d] + P[3*b+d]);
					if (hv) V2[o] = 0.5f * (V[a] + V[b]);
				} else {
					for (int d = 0; d < 3; d++) P2[3*o+d] += 0.1875f * (P[3*a+d] + P[3*b+d]) + 0.125f * P[3*c+d];
					if (hv) V2[o] += 0.1875f * (V[a] + V[b]) + 0.125f * V[c];
				}
			}
		}

		// Children keep the parent's winding.
		std::vector<uint32_t> T2;
		T2.reserve(nt * 12);
		for (size_t f = 0; f < nt; f++) {
			const uint32_t a = m.tri[3*f], b = m.tri[3*f+1], c = m.tri[3*f+2];
			const uint32_t ab = nv + te[3*f], bc = nv + te[3*f+1], ca = nv + te[3*f+2];
			const uint32_t q[12] = { a, ab, ca,  ab, b, bc,  ca, bc, c,  ab, bc, ca };
			T2.insert(T2.end(), q, q + 12);
		}

		m.pos.swap(P2);
		m.tri.swap(T2);
		if (hv) m.val.swap(V2);
	}
}

// Catmull-Clark (1978) applied to the triangle mesh: each triangle -> 3 quads
// (face point, edge points, vertex), each quad split into 2 triangles along a
// consistent diagonal -> 6 triangles per input triangle per level.
//
// The quad (v, e_next, f, e_prev) is always split along the diagonal from the
// old vertex v to the face point f.
inline void iso_subdivide_catmull_clark(IsoMesh& m, int levels) {
	const bool hv = m.has_val && m.val.size() == m.pos.size() / 3;
	for (int level = 0; level < levels; level++) {
		const uint32_t nv = m.vertex_count();
		const size_t nt = m.triangle_count();
		if (nt == 0 || nv == 0) return;
		if (nt * 6 > ISO_SUBDIVIDE_MAX_TRIANGLES) return;

		const MeshTopology t = build_topology(m);
		const size_t ne = t.edges.size();
		const std::vector<uint32_t> te = iso_triangle_edges(m, t);
		const size_t fbase = (size_t)nv + ne;   // face points follow edge points
		const size_t nv2 = fbase + nt;

		const std::vector<float>& P = m.pos;
		const std::vector<float>& V = m.val;
		std::vector<float> P2(nv2 * 3, 0.0f);
		std::vector<float> V2(hv ? nv2 : 0, 0.0f);

		// Face points (centroids).
		for (size_t f = 0; f < nt; f++) {
			const uint32_t a = m.tri[3*f], b = m.tri[3*f+1], c = m.tri[3*f+2];
			const size_t o = fbase + f;
			for (int d = 0; d < 3; d++) P2[3*o+d] = (P[3*a+d] + P[3*b+d] + P[3*c+d]) * (1.0f / 3.0f);
			if (hv) V2[o] = (V[a] + V[b] + V[c]) * (1.0f / 3.0f);
		}

		// Edge points: interior (a + b + F1 + F2)/4, accumulated per triangle as
		// (a+b)/8 + F/4; boundary edges get the midpoint.
		for (size_t f = 0; f < nt; f++) {
			const size_t fo = fbase + f;
			for (int k = 0; k < 3; k++) {
				const uint32_t a = m.tri[3*f + k], b = m.tri[3*f + (k+1)%3];
				const uint32_t e = te[3*f + k];
				const size_t o = (size_t)nv + e;
				if (t.is_boundary_edge(e)) {
					for (int d = 0; d < 3; d++) P2[3*o+d] = 0.5f * (P[3*a+d] + P[3*b+d]);
					if (hv) V2[o] = 0.5f * (V[a] + V[b]);
				} else {
					for (int d = 0; d < 3; d++) P2[3*o+d] += 0.125f * (P[3*a+d] + P[3*b+d]) + 0.25f * P2[3*fo+d];
					if (hv) V2[o] += 0.125f * (V[a] + V[b]) + 0.25f * V2[fo];
				}
			}
		}

		// Vertex -> face incidence (CSR) for the averaged face points F.
		std::vector<uint32_t> vf_start(nv + 1, 0);
		for (size_t i = 0; i < nt * 3; i++) vf_start[m.tri[i] + 1]++;
		for (uint32_t v = 0; v < nv; v++) vf_start[v + 1] += vf_start[v];
		std::vector<uint32_t> vf(nt * 3);
		{
			std::vector<uint32_t> fill(vf_start.begin(), vf_start.end() - 1);
			for (size_t i = 0; i < nt * 3; i++) vf[fill[m.tri[i]]++] = (uint32_t)(i / 3);
		}

		// Old vertices: (F + 2R + (n-3) P) / n, with R = avg edge midpoint
		// = (P + avg ring)/2, so 2R = P + avg ring.
		for (uint32_t v = 0; v < nv; v++) {
			if (t.boundary[v] || vf_start[v] == vf_start[v + 1]) {
				P2[3*v] = P[3*v]; P2[3*v+1] = P[3*v+1]; P2[3*v+2] = P[3*v+2];
				if (hv) V2[v] = V[v];
				continue;
			}
			const uint32_t n = t.valence(v);
			double fx = 0, fy = 0, fz = 0, fv = 0;
			for (uint32_t k = vf_start[v]; k < vf_start[v + 1]; k++) {
				const size_t o = fbase + vf[k];
				fx += P2[3*o]; fy += P2[3*o+1]; fz += P2[3*o+2];
				if (hv) fv += V2[o];
			}
			const double fi = 1.0 / (vf_start[v + 1] - vf_start[v]);
			double rx = 0, ry = 0, rz = 0, rv = 0;
			for (uint32_t k = t.ring_start[v]; k < t.ring_start[v + 1]; k++) {
				const uint32_t u = t.ring[k];
				rx += P[3*u]; ry += P[3*u+1]; rz += P[3*u+2];
				if (hv) rv += V[u];
			}
			const double ri = 1.0 / n, ni = 1.0 / n;
			// (F + P + ring_avg + (n-3) P) / n = (F + ring_avg + (n-2) P) / n
			const double wp = (double)(n - 2) * ni;
			P2[3*v]   = (float)((fx * fi + rx * ri) * ni + wp * P[3*v]);
			P2[3*v+1] = (float)((fy * fi + ry * ri) * ni + wp * P[3*v+1]);
			P2[3*v+2] = (float)((fz * fi + rz * ri) * ni + wp * P[3*v+2]);
			if (hv) V2[v] = (float)((fv * fi + rv * ri) * ni + wp * V[v]);
		}

		// Children: quad (v, e_next, f, e_prev) -> (v, e_next, f) + (v, f, e_prev),
		// same winding as the parent.
		std::vector<uint32_t> T2;
		T2.reserve(nt * 18);
		for (size_t f = 0; f < nt; f++) {
			const uint32_t fp = (uint32_t)(fbase + f);
			for (int k = 0; k < 3; k++) {
				const uint32_t v = m.tri[3*f + k];
				const uint32_t en = nv + te[3*f + k];           // edge (v, next)
				const uint32_t ep = nv + te[3*f + (k+2)%3];     // edge (prev, v)
				const uint32_t q[6] = { v, en, fp,  v, fp, ep };
				T2.insert(T2.end(), q, q + 6);
			}
		}

		m.pos.swap(P2);
		m.tri.swap(T2);
		if (hv) m.val.swap(V2);
	}
}

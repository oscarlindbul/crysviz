// Vertex-relaxation smoothers for welded isosurface meshes (IsoMesh).
//
// All smoothers use uniform ("umbrella") weights over the one-ring from
// MeshTopology and Jacobi updates (every vertex reads the positions of the
// previous pass). Boundary vertices (t.boundary[v] != 0: periodic seams, mask
// rims, non-manifold and isolated vertices) are never moved, but still act as
// neighbours of the interior vertices. Only m.pos changes.
//
// Scratch buffers are allocated once per call and swapped with m.pos between
// passes, so there is no allocation inside the iteration loops.
#pragma once

#include "iso_mesh.hpp"

#include <cmath>
#include <utility>

// One uniform Laplacian pass: dst = src + lambda * (avg(neighbours) - src)
// for interior vertices, dst = src for boundary vertices.
inline void iso_laplacian_pass(const float* src, float* dst, const MeshTopology& t, float lambda) {
	const uint32_t nv = t.nv;
	const uint32_t* ring = t.ring.data();
	const uint32_t* rs = t.ring_start.data();
	const uint8_t* bnd = t.boundary.data();
	for (uint32_t v = 0; v < nv; v++) {
		const float x = src[3*v], y = src[3*v+1], z = src[3*v+2];
		const uint32_t r0 = rs[v], r1 = rs[v + 1];
		if (bnd[v] || r1 == r0) {
			dst[3*v] = x; dst[3*v+1] = y; dst[3*v+2] = z;
			continue;
		}
		float sx = 0.0f, sy = 0.0f, sz = 0.0f;
		for (uint32_t k = r0; k < r1; k++) {
			const float* q = src + 3 * (size_t)ring[k];
			sx += q[0]; sy += q[1]; sz += q[2];
		}
		const float s = lambda / (float)(r1 - r0);
		dst[3*v]   = x + (sx * s - lambda * x);
		dst[3*v+1] = y + (sy * s - lambda * y);
		dst[3*v+2] = z + (sz * s - lambda * z);
	}
}

// True when the mesh/topology pair can be smoothed at all.
inline bool iso_smooth_usable(const IsoMesh& m, const MeshTopology& t) {
	const uint32_t nv = m.vertex_count();
	return nv > 0 && t.nv == nv && t.ring_start.size() == (size_t)nv + 1 && t.boundary.size() == nv;
}

// plain umbrella Laplacian: x += lambda * (avg(neighbours) - x)
inline void iso_smooth_laplacian(IsoMesh& m, const MeshTopology& t, int iterations, float lambda) {
	if (iterations <= 0 || !iso_smooth_usable(m, t)) return;
	std::vector<float> tmp(m.pos.size());
	for (int it = 0; it < iterations; it++) {
		iso_laplacian_pass(m.pos.data(), tmp.data(), t, lambda);
		m.pos.swap(tmp);
	}
}

// Taubin lambda|mu: per iteration one pass with +lambda then one with mu (mu negative, |mu| > lambda)
inline void iso_smooth_taubin(IsoMesh& m, const MeshTopology& t, int iterations, float lambda, float mu) {
	if (iterations <= 0 || !iso_smooth_usable(m, t)) return;
	std::vector<float> tmp(m.pos.size());
	for (int it = 0; it < iterations; it++) {
		iso_laplacian_pass(m.pos.data(), tmp.data(), t, lambda);
		iso_laplacian_pass(tmp.data(), m.pos.data(), t, mu);
	}
}

// Vollmer, Mencl & Mueller 1999 "Improved Laplacian smoothing of noisy surface meshes" (HC algorithm), alpha in [0,1], beta in [0,1]
//
// o: original positions, q: positions of the previous iteration.
//   p = avg(neighbours of q)                       (interior; p = q on boundary)
//   b = p - (alpha * o + (1 - alpha) * q)          (b = 0 on boundary)
//   p -= beta * b + (1 - beta) / valence * sum_neighbours(b)
inline void iso_smooth_hc(IsoMesh& m, const MeshTopology& t, int iterations, float alpha, float beta) {
	if (iterations <= 0 || !iso_smooth_usable(m, t)) return;
	alpha = std::min(std::max(alpha, 0.0f), 1.0f);
	beta = std::min(std::max(beta, 0.0f), 1.0f);
	const uint32_t nv = t.nv;
	const uint32_t* ring = t.ring.data();
	const uint32_t* rs = t.ring_start.data();
	const uint8_t* bnd = t.boundary.data();

	const std::vector<float> o = m.pos;
	std::vector<float> p(m.pos.size());
	std::vector<float> b(m.pos.size());
	const float ia = 1.0f - alpha, ib = 1.0f - beta;

	for (int it = 0; it < iterations; it++) {
		const float* q = m.pos.data();
		iso_laplacian_pass(q, p.data(), t, 1.0f);
		for (uint32_t v = 0; v < nv; v++) {
			const size_t i = 3 * (size_t)v;
			if (bnd[v] || rs[v + 1] == rs[v]) { b[i] = b[i+1] = b[i+2] = 0.0f; continue; }
			b[i]   = p[i]   - (alpha * o[i]   + ia * q[i]);
			b[i+1] = p[i+1] - (alpha * o[i+1] + ia * q[i+1]);
			b[i+2] = p[i+2] - (alpha * o[i+2] + ia * q[i+2]);
		}
		for (uint32_t v = 0; v < nv; v++) {
			const uint32_t r0 = rs[v], r1 = rs[v + 1];
			if (bnd[v] || r1 == r0) continue;   // p == q already
			float sx = 0.0f, sy = 0.0f, sz = 0.0f;
			for (uint32_t k = r0; k < r1; k++) {
				const float* bj = b.data() + 3 * (size_t)ring[k];
				sx += bj[0]; sy += bj[1]; sz += bj[2];
			}
			const float s = ib / (float)(r1 - r0);
			const size_t i = 3 * (size_t)v;
			p[i]   -= beta * b[i]   + s * sx;
			p[i+1] -= beta * b[i+1] + s * sy;
			p[i+2] -= beta * b[i+2] + s * sz;
		}
		m.pos.swap(p);
	}
}

// Helper: Taubin mu from lambda and pass-band frequency k_PB (Taubin 1995: 1/lambda + 1/mu = k_PB). Clamp to sensible negative value.
//
// Valid when 0 < lambda < 1 and 0 < k_PB < 1/lambda; Taubin suggests
// k_PB in [0.01, 0.1]. The result is kept in [-1, -1.001*lambda] so that
// |mu| > lambda (no net shrinking) and the mu pass stays stable.
inline float iso_taubin_mu(float lambda, float k_pb) {
	if (!(lambda > 0.0f)) return 0.0f;
	lambda = std::min(lambda, 0.99f);
	const float lo = -1.0f, hi = -1.001f * lambda;
	const float d = k_pb - 1.0f / lambda;   // 1/mu
	if (!(d < 0.0f) || !std::isfinite(d)) return hi;
	const float mu = 1.0f / d;
	return std::min(std::max(mu, lo), hi);
}

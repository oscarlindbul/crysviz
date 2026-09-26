#pragma once

// Smoothing filters for scalar grids fed to marching cubes (marching_cubes.cpp).
// Grid layout: float32, x fastest, idx = i + nx*(j + ny*k). Grid point i maps to i/(n-1) in the box.
// Periodic grids hold n points per axis, point n is point 0 (not duplicated).
// Optional mask: points with v >= mask_value are "masked" and marching cubes skips any cube touching one.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstddef>
#include <algorithm>
#include <vector>

namespace iso_filters_detail {

inline int64_t wrap_index(int64_t i, int64_t n) {
	i %= n;
	return i < 0 ? i + n : i;
}

inline int64_t clamp_index(int64_t i, int64_t n) {
	return i < 0 ? 0 : (i >= n ? n - 1 : i);
}

inline bool is_masked(float v, bool has_mask, float mask_value) {
	return !std::isfinite(v) || (has_mask && v >= mask_value);
}

/////////////////////////////////////
// Gaussian
/////////////////////////////////////

// One 1-D pass along `axis` (0 = x, 1 = y, 2 = z). m is the original mask (nullptr = none):
// masked points are copied from src, and excluded from their neighbours' sums.
inline void gauss_pass(const float* src, float* dst, const uint8_t* m,
		uint32_t nx, uint32_t ny, uint32_t nz, int axis,
		const std::vector<float>& w, int r, bool periodic) {
	const int64_t taps = 2 * r + 1;
	if (axis == 0) {
		// Pad each row with wrapped / clamped values, then a plain convolution
		const int64_t n = nx;
		std::vector<float> row(n + 2 * r);
		std::vector<uint8_t> mrow(m ? n + 2 * r : 0);
		for (int64_t line = 0; line < (int64_t)ny * nz; line++) {
			const float* s = src + line * n;
			float* d = dst + line * n;
			const uint8_t* ms = m ? m + line * n : nullptr;
			for (int64_t p = 0; p < n + 2 * r; p++) {
				int64_t q = periodic ? wrap_index(p - r, n) : clamp_index(p - r, n);
				row[p] = s[q];
				if (m) mrow[p] = ms[q];
			}
			if (!m) {
				for (int64_t i = 0; i < n; i++) {
					float acc = 0.0f;
					const float* rp = row.data() + i;
					for (int64_t t = 0; t < taps; t++) acc += w[t] * rp[t];
					d[i] = acc;
				}
			} else {
				for (int64_t i = 0; i < n; i++) {
					if (ms[i]) { d[i] = s[i]; continue; }
					float acc = 0.0f, wsum = 0.0f;
					const float* rp = row.data() + i;
					const uint8_t* mp = mrow.data() + i;
					for (int64_t t = 0; t < taps; t++) {
						if (mp[t]) continue;
						acc += w[t] * rp[t];
						wsum += w[t];
					}
					d[i] = acc / wsum; // wsum > 0: the centre sample is unmasked
				}
			}
		}
		return;
	}

	// y / z: accumulate whole x rows (contiguous, cache friendly)
	const int64_t n = axis == 1 ? ny : nz;
	const int64_t outer = axis == 1 ? nz : ny;
	const int64_t step = axis == 1 ? (int64_t)nx : (int64_t)nx * ny;
	const int64_t outer_step = axis == 1 ? (int64_t)nx * ny : (int64_t)nx;
	std::vector<float> acc(nx), wsum(m ? nx : 0);
	for (int64_t o = 0; o < outer; o++) {
		for (int64_t p = 0; p < n; p++) {
			std::fill(acc.begin(), acc.end(), 0.0f);
			if (m) std::fill(wsum.begin(), wsum.end(), 0.0f);
			for (int64_t t = 0; t < taps; t++) {
				int64_t q = periodic ? wrap_index(p + t - r, n) : clamp_index(p + t - r, n);
				const int64_t off = o * outer_step + q * step;
				const float* s = src + off;
				const float wt = w[t];
				if (!m) {
					for (uint32_t i = 0; i < nx; i++) acc[i] += wt * s[i];
				} else {
					const uint8_t* ms = m + off;
					for (uint32_t i = 0; i < nx; i++) {
						if (ms[i]) continue;
						acc[i] += wt * s[i];
						wsum[i] += wt;
					}
				}
			}
			const int64_t off = o * outer_step + p * step;
			float* d = dst + off;
			if (!m) {
				std::memcpy(d, acc.data(), sizeof(float) * nx);
			} else {
				const float* s = src + off;
				const uint8_t* ms = m + off;
				for (uint32_t i = 0; i < nx; i++) d[i] = ms[i] ? s[i] : acc[i] / wsum[i];
			}
		}
	}
}

/////////////////////////////////////
// Refinement
/////////////////////////////////////

// Per-axis precomputed stencil for every output coordinate
struct AxisTaps {
	uint32_t n_out = 0;
	std::vector<int32_t> idx; // n_out*4 source indices (may repeat, unused taps have weight 0)
	std::vector<float> w;     // n_out*4 cubic weights
	std::vector<int32_t> lin0, lin1; // trilinear neighbours
	std::vector<float> lin_t;        // trilinear weight of lin1 (0 => only lin0 contributes)
};

inline void cubic_weights(float t, int kernel, float* ww) {
	const float t2 = t * t, t3 = t2 * t;
	if (kernel == 1) {
		// Uniform cubic B-spline
		const float u = 1.0f - t;
		ww[0] = u * u * u / 6.0f;
		ww[1] = (3.0f * t3 - 6.0f * t2 + 4.0f) / 6.0f;
		ww[2] = (-3.0f * t3 + 3.0f * t2 + 3.0f * t + 1.0f) / 6.0f;
		ww[3] = t3 / 6.0f;
	} else {
		// Catmull-Rom (t = 0 gives exactly 0,1,0,0)
		ww[0] = 0.5f * (-t3 + 2.0f * t2 - t);
		ww[1] = 0.5f * (3.0f * t3 - 5.0f * t2 + 2.0f);
		ww[2] = 0.5f * (-3.0f * t3 + 4.0f * t2 + t);
		ww[3] = 0.5f * (t3 - t2);
	}
}

inline AxisTaps build_axis_taps(uint32_t n, int factor, int kernel, bool periodic) {
	AxisTaps a;
	const uint32_t N = (n - 1) * (uint32_t)factor + 1;
	a.n_out = N;
	a.idx.assign((size_t)N * 4, 0);
	a.w.assign((size_t)N * 4, 0.0f);
	a.lin0.resize(N);
	a.lin1.resize(N);
	a.lin_t.resize(N);
	for (uint32_t o = 0; o < N; o++) {
		const int64_t i0 = o / (uint32_t)factor;
		const float t = (float)(o % (uint32_t)factor) / (float)factor;
		float ww[4];
		cubic_weights(t, kernel, ww);
		// Fold the 4 taps onto distinct in-range indices
		int32_t ci[4];
		float cw[4];
		int cnt = 0;
		auto add = [&](int64_t j, float wt) {
			for (int c = 0; c < cnt; c++) if (ci[c] == j) { cw[c] += wt; return; }
			ci[cnt] = (int32_t)j; cw[cnt] = wt; cnt++;
		};
		for (int s = 0; s < 4; s++) {
			const int64_t j = i0 - 1 + s;
			if (ww[s] == 0.0f) continue;
			if (periodic) add(wrap_index(j, n), ww[s]);
			else if (j >= 0 && j < n) add(j, ww[s]);
			else if (n < 2) add(0, ww[s]);
			else if (j < 0) {
				// Linear extrapolation of the ghost sample: f(-d) = (1+d) f(0) - d f(1)
				const float d = (float)(-j);
				add(0, (1.0f + d) * ww[s]);
				add(1, -d * ww[s]);
			} else {
				const float d = (float)(j - (n - 1));
				add(n - 1, (1.0f + d) * ww[s]);
				add(n - 2, -d * ww[s]);
			}
		}
		for (int c = 0; c < 4; c++) {
			a.idx[(size_t)o * 4 + c] = c < cnt ? ci[c] : (cnt ? ci[0] : 0);
			a.w[(size_t)o * 4 + c] = c < cnt ? cw[c] : 0.0f;
		}
		a.lin0[o] = (int32_t)(periodic ? wrap_index(i0, n) : clamp_index(i0, n));
		a.lin1[o] = (int32_t)(periodic ? wrap_index(i0 + 1, n) : clamp_index(i0 + 1, n));
		a.lin_t[o] = t;
	}
	return a;
}

} // namespace iso_filters_detail

// Separable Gaussian blur, sigma in voxels, radius ceil(3 sigma). Axis ends wrap when periodic, else clamp.
// Masked points (has_mask && v >= mask_value, or non-finite) are copied unchanged and are excluded from
// their neighbours' weights (weights renormalised over unmasked samples), so NCI mask walls don't smear in.
// src and dst may NOT alias. Uses an internal scratch buffer.
inline void iso_gaussian_filter(const float* src, float* dst, uint32_t nx, uint32_t ny, uint32_t nz,
		float sigma, bool periodic, bool has_mask, float mask_value) {
	using namespace iso_filters_detail;
	const size_t total = (size_t)nx * ny * nz;
	if (total == 0) return;
	if (!(sigma > 0.0f)) {
		std::memcpy(dst, src, sizeof(float) * total);
		return;
	}
	const int r = std::max(1, (int)std::ceil(3.0f * sigma));
	std::vector<float> w(2 * r + 1);
	float wsum = 0.0f;
	for (int t = -r; t <= r; t++) {
		w[t + r] = std::exp(-0.5f * (float)(t * t) / (sigma * sigma));
		wsum += w[t + r];
	}
	for (float& v : w) v /= wsum;

	// Mask (also catches non-finite values, which would otherwise poison the neighbourhood)
	std::vector<uint8_t> mask;
	bool any_masked = false;
	for (size_t p = 0; p < total && !any_masked; p++) any_masked = is_masked(src[p], has_mask, mask_value);
	if (any_masked) {
		mask.resize(total);
		for (size_t p = 0; p < total; p++) mask[p] = is_masked(src[p], has_mask, mask_value);
	}
	const uint8_t* m = any_masked ? mask.data() : nullptr;

	// x: src -> dst, y: dst -> scratch, z: scratch -> dst. Masked points keep their original value in
	// every intermediate, so copying them through each pass leaves them unchanged.
	std::vector<float> scratch(total);
	gauss_pass(src, dst, m, nx, ny, nz, 0, w, r, periodic);
	gauss_pass(dst, scratch.data(), m, nx, ny, nz, 1, w, r, periodic);
	gauss_pass(scratch.data(), dst, m, nx, ny, nz, 2, w, r, periodic);
}

// Output size per axis of iso_refine_field: (n-1)*factor + 1 (keeps grid end points => same [0,1] box).
inline uint32_t iso_refined_dim(uint32_t n, int factor) {
	if (factor <= 1 || n < 2) return n;
	return (n - 1) * (uint32_t)factor + 1;
}

// Upsample by an integer factor with a separable 4-tap cubic kernel:
// kernel 0 = Catmull-Rom (interpolating: reproduces original samples exactly at coincident points),
// kernel 1 = uniform cubic B-spline (approximating, smoothing).
// Stencil indices outside [0,n-1] wrap when periodic; otherwise the ghost samples are linearly
// extrapolated from the two end samples (so linear ramps are reproduced up to the grid ends).
// Only taps with a non-zero weight count as support (at coincident points Catmull-Rom uses 1 sample).
// Mask handling: if any of the 4x4x4 support samples is masked, fall back to trilinear of the 8
// surrounding samples; if any of those 8 is masked, output the maximum of them (>= mask_value, so the
// cube stays masked). dst has iso_refined_dim(nx)*iso_refined_dim(ny)*iso_refined_dim(nz) floats.
inline void iso_refine_field(const float* src, uint32_t nx, uint32_t ny, uint32_t nz, int factor, int kernel,
		bool periodic, bool has_mask, float mask_value, float* dst) {
	using namespace iso_filters_detail;
	const size_t total = (size_t)nx * ny * nz;
	if (total == 0) return;
	if (factor <= 1) {
		std::memcpy(dst, src, sizeof(float) * total);
		return;
	}
	const AxisTaps ax = build_axis_taps(nx, factor, kernel, periodic);
	const AxisTaps ay = build_axis_taps(ny, factor, kernel, periodic);
	const AxisTaps az = build_axis_taps(nz, factor, kernel, periodic);
	const size_t NX = ax.n_out, NY = ay.n_out, NZ = az.n_out;

	// x pass: src (nx,ny,nz) -> A (NX,ny,nz)
	std::vector<float> A(NX * ny * nz);
	for (size_t line = 0; line < (size_t)ny * nz; line++) {
		const float* s = src + line * nx;
		float* d = A.data() + line * NX;
		for (size_t I = 0; I < NX; I++) {
			const int32_t* ii = &ax.idx[I * 4];
			const float* ww = &ax.w[I * 4];
			d[I] = ww[0] * s[ii[0]] + ww[1] * s[ii[1]] + ww[2] * s[ii[2]] + ww[3] * s[ii[3]];
		}
	}

	// y pass: A -> B (NX,NY,nz), whole rows at a time
	std::vector<float> B(NX * NY * nz);
	for (size_t k = 0; k < nz; k++) {
		for (size_t J = 0; J < NY; J++) {
			float* d = B.data() + NX * (J + NY * k);
			std::fill(d, d + NX, 0.0f);
			for (int c = 0; c < 4; c++) {
				const float wt = ay.w[J * 4 + c];
				if (wt == 0.0f) continue;
				const float* s = A.data() + NX * (ay.idx[J * 4 + c] + (size_t)ny * k);
				for (size_t I = 0; I < NX; I++) d[I] += wt * s[I];
			}
		}
	}
	A.clear();
	A.shrink_to_fit();

	// z pass: B -> dst (NX,NY,NZ), whole planes at a time
	const size_t plane = NX * NY;
	for (size_t K = 0; K < NZ; K++) {
		float* d = dst + plane * K;
		std::fill(d, d + plane, 0.0f);
		for (int c = 0; c < 4; c++) {
			const float wt = az.w[K * 4 + c];
			if (wt == 0.0f) continue;
			const float* s = B.data() + plane * az.idx[K * 4 + c];
			for (size_t p = 0; p < plane; p++) d[p] += wt * s[p];
		}
	}
	B.clear();
	B.shrink_to_fit();

	if (!has_mask) return;

	// Mask patching. Dilate the source mask separably over the non-zero cubic support,
	// then redo the flagged output points with the trilinear / max fallback.
	std::vector<uint8_t> M0(total);
	bool any_masked = false;
	for (size_t p = 0; p < total; p++) {
		M0[p] = is_masked(src[p], true, mask_value);
		any_masked |= M0[p] != 0;
	}
	if (!any_masked) return;

	std::vector<uint8_t> Mx(NX * ny * nz);
	for (size_t line = 0; line < (size_t)ny * nz; line++) {
		const uint8_t* s = M0.data() + line * nx;
		uint8_t* d = Mx.data() + line * NX;
		for (size_t I = 0; I < NX; I++) {
			uint8_t f = 0;
			for (int c = 0; c < 4; c++) if (ax.w[I * 4 + c] != 0.0f) f |= s[ax.idx[I * 4 + c]];
			d[I] = f;
		}
	}
	std::vector<uint8_t> My(NX * NY * nz);
	for (size_t k = 0; k < nz; k++) {
		for (size_t J = 0; J < NY; J++) {
			uint8_t* d = My.data() + NX * (J + NY * k);
			std::fill(d, d + NX, (uint8_t)0);
			for (int c = 0; c < 4; c++) {
				if (ay.w[J * 4 + c] == 0.0f) continue;
				const uint8_t* s = Mx.data() + NX * (ay.idx[J * 4 + c] + (size_t)ny * k);
				for (size_t I = 0; I < NX; I++) d[I] |= s[I];
			}
		}
	}
	Mx.clear();
	Mx.shrink_to_fit();

	std::vector<uint8_t> Mz(plane);
	for (size_t K = 0; K < NZ; K++) {
		std::fill(Mz.begin(), Mz.end(), (uint8_t)0);
		for (int c = 0; c < 4; c++) {
			if (az.w[K * 4 + c] == 0.0f) continue;
			const uint8_t* s = My.data() + plane * az.idx[K * 4 + c];
			for (size_t p = 0; p < plane; p++) Mz[p] |= s[p];
		}
		// Trilinear corners along z (a zero weight drops the second neighbour)
		const int zc = az.lin_t[K] > 0.0f ? 2 : 1;
		const int32_t zi[2] = { az.lin0[K], az.lin1[K] };
		const float zw[2] = { 1.0f - az.lin_t[K], az.lin_t[K] };
		for (size_t J = 0; J < NY; J++) {
			const int yc = ay.lin_t[J] > 0.0f ? 2 : 1;
			const int32_t yi[2] = { ay.lin0[J], ay.lin1[J] };
			const float yw[2] = { 1.0f - ay.lin_t[J], ay.lin_t[J] };
			for (size_t I = 0; I < NX; I++) {
				if (!Mz[I + NX * J]) continue;
				const int xc = ax.lin_t[I] > 0.0f ? 2 : 1;
				const int32_t xi[2] = { ax.lin0[I], ax.lin1[I] };
				const float xw[2] = { 1.0f - ax.lin_t[I], ax.lin_t[I] };
				float acc = 0.0f, vmax = -INFINITY;
				bool masked = false;
				for (int c = 0; c < zc; c++) {
					for (int b = 0; b < yc; b++) {
						const size_t row = (size_t)nx * (yi[b] + (size_t)ny * zi[c]);
						for (int a = 0; a < xc; a++) {
							const size_t p = row + xi[a];
							const float v = src[p];
							if (M0[p]) {
								masked = true;
								if (!(v <= vmax)) vmax = v; // NaN propagates
							}
							acc += xw[a] * yw[b] * zw[c] * v;
						}
					}
				}
				dst[I + NX * J + plane * K] = masked ? vmax : acc;
			}
		}
	}
}

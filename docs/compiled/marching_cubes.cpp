#include <math.h>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <emscripten/bind.h>

#include "iso_mesh.hpp"
#include "iso_smooth.hpp"
#include "iso_subdivide.hpp"
#include "iso_field_filters.hpp"

/////////////////////////////////////
// Marching cubes lookup tables
/////////////////////////////////////

// These tables are straight from Paul Bourke's page:
// http://paulbourke.net/geometry/polygonise/
// who in turn got them from Cory Gene Bloyd.

const unsigned short edgeTable[256] = {
	0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33, 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99, 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};

const unsigned char triTable[256][16] = {
	{ 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 8, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 1, 9, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 8, 3, 9, 8, 1, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 2, 10, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 8, 3, 1, 2, 10, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 2, 10, 0, 2, 9, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 2, 8, 3, 2, 10, 8, 10, 9, 8, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 11, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 11, 2, 8, 11, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 9, 0, 2, 3, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 11, 2, 1, 9, 11, 9, 8, 11, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 10, 1, 11, 10, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 10, 1, 0, 8, 10, 8, 11, 10, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 9, 0, 3, 11, 9, 11, 10, 9, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 8, 10, 10, 8, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 7, 8, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 3, 0, 7, 3, 4, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 1, 9, 8, 4, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 1, 9, 4, 7, 1, 7, 3, 1, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 2, 10, 8, 4, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 4, 7, 3, 0, 4, 1, 2, 10, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 2, 10, 9, 0, 2, 8, 4, 7, 255, 255, 255, 255, 255, 255, 255 },
	{ 2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, 255, 255, 255, 255 },
	{ 8, 4, 7, 3, 11, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 11, 4, 7, 11, 2, 4, 2, 0, 4, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 0, 1, 8, 4, 7, 2, 3, 11, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, 255, 255, 255, 255 },
	{ 3, 10, 1, 3, 11, 10, 7, 8, 4, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, 255, 255, 255, 255 },
	{ 4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, 255, 255, 255, 255 },
	{ 4, 7, 11, 4, 11, 9, 9, 11, 10, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 5, 4, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 5, 4, 0, 8, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 5, 4, 1, 5, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 8, 5, 4, 8, 3, 5, 3, 1, 5, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 2, 10, 9, 5, 4, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 0, 8, 1, 2, 10, 4, 9, 5, 255, 255, 255, 255, 255, 255, 255 },
	{ 5, 2, 10, 5, 4, 2, 4, 0, 2, 255, 255, 255, 255, 255, 255, 255 },
	{ 2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, 255, 255, 255, 255 },
	{ 9, 5, 4, 2, 3, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 11, 2, 0, 8, 11, 4, 9, 5, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 5, 4, 0, 1, 5, 2, 3, 11, 255, 255, 255, 255, 255, 255, 255 },
	{ 2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, 255, 255, 255, 255 },
	{ 10, 3, 11, 10, 1, 3, 9, 5, 4, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, 255, 255, 255, 255 },
	{ 5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, 255, 255, 255, 255 },
	{ 5, 4, 8, 5, 8, 10, 10, 8, 11, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 7, 8, 5, 7, 9, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 3, 0, 9, 5, 3, 5, 7, 3, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 7, 8, 0, 1, 7, 1, 5, 7, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 5, 3, 3, 5, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 7, 8, 9, 5, 7, 10, 1, 2, 255, 255, 255, 255, 255, 255, 255 },
	{ 10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, 255, 255, 255, 255 },
	{ 8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, 255, 255, 255, 255 },
	{ 2, 10, 5, 2, 5, 3, 3, 5, 7, 255, 255, 255, 255, 255, 255, 255 },
	{ 7, 9, 5, 7, 8, 9, 3, 11, 2, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, 255, 255, 255, 255 },
	{ 2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, 255, 255, 255, 255 },
	{ 11, 2, 1, 11, 1, 7, 7, 1, 5, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, 255, 255, 255, 255 },
	{ 5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, 255 },
	{ 11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, 255 },
	{ 11, 10, 5, 7, 11, 5, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 10, 6, 5, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 8, 3, 5, 10, 6, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 0, 1, 5, 10, 6, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 8, 3, 1, 9, 8, 5, 10, 6, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 6, 5, 2, 6, 1, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 6, 5, 1, 2, 6, 3, 0, 8, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 6, 5, 9, 0, 6, 0, 2, 6, 255, 255, 255, 255, 255, 255, 255 },
	{ 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, 255, 255, 255, 255 },
	{ 2, 3, 11, 10, 6, 5, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 11, 0, 8, 11, 2, 0, 10, 6, 5, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 1, 9, 2, 3, 11, 5, 10, 6, 255, 255, 255, 255, 255, 255, 255 },
	{ 5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, 255, 255, 255, 255 },
	{ 6, 3, 11, 6, 5, 3, 5, 1, 3, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, 255, 255, 255, 255 },
	{ 3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, 255, 255, 255, 255 },
	{ 6, 5, 9, 6, 9, 11, 11, 9, 8, 255, 255, 255, 255, 255, 255, 255 },
	{ 5, 10, 6, 4, 7, 8, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 3, 0, 4, 7, 3, 6, 5, 10, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 9, 0, 5, 10, 6, 8, 4, 7, 255, 255, 255, 255, 255, 255, 255 },
	{ 10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, 255, 255, 255, 255 },
	{ 6, 1, 2, 6, 5, 1, 4, 7, 8, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, 255, 255, 255, 255 },
	{ 8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, 255, 255, 255, 255 },
	{ 7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, 255 },
	{ 3, 11, 2, 7, 8, 4, 10, 6, 5, 255, 255, 255, 255, 255, 255, 255 },
	{ 5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, 255, 255, 255, 255 },
	{ 0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, 255, 255, 255, 255 },
	{ 9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, 255 },
	{ 8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, 255, 255, 255, 255 },
	{ 5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, 255 },
	{ 0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, 255 },
	{ 6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, 255, 255, 255, 255 },
	{ 10, 4, 9, 6, 4, 10, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 10, 6, 4, 9, 10, 0, 8, 3, 255, 255, 255, 255, 255, 255, 255 },
	{ 10, 0, 1, 10, 6, 0, 6, 4, 0, 255, 255, 255, 255, 255, 255, 255 },
	{ 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, 255, 255, 255, 255 },
	{ 1, 4, 9, 1, 2, 4, 2, 6, 4, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, 255, 255, 255, 255 },
	{ 0, 2, 4, 4, 2, 6, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 8, 3, 2, 8, 2, 4, 4, 2, 6, 255, 255, 255, 255, 255, 255, 255 },
	{ 10, 4, 9, 10, 6, 4, 11, 2, 3, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, 255, 255, 255, 255 },
	{ 3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, 255, 255, 255, 255 },
	{ 6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, 255 },
	{ 9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, 255, 255, 255, 255 },
	{ 8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, 255 },
	{ 3, 11, 6, 3, 6, 0, 0, 6, 4, 255, 255, 255, 255, 255, 255, 255 },
	{ 6, 4, 8, 11, 6, 8, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 7, 10, 6, 7, 8, 10, 8, 9, 10, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, 255, 255, 255, 255 },
	{ 10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, 255, 255, 255, 255 },
	{ 10, 6, 7, 10, 7, 1, 1, 7, 3, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, 255, 255, 255, 255 },
	{ 2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, 255 },
	{ 7, 8, 0, 7, 0, 6, 6, 0, 2, 255, 255, 255, 255, 255, 255, 255 },
	{ 7, 3, 2, 6, 7, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, 255, 255, 255, 255 },
	{ 2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, 255 },
	{ 1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, 255 },
	{ 11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, 255, 255, 255, 255 },
	{ 8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, 255 },
	{ 0, 9, 1, 11, 6, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, 255, 255, 255, 255 },
	{ 7, 11, 6, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 7, 6, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 0, 8, 11, 7, 6, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 1, 9, 11, 7, 6, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 8, 1, 9, 8, 3, 1, 11, 7, 6, 255, 255, 255, 255, 255, 255, 255 },
	{ 10, 1, 2, 6, 11, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 2, 10, 3, 0, 8, 6, 11, 7, 255, 255, 255, 255, 255, 255, 255 },
	{ 2, 9, 0, 2, 10, 9, 6, 11, 7, 255, 255, 255, 255, 255, 255, 255 },
	{ 6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, 255, 255, 255, 255 },
	{ 7, 2, 3, 6, 2, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 7, 0, 8, 7, 6, 0, 6, 2, 0, 255, 255, 255, 255, 255, 255, 255 },
	{ 2, 7, 6, 2, 3, 7, 0, 1, 9, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, 255, 255, 255, 255 },
	{ 10, 7, 6, 10, 1, 7, 1, 3, 7, 255, 255, 255, 255, 255, 255, 255 },
	{ 10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, 255, 255, 255, 255 },
	{ 0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, 255, 255, 255, 255 },
	{ 7, 6, 10, 7, 10, 8, 8, 10, 9, 255, 255, 255, 255, 255, 255, 255 },
	{ 6, 8, 4, 11, 8, 6, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 6, 11, 3, 0, 6, 0, 4, 6, 255, 255, 255, 255, 255, 255, 255 },
	{ 8, 6, 11, 8, 4, 6, 9, 0, 1, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, 255, 255, 255, 255 },
	{ 6, 8, 4, 6, 11, 8, 2, 10, 1, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, 255, 255, 255, 255 },
	{ 4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, 255, 255, 255, 255 },
	{ 10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, 255 },
	{ 8, 2, 3, 8, 4, 2, 4, 6, 2, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 4, 2, 4, 6, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, 255, 255, 255, 255 },
	{ 1, 9, 4, 1, 4, 2, 2, 4, 6, 255, 255, 255, 255, 255, 255, 255 },
	{ 8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, 255, 255, 255, 255 },
	{ 10, 1, 0, 10, 0, 6, 6, 0, 4, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, 255 },
	{ 10, 9, 4, 6, 10, 4, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 9, 5, 7, 6, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 8, 3, 4, 9, 5, 11, 7, 6, 255, 255, 255, 255, 255, 255, 255 },
	{ 5, 0, 1, 5, 4, 0, 7, 6, 11, 255, 255, 255, 255, 255, 255, 255 },
	{ 11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, 255, 255, 255, 255 },
	{ 9, 5, 4, 10, 1, 2, 7, 6, 11, 255, 255, 255, 255, 255, 255, 255 },
	{ 6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, 255, 255, 255, 255 },
	{ 7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, 255, 255, 255, 255 },
	{ 3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, 255 },
	{ 7, 2, 3, 7, 6, 2, 5, 4, 9, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, 255, 255, 255, 255 },
	{ 3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, 255, 255, 255, 255 },
	{ 6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, 255 },
	{ 9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, 255, 255, 255, 255 },
	{ 1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, 255 },
	{ 4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, 255 },
	{ 7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, 255, 255, 255, 255 },
	{ 6, 9, 5, 6, 11, 9, 11, 8, 9, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, 255, 255, 255, 255 },
	{ 0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, 255, 255, 255, 255 },
	{ 6, 11, 3, 6, 3, 5, 5, 3, 1, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, 255, 255, 255, 255 },
	{ 0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, 255 },
	{ 11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, 255 },
	{ 6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, 255, 255, 255, 255 },
	{ 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, 255, 255, 255, 255 },
	{ 9, 5, 6, 9, 6, 0, 0, 6, 2, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, 255 },
	{ 1, 5, 6, 2, 1, 6, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, 255 },
	{ 10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, 255, 255, 255, 255 },
	{ 0, 3, 8, 5, 6, 10, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 10, 5, 6, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 11, 5, 10, 7, 5, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 11, 5, 10, 11, 7, 5, 8, 3, 0, 255, 255, 255, 255, 255, 255, 255 },
	{ 5, 11, 7, 5, 10, 11, 1, 9, 0, 255, 255, 255, 255, 255, 255, 255 },
	{ 10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, 255, 255, 255, 255 },
	{ 11, 1, 2, 11, 7, 1, 7, 5, 1, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, 255, 255, 255, 255 },
	{ 9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, 255, 255, 255, 255 },
	{ 7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, 255 },
	{ 2, 5, 10, 2, 3, 5, 3, 7, 5, 255, 255, 255, 255, 255, 255, 255 },
	{ 8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, 255, 255, 255, 255 },
	{ 9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, 255, 255, 255, 255 },
	{ 9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, 255 },
	{ 1, 3, 5, 3, 7, 5, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 8, 7, 0, 7, 1, 1, 7, 5, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 0, 3, 9, 3, 5, 5, 3, 7, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 8, 7, 5, 9, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 5, 8, 4, 5, 10, 8, 10, 11, 8, 255, 255, 255, 255, 255, 255, 255 },
	{ 5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, 255, 255, 255, 255 },
	{ 0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, 255, 255, 255, 255 },
	{ 10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, 255 },
	{ 2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, 255, 255, 255, 255 },
	{ 0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, 255 },
	{ 0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, 255 },
	{ 9, 4, 5, 2, 11, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, 255, 255, 255, 255 },
	{ 5, 10, 2, 5, 2, 4, 4, 2, 0, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, 255 },
	{ 5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, 255, 255, 255, 255 },
	{ 8, 4, 5, 8, 5, 3, 3, 5, 1, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 4, 5, 1, 0, 5, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, 255, 255, 255, 255 },
	{ 9, 4, 5, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 11, 7, 4, 9, 11, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, 255, 255, 255, 255 },
	{ 1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, 255, 255, 255, 255 },
	{ 3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, 255 },
	{ 4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, 255, 255, 255, 255 },
	{ 9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, 255 },
	{ 11, 7, 4, 11, 4, 2, 2, 4, 0, 255, 255, 255, 255, 255, 255, 255 },
	{ 11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, 255, 255, 255, 255 },
	{ 2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, 255, 255, 255, 255 },
	{ 9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, 255 },
	{ 3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, 255 },
	{ 1, 10, 2, 8, 7, 4, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 9, 1, 4, 1, 7, 7, 1, 3, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, 255, 255, 255, 255 },
	{ 4, 0, 3, 7, 4, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 4, 8, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 10, 8, 10, 11, 8, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 0, 9, 3, 9, 11, 11, 9, 10, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 1, 10, 0, 10, 8, 8, 10, 11, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 1, 10, 11, 3, 10, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 2, 11, 1, 11, 9, 9, 11, 8, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, 255, 255, 255, 255 },
	{ 0, 2, 11, 8, 0, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 3, 2, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 2, 3, 8, 2, 8, 10, 10, 8, 9, 255, 255, 255, 255, 255, 255, 255 },
	{ 9, 10, 2, 0, 9, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, 255, 255, 255, 255 },
	{ 1, 10, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 1, 3, 8, 9, 1, 8, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 9, 1, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 0, 3, 8, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 },
	{ 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 }
};

enum class EdgeDirection {
	X = 0,
	Y = 1,
	Z = 2
};

const unsigned char edge2vertex[12][2] = {
	{0,1}, {1,2}, {3,2}, {0,3},
	{4,5}, {5,6}, {7,6}, {4,7},
	{0,4}, {1,5}, {2,6}, {3,7}
};

const EdgeDirection edge_directions[12] = {
	EdgeDirection::Y, EdgeDirection::X, EdgeDirection::Y, EdgeDirection::X,
	EdgeDirection::Y, EdgeDirection::X, EdgeDirection::Y, EdgeDirection::X,
	EdgeDirection::Z, EdgeDirection::Z, EdgeDirection::Z, EdgeDirection::Z
};

const bool cube_vertex_pos[8][3] = {
	{0,0,0}, {0,1,0}, {1,1,0}, {1,0,0},
	{0,0,1}, {0,1,1}, {1,1,1}, {1,0,1}
};

///////////////////////////////////////////////
//              Cubes convention             //
///////////////////////////////////////////////
//                  4_______e4_____________5
//                  /|                    /|
//                 / |                   / |
//              e7/  |                e5/  |
//               /___|______e6_________/   |
//             7|    |                 |6  |e9
//              |    |                 |   |
//              |    |e8            e10|   |
//           e11|    |                 |   |
//              |    |_________________|___|
//              |   / 0      e0        |   /1
//              |  /                   |  /
//              | /e3                  | /e1
//              |/_____________________|/
//              3         e2          2
//               ----> y
//			   /
//			  /
//           v
//			x

using namespace std;

struct TriangleDistanceIndex {
	unsigned int index;
	float distance_sq;
};

uintptr_t mallocFloatArray(unsigned int length) {
	if (length == 0) return 0;
	return reinterpret_cast<uintptr_t>(std::malloc(sizeof(float) * length));
}

uintptr_t mallocUIntArray(unsigned int length) {
	if (length == 0) return 0;
	return reinterpret_cast<uintptr_t>(std::malloc(sizeof(unsigned int) * length));
}

void freeArray(uintptr_t ptr) {
	if (ptr == 0) return;
	std::free(reinterpret_cast<void*>(ptr));
}

void buildTriangleDistancePermutation(uintptr_t positions_ptr,
		unsigned int vertex_count,
		float point_x,
		float point_y,
		float point_z,
		uintptr_t permutation_ptr) {
	if (positions_ptr == 0 || permutation_ptr == 0 || vertex_count < 3) return;

	const unsigned int triangle_count = vertex_count / 3;
	if (triangle_count == 0) return;

	const float* positions = reinterpret_cast<const float*>(positions_ptr);
	unsigned int* permutation = reinterpret_cast<unsigned int*>(permutation_ptr);
	vector<TriangleDistanceIndex> triangle_distances(triangle_count);

	for (unsigned int i = 0; i < triangle_count; ++i) {
		const unsigned int base = i * 9;
		const float center_x = (positions[base] + positions[base + 3] + positions[base + 6]) / 3.0f;
		const float center_y = (positions[base + 1] + positions[base + 4] + positions[base + 7]) / 3.0f;
		const float center_z = (positions[base + 2] + positions[base + 5] + positions[base + 8]) / 3.0f;
		const float dx = center_x - point_x;
		const float dy = center_y - point_y;
		const float dz = center_z - point_z;
		triangle_distances[i] = { i, dx * dx + dy * dy + dz * dz };
	}

	stable_sort(triangle_distances.begin(), triangle_distances.end(),
		[](const TriangleDistanceIndex& left, const TriangleDistanceIndex& right) {
			return left.distance_sq > right.distance_sq;
		});

	for (unsigned int i = 0; i < triangle_count; ++i) {
		permutation[i] = triangle_distances[i].index;
	}
}

void reorderArrayByPermutation(uintptr_t array_ptr,
		unsigned int vertex_count,
		unsigned int item_size,
		uintptr_t permutation_ptr) {
	if (array_ptr == 0 || permutation_ptr == 0 || vertex_count < 3 || item_size == 0) return;

	const unsigned int triangle_count = vertex_count / 3;
	if (triangle_count == 0) return;

	const unsigned int triangle_stride = item_size * 3;
	float* array = reinterpret_cast<float*>(array_ptr);
	const unsigned int* permutation = reinterpret_cast<const unsigned int*>(permutation_ptr);
	vector<float> temp(array, array + vertex_count * item_size);

	for (unsigned int i = 0; i < triangle_count; ++i) {
		const unsigned int src_base = permutation[i] * triangle_stride;
		const unsigned int dst_base = i * triangle_stride;
		for (unsigned int j = 0; j < triangle_stride; ++j) {
			array[dst_base + j] = temp[src_base + j];
		}
	}
}

inline float linear_interp(float isoval, float x1, float x2, float val1, float val2) {
	if (fabs(val1 - val2) < 0.000000000001f) return x1; // avoid division by zero
	return x1 + (x2 - x1) * (isoval - val1) / (val2 - val1); // interpolation
}

inline float linear_interp(float x1, float x2, float factor) {
	return x1 + (x2 - x1) * factor; // interpolation
}



// Isosurface smoothing method ids, shared with model/MarchingCubesWrapper.js.
// Mesh methods post-process the welded marching-cubes mesh; field methods
// replace the grid that is meshed (see apply_field_filter).
enum MeshSmoothing { MESH_OFF = 0, MESH_LAPLACIAN = 1, MESH_TAUBIN = 2, MESH_HC = 3, MESH_LOOP = 4, MESH_CATMULL_CLARK = 5 };
enum FieldFilter { FIELD_OFF = 0, FIELD_GAUSSIAN = 1, FIELD_TRICUBIC = 2 };

class MarchingCubes {
public:
	// The grid being MESHED. Equal to the source grid unless a field filter
	// (Gaussian blur, tricubic refinement) is active, in which case these
	// describe the filtered copy -- possibly with more points.
	unsigned int x_step, y_step, z_step;
	unsigned int nx, ny, nz;
	size_t field_size;
	size_t vertex_count;
	float dx, dy, dz;
	float* field;

	// The source grid as JS filled it (getField()). Always owned.
	unsigned int src_nx, src_ny, src_nz;
	size_t src_size;
	float* src_field;
	// Filtered copy of the field (and of the colour field for a refinement);
	// empty while no field filter is active.
	vector<float> work_field;
	vector<float> work_color;

	vector<float> vertex_list;
	vector<float> vnormal_list;
	vector<float> vnormal_cache; // 3 per grid point of the meshed grid, 0,0,0 = not computed yet

	// Optional secondary ("colour") field, e.g. sign(lambda2)*rho for an NCI
	// surface of s. Same grid as `field`; interpolated along each crossed edge
	// with the SAME `mu` as the vertex position (Jmol's process(vA, vB, f)), so
	// every output vertex carries the colour field's value at that point.
	// Both buffers exist only while enabled, so the default (single-field) path
	// pays no memory, and update_vertices dispatches to a template instance
	// without the extra work.
	float* src_color = nullptr;  // JS-filled, source grid
	float* color_field = nullptr; // meshed grid (src_color or work_color)
	vector<float> vvalue_list;

	// Optional mask: a cube with any corner value >= mask_value produces no
	// triangles. Jmol marks NCI grid points outside the density window as NaN
	// and skips every cube touching one; without the same rule the surface of s
	// closes off in artificial walls wherever the window ends (s jumps from its
	// real value to the excluded marker). Off (+inf) by default.
	float mask_value = INFINITY;
	bool has_mask = false;

	// Mesh smoothing settings (set_smoothing). The grid edge each soup vertex
	// was cut from is recorded while a method is active: it is the weld key.
	int mesh_method = MESH_OFF;
	int mesh_iterations = 0;
	float mesh_a = 0.0f, mesh_b = 0.0f;
	int mesh_levels = 0;
	bool mesh_project = false;
	vector<uint64_t> edge_keys;

	MarchingCubes(unsigned int resx, unsigned int resy, unsigned int resz)
			: vertex_count(0), src_nx(resx), src_ny(resy), src_nz(resz) {
		src_size = (size_t)resx * resy * resz;
		src_field = new float[src_size];
		use_source_grid();
	}

	~MarchingCubes() {
		delete[] src_field;
		// The colour buffer is always owned.
		delete[] src_color;
	}

	// Point the meshed grid at the source grid (no field filter).
	void use_source_grid() {
		work_field.clear(); work_field.shrink_to_fit();
		work_color.clear(); work_color.shrink_to_fit();
		set_mesh_grid(src_nx, src_ny, src_nz, src_field, src_color);
	}

	void set_mesh_grid(unsigned int resx, unsigned int resy, unsigned int resz, float* data, float* color) {
		nx = resx; ny = resy; nz = resz;
		x_step = 1;
		y_step = resx;
		z_step = resy * resx;
		dx = 1.0 / (float)(resx-1);
		dy = 1.0 / (float)(resy-1);
		dz = 1.0 / (float)(resz-1);
		field_size = (size_t)resx * resy * resz;
		field = data;
		color_field = color;
		// the gradients belong to the grid, so any change of grid resets them
		vnormal_cache.assign(field_size * 3, 0.0f);
	}

	/*
	 * Mesh a filtered copy of the source grid instead of the grid itself.
	 * Call again whenever JS rewrote getField() / getColorField(), since the
	 * copy is not kept in sync.
	 *   FIELD_OFF:      mesh the source grid
	 *   FIELD_GAUSSIAN: a = sigma in voxels. The colour field is left as is:
	 *                   it keeps reporting the true values at each vertex.
	 *   FIELD_TRICUBIC: factor = upsampling factor, a = kernel (0 Catmull-Rom,
	 *                   1 B-spline). The colour field is refined alongside.
	 * `periodic` wraps the filter stencils (Field.periodic).
	 */
	void apply_field_filter(int method, float a, int factor, bool periodic) {
		if (method == FIELD_GAUSSIAN && a > 0.0f) {
			work_color.clear(); work_color.shrink_to_fit();
			work_field.resize(src_size);
			iso_gaussian_filter(src_field, work_field.data(), src_nx, src_ny, src_nz,
				a, periodic, has_mask, mask_value);
			set_mesh_grid(src_nx, src_ny, src_nz, work_field.data(), src_color);
		} else if (method == FIELD_TRICUBIC && factor > 1) {
			const int kernel = a >= 0.5f ? 1 : 0;
			const unsigned int rx = iso_refined_dim(src_nx, factor);
			const unsigned int ry = iso_refined_dim(src_ny, factor);
			const unsigned int rz = iso_refined_dim(src_nz, factor);
			work_field.resize((size_t)rx * ry * rz);
			iso_refine_field(src_field, src_nx, src_ny, src_nz, factor, kernel,
				periodic, has_mask, mask_value, work_field.data());
			float* color = nullptr;
			if (src_color != nullptr) {
				work_color.resize(work_field.size());
				iso_refine_field(src_color, src_nx, src_ny, src_nz, factor, kernel,
					periodic, false, INFINITY, work_color.data());
				color = work_color.data();
			} else {
				work_color.clear(); work_color.shrink_to_fit();
			}
			set_mesh_grid(rx, ry, rz, work_field.data(), color);
		} else {
			use_source_grid();
		}
	}
	unsigned int get_mesh_nx() { return nx; }
	unsigned int get_mesh_ny() { return ny; }
	unsigned int get_mesh_nz() { return nz; }

	/*
	 * Mesh smoothing applied after every updateVertices():
	 *   MESH_LAPLACIAN:     iterations, a = lambda
	 *   MESH_TAUBIN:        iterations, a = lambda, b = pass-band k_PB
	 *   MESH_HC:            iterations, a = alpha, b = beta
	 *   MESH_LOOP:          levels
	 *   MESH_CATMULL_CLARK: levels
	 * `project` pulls the moved vertices back onto the isosurface afterwards.
	 */
	void set_smoothing(int method, int iterations, float a, float b, int levels, bool project) {
		mesh_method = (method >= MESH_LAPLACIAN && method <= MESH_CATMULL_CLARK) ? method : MESH_OFF;
		mesh_iterations = iterations;
		mesh_a = a;
		mesh_b = b;
		mesh_levels = levels;
		mesh_project = project;
		if (mesh_method == MESH_OFF) { edge_keys.clear(); edge_keys.shrink_to_fit(); }
	}

	/*
	 * Allocate (true) or free (false) the colour field. Idempotent. After
	 * enabling, JS fills getColorField() (one float per source grid point, same
	 * layout as getField()) and every updateVertices() also writes one
	 * interpolated value per vertex to getVertexValues().
	 */
	void set_color_field_enabled(bool enabled) {
		if (enabled == (src_color != nullptr)) return;
		if (enabled) {
			src_color = new float[src_size]();
		} else {
			delete[] src_color;
			src_color = nullptr;
			vvalue_list.clear(); vvalue_list.shrink_to_fit();
		}
		// a stale refined colour copy would outlive the change: back to the
		// source grid until JS re-applies its field filter
		use_source_grid();
	}
	bool get_color_field_enabled() {
		return src_color != nullptr;
	}
	uintptr_t get_color_field() {
		return reinterpret_cast<uintptr_t>(src_color);
	}
	uintptr_t get_vertex_values() {
		return reinterpret_cast<uintptr_t>(vvalue_list.data());
	}

	void set_mask_value(float value) {
		has_mask = std::isfinite(value);
		mask_value = has_mask ? value : INFINITY;
	}

	uintptr_t get_field() {
		return reinterpret_cast<uintptr_t>(src_field);
	}
	uintptr_t get_vertex_list() {
		return reinterpret_cast<uintptr_t>(vertex_list.data());
	}
	uintptr_t get_vnormal_list() {
		return reinterpret_cast<uintptr_t>(vnormal_list.data());
	}
	uintptr_t get_vnormal_cache() {
		return reinterpret_cast<uintptr_t>(vnormal_cache.data());
	}
	float get_vertex_count() {
		return vertex_count;
	}

	/*
	 * A starting isosurface level: the magnitude exceeded by `fraction` of the
	 * grid points.
	 *
	 * Picking the level by how much of the cell it encloses is what makes one
	 * rule work for a charge density, an ELF grid and a wavefunction alike --
	 * see model/CompositeField.js defaultIsoValue() for why the old midpoint of
	 * [min, max] does not. The JS side computed it by sampling 100k points and
	 * sorting them, which is a visible pause on every field of a large file.
	 * Here the field is already resident, so it costs one linear pass.
	 *
	 * The histogram bins on the FLOAT BIT PATTERN rather than on log10(|v|).
	 * Positive IEEE-754 floats compare in the same order as their bit patterns
	 * read as integers, so the top bits are already a log-spaced bucket index:
	 * the 8 exponent bits give the binade and the next 7 mantissa bits divide it
	 * into 128, which is ~0.8% relative resolution everywhere from 1e-38 to
	 * 1e38. That matters because a log10() per point cost more on a 120^3 grid
	 * (~40ms) than the sampling it was meant to replace; this way it is integer
	 * work only.
	 *
	 * Returns 0 for an empty or all-zero field.
	 */
	float default_isovalue(float fraction) {
		/* Keep sign+exponent+7 mantissa bits. Positive floats have the sign bit
		 * clear, so they occupy the low half of the index space. */
		const int SHIFT = 16;
		const int BIN_COUNT = 1 << (31 - SHIFT); /* 32768 */

		if (src_size == 0) return 0.0f;
		if (fraction < 0.001f) fraction = 0.001f;
		if (fraction > 0.999f) fraction = 0.999f;

		std::vector<size_t> bins(BIN_COUNT, 0);
		size_t counted = 0;
		for (size_t i = 0; i < src_size; i++) {
			const float a = std::fabs(src_field[i]);
			if (!(a > 0.0f) || !std::isfinite(a)) continue; /* zeros, NaN, inf */
			uint32_t bits;
			std::memcpy(&bits, &a, sizeof(bits));
			bins[bits >> SHIFT]++;
			counted++;
		}
		if (counted == 0) return 0.0f;

		/* Walk down from the largest values until `fraction` of the WHOLE grid
		 * is accounted for -- the zeros skipped above are part of the cell and
		 * have to count towards the enclosed volume, or a mostly-empty field
		 * would report a level enclosing far more than asked. */
		const double target = (double)fraction * (double)src_size;
		double running = 0.0;
		for (int bin = BIN_COUNT - 1; bin >= 0; bin--) {
			const double in_bin = (double)bins[bin];
			if (in_bin <= 0.0) continue;
			if (running + in_bin >= target) {
				/* Interpolate across the bin's bit range, assuming its points are
				 * spread evenly through it; without this the answer quantises to
				 * bin edges. Bit patterns are linear in the value within a
				 * binade, so this is a linear interpolation there. */
				const double needed = target - running;
				const double within = 1.0 - (needed / in_bin);
				uint32_t bits = ((uint32_t)bin << SHIFT)
					+ (uint32_t)(within * (double)((uint32_t)1 << SHIFT));
				float level;
				std::memcpy(&level, &bits, sizeof(level));
				return std::isfinite(level) ? level : 0.0f;
			}
			running += in_bin;
		}

		/* Fewer non-zero points than the requested fraction: the field is mostly
		 * zeros, so the smallest non-zero level is the honest answer. */
		for (int bin = 0; bin < BIN_COUNT; bin++) {
			if (bins[bin] > 0) {
				uint32_t bits = (uint32_t)bin << SHIFT;
				float level;
				std::memcpy(&level, &bits, sizeof(level));
				return level;
			}
		}
		return 0.0f;
	}


	inline void calc_norm(size_t v, float* norm) {
		// takes a grid point index, and 3-element array to write the normal to.
		// Computes the normal using central differences, with forward/backward
		// differences at the ends of each axis (tested per axis: at i=0 the flat
		// index v-1 would be the previous row's last point, not a neighbour).
		const size_t k = v / z_step;
		const size_t rem = v - k * z_step;
		const size_t j = rem / y_step;
		const size_t i = rem - j * y_step;
		norm[0] = axis_diff(v, x_step, i, nx, dx);
		norm[1] = axis_diff(v, y_step, j, ny, dy);
		norm[2] = axis_diff(v, z_step, k, nz, dz);
	}

	inline float axis_diff(size_t v, size_t step, size_t pos, unsigned int n, float h) const {
		if (n < 2) return 0.0f;
		if (pos == 0) return (field[v + step] - field[v]) / h;          // forward difference at start
		if (pos == n - 1) return (field[v] - field[v - step]) / h;     // backward difference at end
		return (field[v + step] - field[v - step]) / h * 0.5f;         // central difference
	}

	inline void _update_norm_cache(size_t v) {
		const size_t cache_ind = 3*v;
		if (vnormal_cache[cache_ind] == 0 && vnormal_cache[cache_ind + 1] == 0 && vnormal_cache[cache_ind + 2] == 0) {
			calc_norm(v, &vnormal_cache[cache_ind]);
		}
	}

	template <bool WithColor>
	inline void _add_edge_vertex(float isoval,
			unsigned char edge_ind,
			size_t v1_ind, size_t v2_ind,
			float x1, float y1, float z1,
			float field1, float field2,
			float* pos_out, float* normal_out, float* value_out) {

		float mu = 0;
		if (abs(field1 - field2) > 1e-16) { // avoid division by zero
			mu = (isoval - field1) / (field2 - field1);
		}

		// set position
		pos_out[0] = x1;
		pos_out[1] = y1;
		pos_out[2] = z1;
		switch (edge_directions[edge_ind]) {
			case EdgeDirection::X:
				pos_out[0] = linear_interp(x1, x1 + dx, mu);
				break;
			case EdgeDirection::Y:
				pos_out[1] = linear_interp(y1, y1 + dy, mu);
				break;
			case EdgeDirection::Z:
				pos_out[2] = linear_interp(z1, z1 + dz, mu);
				break;
		}

		// set normal via interpolation of the cached normals at the edge endpoints
		for (int i = 0; i < 3; i++) {
			normal_out[i] = linear_interp(vnormal_cache[3*v1_ind + i], vnormal_cache[3*v2_ind + i], mu);
		}

		// colour value: the colour field is read straight from the grid (no
		// per-point cache needed, unlike the normals), interpolated like the
		// position
		if (WithColor) {
			*value_out = linear_interp(color_field[v1_ind], color_field[v2_ind], mu);
		}
	}

	void update_vertices(float isovalue) {
		const bool record = mesh_method != MESH_OFF;
		if (color_field != nullptr) {
			if (record) _update_vertices<true, true>(isovalue);
			else _update_vertices<true, false>(isovalue);
		} else {
			if (record) _update_vertices<false, true>(isovalue);
			else _update_vertices<false, false>(isovalue);
		}
		if (record && vertex_count > 0) smooth_mesh(isovalue);
	}

	template <bool WithColor, bool RecordEdges>
	void _update_vertices(float isovalue) {
		float vertices_on_edge[12*3]; // flattened vertex positions for the 12 edges, each with x,y,z components
		float vnormals_on_edge[12*3]; // flattened vertex normals for the 12 edges, each with x,y,z components
		float values_on_edge[12]; // colour-field value per edge (WithColor only)
		uint64_t keys_on_edge[12]; // grid edge id per edge (RecordEdges only)
		int z_ind, y_ind, x_ind;
		size_t vertex_index = (size_t)(-1);
		this->vertex_count = 0;
		vertex_list.clear();
		vnormal_list.clear();
		vvalue_list.clear();
		edge_keys.clear();
		unsigned char cube_index;
		unsigned short edge_flags;
		size_t cube_points[8];

		if (nx < 2 || ny < 2 || nz < 2) {
			return;
		}

		for (z_ind = 0; z_ind < nz-1; z_ind++) {
			for (y_ind = 0; y_ind < ny-1; y_ind++) {
				for (x_ind = 0; x_ind < nx-1; x_ind++) {
					vertex_index = z_ind * z_step + y_ind * y_step + x_ind * x_step; // index of v0 vertex for current cube in flattened field array

					cube_points[0] = vertex_index; 								// v0
					cube_points[1] = vertex_index 		   + y_step			;	// v1
					cube_points[2] = vertex_index + x_step + y_step			; 	// v2
					cube_points[3] = vertex_index + x_step					; 	// v3
					cube_points[4] = vertex_index   				+ z_step; 	// v4
					cube_points[5] = vertex_index 		   + y_step + z_step; 	// v5
					cube_points[6] = vertex_index + x_step + y_step + z_step; 	// v6
					cube_points[7] = vertex_index + x_step 			+ z_step; 	// v7

					// map the isosurface encapsulation to byte-formatted table index
					cube_index = 0;
					if (field[cube_points[0]] < isovalue) cube_index |= 1;   // v0
					if (field[cube_points[1]] < isovalue) cube_index |= 2;   // v1
					if (field[cube_points[2]] < isovalue) cube_index |= 4;   // v2
					if (field[cube_points[3]] < isovalue) cube_index |= 8;   // v3
					if (field[cube_points[4]] < isovalue) cube_index |= 16;  // v4
					if (field[cube_points[5]] < isovalue) cube_index |= 32;  // v5
					if (field[cube_points[6]] < isovalue) cube_index |= 64;  // v6
					if (field[cube_points[7]] < isovalue) cube_index |= 128; // v7

					if (cube_index == 0 || cube_index == 255) continue; // skip empty cubes

					if (has_mask) {
						bool masked = false;
						for (int c = 0; c < 8; c++) {
							if (field[cube_points[c]] >= mask_value) { masked = true; break; }
						}
						if (masked) continue;
					}

					edge_flags = edgeTable[cube_index];

					// add new vertices and normals to list for current cube
					for (int i = 0; i < 12; i++) {
						if (edge_flags & (1 << i)) { // iterate through bit flags for each edge, and if set, compute vertex position and normal for that edge
							const size_t v1_ind = cube_points[edge2vertex[i][0]];
							const size_t v2_ind = cube_points[edge2vertex[i][1]];
							// pre-calculate normal vectors before interpolation, and cache them in vnormal_cache to avoid redundant calculations for shared vertices across edges
							_update_norm_cache(v1_ind);
							_update_norm_cache(v2_ind);
							// compute vertex position and normal for this edge
							// and store in vertex_list and vnormal_list at the appropriate index
							const float pos1[3] = {
								(x_ind + (int)cube_vertex_pos[edge2vertex[i][0]][0])*dx,
								(y_ind + (int)cube_vertex_pos[edge2vertex[i][0]][1])*dy,
								(z_ind + (int)cube_vertex_pos[edge2vertex[i][0]][2])*dz
							};
							const float field1 = field[v1_ind];
							const float field2 = field[v2_ind];
							_add_edge_vertex<WithColor>(isovalue,
								i, // edge index
								v1_ind, v2_ind, // vertex indices for the edge endpoints
								pos1[0], pos1[1], pos1[2], // starting position (assuming v1 and v2 are at the endpoints of same edge)
								field1, field2, // field values for interpolation
								&vertices_on_edge[i*3], // position output
								&vnormals_on_edge[i*3], // normal vector output
								&values_on_edge[i] // colour value output (WithColor only)
							);
							if (RecordEdges) {
								// a grid edge is its lower end point plus its axis;
								// the cubes sharing it all produce the same key
								keys_on_edge[i] = (uint64_t)std::min(v1_ind, v2_ind) * 3
									+ (uint64_t)edge_directions[i];
							}
						}
					}

					for (int i = 0; triTable[cube_index][i] != 255; i++) {
						const int edge = triTable[cube_index][i];
						vertex_list.push_back(vertices_on_edge[edge*3]);
						vertex_list.push_back(vertices_on_edge[edge*3 + 1]);
						vertex_list.push_back(vertices_on_edge[edge*3 + 2]);

						vnormal_list.push_back(vnormals_on_edge[edge*3]);
						vnormal_list.push_back(vnormals_on_edge[edge*3 + 1]);
						vnormal_list.push_back(vnormals_on_edge[edge*3 + 2]);
						if (WithColor) {
							vvalue_list.push_back(values_on_edge[edge]);
						}
						if (RecordEdges) {
							edge_keys.push_back(keys_on_edge[edge]);
						}
						this->vertex_count++;
					}
				}
			}
		}
	}

	/////////////////////////////////////
	// Mesh smoothing
	/////////////////////////////////////

	// Cell of the meshed grid containing normalised-box point p, and the
	// fractional position inside it. Points outside the box are clamped.
	inline void locate(const float* p, size_t* base, float* t) const {
		const unsigned int n[3] = {nx, ny, nz};
		const float h[3] = {dx, dy, dz};
		size_t idx[3];
		for (int a = 0; a < 3; a++) {
			float g = p[a] / h[a];
			const float gmax = (float)(n[a] - 1);
			if (!(g > 0.0f)) g = 0.0f;
			if (g > gmax) g = gmax;
			size_t i0 = (size_t)g;
			if (i0 >= n[a] - 1) i0 = n[a] - 2;
			idx[a] = i0;
			t[a] = g - (float)i0;
		}
		*base = idx[0] * x_step + idx[1] * y_step + idx[2] * z_step;
	}

	// Trilinear sample at normalised-box point p of the field value, its
	// gradient (from the cached central differences) and, when `color` is
	// given, the colour field.
	void sample(const float* p, float* value, float* grad, float* color) {
		size_t base;
		float t[3];
		locate(p, &base, t);
		*value = 0.0f;
		grad[0] = grad[1] = grad[2] = 0.0f;
		if (color) *color = 0.0f;
		for (int c = 0; c < 8; c++) {
			const int ox = c & 1, oy = (c >> 1) & 1, oz = (c >> 2) & 1;
			const size_t v = base + ox * x_step + oy * y_step + oz * z_step;
			const float w = (ox ? t[0] : 1.0f - t[0]) * (oy ? t[1] : 1.0f - t[1]) * (oz ? t[2] : 1.0f - t[2]);
			if (w == 0.0f) continue;
			_update_norm_cache(v);
			*value += w * field[v];
			grad[0] += w * vnormal_cache[3*v];
			grad[1] += w * vnormal_cache[3*v + 1];
			grad[2] += w * vnormal_cache[3*v + 2];
			if (color) *color += w * color_field[v];
		}
	}

	/*
	 * Post-process the triangle soup of the last update_vertices() with the
	 * active mesh method, then write it back as a soup (vertex_list etc.), so
	 * the JS side sees the same layout either way.
	 */
	void smooth_mesh(float isovalue) {
		IsoMesh m;
		weld(m);

		if (mesh_method == MESH_LOOP) {
			iso_subdivide_loop(m, mesh_levels);
		} else if (mesh_method == MESH_CATMULL_CLARK) {
			iso_subdivide_catmull_clark(m, mesh_levels);
		}
		const MeshTopology topo = build_topology(m);
		if (mesh_method == MESH_LAPLACIAN) {
			iso_smooth_laplacian(m, topo, mesh_iterations, mesh_a);
		} else if (mesh_method == MESH_TAUBIN) {
			iso_smooth_taubin(m, topo, mesh_iterations, mesh_a, iso_taubin_mu(mesh_a, mesh_b));
		} else if (mesh_method == MESH_HC) {
			iso_smooth_hc(m, topo, mesh_iterations, mesh_a, mesh_b);
		}

		const uint32_t nv = m.vertex_count();
		const bool with_color = color_field != nullptr;
		const float max_step = std::min(dx, std::min(dy, dz));
		vector<float> normals((size_t)nv * 3);
		vector<float> values(with_color ? nv : 0);
		for (uint32_t v = 0; v < nv; v++) {
			float* p = &m.pos[3*v];
			float f, g[3], c;
			sample(p, &f, g, with_color ? &c : nullptr);
			// Newton steps back onto f = iso along the gradient. Rim vertices
			// stay where marching cubes put them (on the isovalue already), so
			// the periodic copies keep meeting.
			if (mesh_project && !topo.boundary[v]) {
				for (int it = 0; it < 2; it++) {
					const float g2 = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
					if (!(g2 > 1e-30f)) break;
					const float s = (f - isovalue) / g2;
					float step[3] = {s * g[0], s * g[1], s * g[2]};
					const float len = std::sqrt(step[0]*step[0] + step[1]*step[1] + step[2]*step[2]);
					if (!std::isfinite(len)) break;
					if (len > max_step) {
						const float k = max_step / len;
						step[0] *= k; step[1] *= k; step[2] *= k;
					}
					for (int a = 0; a < 3; a++) p[a] = std::min(1.0f, std::max(0.0f, p[a] - step[a]));
					sample(p, &f, g, with_color ? &c : nullptr);
				}
			}
			normals[3*v] = g[0];
			normals[3*v + 1] = g[1];
			normals[3*v + 2] = g[2];
			if (with_color) values[v] = c;
		}

		// back to a soup
		const size_t nt = m.triangle_count();
		vertex_count = nt * 3;
		vertex_list.resize(vertex_count * 3);
		vnormal_list.resize(vertex_count * 3);
		vvalue_list.resize(with_color ? vertex_count : 0);
		for (size_t k = 0; k < vertex_count; k++) {
			const uint32_t v = m.tri[k];
			for (int a = 0; a < 3; a++) {
				vertex_list[3*k + a] = m.pos[3*v + a];
				vnormal_list[3*k + a] = normals[3*v + a];
			}
			if (with_color) vvalue_list[k] = values[v];
		}
	}

	// Soup -> indexed mesh: soup vertices cut from the same grid edge are
	// the same vertex.
	void weld(IsoMesh& m) {
		const size_t n = vertex_count;
		vector<uint32_t> order(n);
		for (size_t i = 0; i < n; i++) order[i] = (uint32_t)i;
		std::sort(order.begin(), order.end(), [this](uint32_t a, uint32_t b) {
			return edge_keys[a] < edge_keys[b];
		});
		vector<uint32_t> remap(n);
		m.pos.reserve(n / 3 + 16);
		for (size_t i = 0; i < n;) {
			const uint32_t id = m.vertex_count();
			const uint32_t first = order[i];
			m.pos.push_back(vertex_list[3*first]);
			m.pos.push_back(vertex_list[3*first + 1]);
			m.pos.push_back(vertex_list[3*first + 2]);
			size_t j = i;
			while (j < n && edge_keys[order[j]] == edge_keys[first]) remap[order[j++]] = id;
			i = j;
		}
		// Colour values are re-sampled at the final positions (smooth_mesh),
		// so the mesh does not need to carry them.
		m.has_val = false;
		m.tri.reserve(n);
		for (size_t t = 0; t + 2 < n; t += 3) {
			const uint32_t a = remap[t], b = remap[t + 1], c = remap[t + 2];
			if (a == b || b == c || c == a) continue;
			m.tri.push_back(a);
			m.tri.push_back(b);
			m.tri.push_back(c);
		}
	}
};

EMSCRIPTEN_BINDINGS(marching_cubes_module) {
	emscripten::function("mallocFloatArray", &mallocFloatArray);
	emscripten::function("mallocUIntArray", &mallocUIntArray);
	emscripten::function("freeArray", &freeArray);
	emscripten::function("buildTriangleDistancePermutation", &buildTriangleDistancePermutation);
	emscripten::function("reorderArrayByPermutation", &reorderArrayByPermutation);
	emscripten::class_<MarchingCubes>("MarchingCubes")
		.constructor<unsigned int, unsigned int, unsigned int>()
		.function("getField", &MarchingCubes::get_field)
		.function("getVertices", &MarchingCubes::get_vertex_list)
		.function("getNormals", &MarchingCubes::get_vnormal_list)
		.function("getVNormalCache", &MarchingCubes::get_vnormal_cache)
		.function("getVertexCount", &MarchingCubes::get_vertex_count)
		.function("defaultIsoValue", &MarchingCubes::default_isovalue)
		.function("setColorFieldEnabled", &MarchingCubes::set_color_field_enabled)
		.function("getColorFieldEnabled", &MarchingCubes::get_color_field_enabled)
		.function("getColorField", &MarchingCubes::get_color_field)
		.function("getVertexValues", &MarchingCubes::get_vertex_values)
		.function("setMaskValue", &MarchingCubes::set_mask_value)
		.function("applyFieldFilter", &MarchingCubes::apply_field_filter)
		.function("getMeshNx", &MarchingCubes::get_mesh_nx)
		.function("getMeshNy", &MarchingCubes::get_mesh_ny)
		.function("getMeshNz", &MarchingCubes::get_mesh_nz)
		.function("setSmoothing", &MarchingCubes::set_smoothing)
		.function("updateVertices", &MarchingCubes::update_vertices);
}

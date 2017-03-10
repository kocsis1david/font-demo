#pragma once

#include <stdint.h>
#include "linmath.h"
#include "geometry.h"

#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_BBOX_H 
#include FT_OUTLINE_H

#define FT_CHECK(r) { FT_Error err = (r); assert(!err); } while (0)

#define FD_OUTLINE_MAX_POINTS (255 * 2)

typedef struct fd_ContourRange
{
	uint32_t begin, end;
} fd_ContourRange;

typedef struct fd_Outline
{
	fd_Rect bbox;

	vec2 *points;
	uint32_t num_of_points;
	uint32_t point_capacity;

	fd_ContourRange *contours;
	uint32_t num_of_contours;
	uint32_t contour_capacity;

	uint32_t *cells;
	uint32_t cell_count_x;
	uint32_t cell_count_y;

	uint32_t corner_fix_begin;
} fd_Outline;

typedef struct fd_PointU16
{
	uint16_t x, y;
} fd_PointU16;

void fd_outline_convert(FT_Outline *outline, fd_Outline *o, char c);
void fd_outline_decompose(FT_Outline *outline, fd_Outline *o);
void fd_outline_make_cells(fd_Outline *o);
void fd_outline_subdivide(fd_Outline *o);
//void fd_outline_fix_corners(fd_Outline *o);
void fd_outline_destroy(fd_Outline *o);
void fd_outline_cbox(fd_Outline *o, fd_Rect *cbox);
void fd_outline_u16_points(fd_Outline *o, fd_Rect *cbox, fd_PointU16 *pout);
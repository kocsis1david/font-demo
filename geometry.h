#pragma once

#include <stdint.h>
#include "linmath.h"

typedef struct fd_Rect
{
	float min_x;
	float min_y;
	float max_x;
	float max_y;
} fd_Rect;

bool fd_bbox_bezier2_intersect(const fd_Rect *bbox, const vec2 bezier[3]);
float fd_line_signed_distance(const vec2 a, const vec2 b, const vec2 p);
float fd_line_calculate_t(const vec2 a, const vec2 b, const  vec2 p);
void fd_bezier2_point(vec2 r, const vec2 bezier[3], float t);
void fd_bezier2_split_lr(vec2 left[3], vec2 right[3], const vec2 bezier[3], float t);
void fd_bezier2_split_5p(vec2 ret[5], const vec2 bezier[3], float t);
void fd_bezier2_split_3p(vec2 ret[3], const vec2 bezier[3], float t);

void fd_bezier2_derivative(const vec2 bezier[3], vec2 derivative[2]);
void fd_bezier2_bbox(const vec2 bezier[3], fd_Rect *bbox);
void fd_bezier2_align_to_self(vec2 r[3], const vec2 bezier[3]);
void fd_bezier2_align_to_line(vec2 r[3], const vec2 bezier[3], const vec2 line0, const vec2 line1);
bool fd_bezier2_line_is_intersecting(const vec2 bezier[3], const vec2 line0, const vec2 line1);
#include "geometry.h"
#include <assert.h>
#include <stdlib.h>

enum
{
	FD_QUADRATIC_SOLUTION_NONE,
	FD_QUADRATIC_SOLUTION_ALL,
	FD_QUADRATIC_SOLUTION_TOUCH,
	FD_QUADRATIC_SOLUTION_ONE,
	FD_QUADRATIC_SOLUTION_TWO,
};

static uint32_t solve_quadratic(float a, float b, float c, float *x1, float *x2)
{
	float discriminant = b * b - 4.0f * a * c;

	if (discriminant > 0.0f)
	{
		float sqrt_d = sqrtf(discriminant);
		float common = b >= 0.0f ? -b - sqrt_d : -b + sqrt_d;

		*x1 = 2.0f * c / common;

		if (a == 0.0f)
			return FD_QUADRATIC_SOLUTION_ONE;

		*x2 = common / (2.0f * a);
		return FD_QUADRATIC_SOLUTION_TWO;
	}
	else if (discriminant == 0.0f)
	{
		if (b == 0.0f)
		{
			if (a == 0.0f)
			{
				if (c == 0.0f) return FD_QUADRATIC_SOLUTION_ALL;
				else return FD_QUADRATIC_SOLUTION_NONE;
			}
			else
			{
				*x1 = 0.0f;
				return FD_QUADRATIC_SOLUTION_TOUCH;
			}
		}
		else
		{
			*x1 = 2.0f * c / -b;
			return FD_QUADRATIC_SOLUTION_TOUCH;
		}
	}

	return FD_QUADRATIC_SOLUTION_NONE;
}

static float line_vertical_intersect(float x, const vec2 p1, const vec2 p2)
{
	float m = (p2[1] - p1[1]) / (p2[0] - p1[0]);
	return p1[1] - m * (p1[0] - x);
}

static float line_horizontal_intersect(float y, const vec2 p1, const vec2 p2)
{
	float n = (p2[0] - p1[0]) / (p2[1] - p1[1]);
	return p1[0] - n * (p1[1] - y);
}

static inline bool is_between(float value, float min, float max)
{
	return value >= min && value <= max;
}

static inline bool is_between_exclusive(float value, float min, float max)
{
	return value > min && value < max;
}

static inline bool is_point_inside_bbox(const fd_Rect *bbox, const vec2 p)
{
	return is_between(p[0], bbox->min_x, bbox->max_x) && is_between(p[1], bbox->min_y, bbox->max_y);
}

static inline bool is_point_inside_bbox_exclusive(const fd_Rect *bbox, const vec2 p)
{
	return is_between_exclusive(p[0], bbox->min_x, bbox->max_x) && is_between_exclusive(p[1], bbox->min_y, bbox->max_y);
}

static bool is_intersection_in_line_segment(const vec2 p1, const  vec2 p2, const vec2 i)
{
	float px_min = min(p1[0], p2[0]);
	float px_max = max(p1[0], p2[0]);
	float py_min = min(p1[1], p2[1]);
	float py_max = max(p1[1], p2[1]);

	return is_between(i[0], px_min, px_max) && is_between(i[1], py_min, py_max);
}

static bool is_line_segment_intersecting_bbox(const fd_Rect *bbox, const vec2 p1, const vec2 p2)
{
	if (is_point_inside_bbox_exclusive(bbox, p1)) return true;
	if (is_point_inside_bbox_exclusive(bbox, p2)) return true;

	float x_top = line_horizontal_intersect(bbox->max_y, p1, p2);
	float x_bottom = line_horizontal_intersect(bbox->min_y, p1, p2);
	float y_left = line_vertical_intersect(bbox->min_x, p1, p2);
	float y_right = line_vertical_intersect(bbox->max_x, p1, p2);

	vec2 top = { x_top, bbox->max_y };
	vec2 bottom = { x_bottom, bbox->min_y };
	vec2 left = { bbox->min_x, y_left };
	vec2 right = { bbox->max_x, y_right };

	if (is_between(x_top, bbox->min_x, bbox->max_x) &&
		is_intersection_in_line_segment(p1, p2, top))
	{
		return true;
	}

	if (is_between(x_bottom, bbox->min_x, bbox->max_x) &&
		is_intersection_in_line_segment(p1, p2, bottom))
	{
		return true;
	}

	if (is_between(y_left, bbox->min_y, bbox->max_y) &&
		is_intersection_in_line_segment(p1, p2, left))
	{
		return true;
	}

	if (is_between(y_right, bbox->min_y, bbox->max_y) &&
		is_intersection_in_line_segment(p1, p2, right))
	{
		return true;
	}

	return false;
}

bool fd_bbox_bezier2_intersect(const fd_Rect *bbox, const vec2 bezier[3])
{
	if (is_point_inside_bbox_exclusive(bbox, bezier[0])) return true;
	if (is_point_inside_bbox_exclusive(bbox, bezier[2])) return true;

	vec2 bl = { bbox->min_x, bbox->min_y };
	vec2 br = { bbox->max_x, bbox->min_y };
	vec2 tl = { bbox->min_x, bbox->max_y };
	vec2 tr = { bbox->max_x, bbox->max_y };

	return fd_bezier2_line_is_intersecting(bezier, bl, br) ||
		   fd_bezier2_line_is_intersecting(bezier, br, tr) ||
		   fd_bezier2_line_is_intersecting(bezier, tr, tl) ||
		   fd_bezier2_line_is_intersecting(bezier, tl, bl);
}

float fd_line_signed_distance(const vec2 a, const vec2 b, const vec2 p)
{
	vec2 line_dir;
	vec2_sub(line_dir, b, a);
	assert(vec2_len(line_dir) > 0.0f);

	vec2 perp_dir = { -line_dir[1], line_dir[0] };
	vec2_norm_(perp_dir);

	vec2 dir_to_a;
	vec2_sub(dir_to_a, a, p);

	return vec2_dot(perp_dir, dir_to_a);
}

float fd_line_calculate_t(const vec2 a, const vec2 b, const vec2 p)
{
	vec2 ab, ap;
	vec2_sub(ab, b, a);
	vec2_sub(ap, p, a);

	float t = vec2_dot(ap, ab) / vec2_dot(ab, ab);
	return clampf(t, 0.0, 1.0);
}

static inline void bezier2_points(vec2 q0, vec2 q1, vec2 r, const vec2 bezier[3], float t)
{
	vec2_lerp(q0, bezier[0], bezier[1], t);
	vec2_lerp(q1, bezier[1], bezier[2], t);
	vec2_lerp(r, q0, q1, t);
}

void fd_bezier2_point(vec2 r, const vec2 bezier[3], float t)
{
	vec2 q0, q1;
	bezier2_points(q0, q1, r, bezier, t);
}

void fd_bezier2_split_lr(vec2 left[3], vec2 right[3], const vec2 bezier[3], float t)
{
	vec2 q0, q1, r;
	bezier2_points(q0, q1, r, bezier, t);

	vec2_dup(left[0], bezier[0]);
	vec2_dup(left[1], q0);
	vec2_dup(left[2], r);

	vec2_dup(right[0], r);
	vec2_dup(right[1], q1);
	vec2_dup(right[2], bezier[2]);
}

void fd_bezier2_split_5p(vec2 ret[5], const vec2 bezier[3], float t)
{
	vec2 q0, q1, r;
	bezier2_points(q0, q1, r, bezier, t);

	vec2_dup(ret[0], bezier[0]);
	vec2_dup(ret[1], q0);
	vec2_dup(ret[2], r);
	vec2_dup(ret[3], q1);
	vec2_dup(ret[4], bezier[2]);
}

void fd_bezier2_split_3p(vec2 ret[3], const vec2 bezier[3], float t)
{
	vec2 q0, q1, r;
	bezier2_points(q0, q1, r, bezier, t);

	vec2_dup(ret[0], q0);
	vec2_dup(ret[1], r);
	vec2_dup(ret[2], q1);
}

void fd_bezier2_derivative(const vec2 bezier[3], vec2 derivative[2])
{
	vec2_sub(derivative[0], bezier[1], bezier[0]);
	vec2_scale_(derivative[0], 2);

	vec2_sub(derivative[1], bezier[2], bezier[1]);
	vec2_scale_(derivative[1], 2);
}

static inline float bezier2_component(float p0, float p1, float p2, float t)
{
	return lerpf(lerpf(p0, p1, t), lerpf(p1, p2, t), t);
}

void fd_bezier2_bbox(const vec2 bezier[3], fd_Rect *bbox)
{
	vec2 deriv[2];
	fd_bezier2_derivative(bezier, deriv);

	float tx = deriv[0][0] / (deriv[0][0] - deriv[1][0]);
	float ty = deriv[0][1] / (deriv[0][1] - deriv[1][1]);

	bbox->min_x = min(bezier[0][0], bezier[2][0]);
	bbox->min_y = min(bezier[0][1], bezier[2][1]);
	bbox->max_x = max(bezier[0][0], bezier[2][0]);
	bbox->max_y = max(bezier[0][1], bezier[2][1]);

	if (0.0f <= tx && tx <= 1.0f)
	{
		float x = bezier2_component(bezier[0][0], bezier[1][0], bezier[2][0], tx);

		if (deriv[0][0] < deriv[1][0])
			bbox->min_x = min(bbox->min_x, x);
		else bbox->max_x = max(bbox->max_x, x);
	}

	if (0.0f <= ty && ty <= 1.0f)
	{
		float y = bezier2_component(bezier[0][1], bezier[1][1], bezier[2][1], ty);

		if (deriv[0][1] < deriv[1][1])
			bbox->min_y = min(bbox->min_y, y);
		else bbox->max_y = max(bbox->max_y, y);
	}
}

static inline void align_point(vec2 r, const vec2 p, const vec2 t, float s, float c)
{
	vec2 tmp;
	vec2_sub(tmp, p, t);

	r[0] = tmp[0] * c - tmp[1] * s;
	r[1] = tmp[0] * s + tmp[1] * c;
}

static inline void align_lsc(const vec2 p0, const vec2 p1, float *l, float *s, float *c)
{
	vec2 v;
	vec2_sub(v, p1, p0);

	*l = vec2_len(v);
	*s = -v[1] / *l;
	*c = v[0] / *l;
}

void fd_bezier2_align_to_self(vec2 r[3], const vec2 bezier[3])
{
	float l, s, c;
	align_lsc(bezier[0], bezier[2], &l, &s, &c);

	vec2_set(r[0], 0.0f, 0.0f);
	align_point(r[1], bezier[1], bezier[0], s, c);
	vec2_set(r[2], l, 0.0f);
}

void fd_bezier2_align_to_line(vec2 r[3], const vec2 bezier[3], const vec2 line0, const vec2 line1)
{
	float l, s, c;
	align_lsc(line0, line1, &l, &s, &c);

	align_point(r[0], bezier[0], line0, s, c);
	align_point(r[1], bezier[1], line0, s, c);
	align_point(r[2], bezier[2], line0, s, c);
}

bool fd_bezier2_line_is_intersecting(const vec2 bezier[3], const vec2 line0, const vec2 line1)
{
	float l, si, co;
	align_lsc(line0, line1, &l, &si, &co);

	vec2 bez[3];
	align_point(bez[0], bezier[0], line0, si, co);
	align_point(bez[1], bezier[1], line0, si, co);
	align_point(bez[2], bezier[2], line0, si, co);

	float x0 = bez[0][0], y0 = bez[0][1];
	float x1 = bez[1][0], y1 = bez[1][1];
	float x2 = bez[2][0], y2 = bez[2][1];

	float a = y0 - 2 * y1 + y2;
	float b = 2 * (y1 - y0);
	float c = y0;

	float t0, t1, xt0, xt1;
	uint32_t sol = solve_quadratic(a, b, c, &t0, &t1);

	switch (sol)
	{
	case FD_QUADRATIC_SOLUTION_NONE:
	case FD_QUADRATIC_SOLUTION_ALL:
		return false;

	case FD_QUADRATIC_SOLUTION_TOUCH:
	case FD_QUADRATIC_SOLUTION_ONE:
		xt0 = bezier2_component(x0, x1, x2, t0);
		return is_between(t0, 0, 1) && is_between(xt0, 0, l);

	case FD_QUADRATIC_SOLUTION_TWO:
		xt0 = bezier2_component(x0, x1, x2, t0);
		xt1 = bezier2_component(x0, x1, x2, t1);

		return is_between(t0, 0, 1) && is_between(xt0, 0, l) ||
			   is_between(t1, 0, 1) && is_between(xt1, 0, l);

	default:
		assert(false);
		return false;
	}
}
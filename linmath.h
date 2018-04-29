#ifndef LINMATH_H
#define LINMATH_H

#include <math.h>
#include <stdbool.h>

#ifndef max
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#endif
#ifndef min
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

static inline float clampf(float x, float min, float max)
{
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

static inline float lerpf(float x, float y, float t)
{
    return x * (1 - t) + y * t;
}

#define LINMATH_H_DEFINE_VEC(n) \
typedef float vec##n[n]; \
static inline void vec##n##_dup(vec##n r, vec##n const a) \
{ \
    for (int i = 0; i < n; ++i) \
        r[i] = a[i]; \
} \
static inline bool vec##n##_eq(vec##n const a, vec##n const b) \
{ \
    for (int i = 0; i < n; ++i) \
        if (b[i] != a[i]) return false; \
    return true; \
} \
static inline void vec##n##_add(vec##n r, vec##n const a, vec##n const b) \
{ \
    for (int i = 0; i < n; ++i) \
        r[i] = a[i] + b[i]; \
} \
static inline void vec##n##_add_(vec##n r, vec##n const a) \
{ \
    for (int i = 0; i < n; ++i) \
        r[i] += a[i]; \
} \
static inline void vec##n##_sub(vec##n r, vec##n const a, vec##n const b) \
{ \
    for (int i = 0; i < n; ++i) \
        r[i] = a[i] - b[i]; \
} \
static inline void vec##n##_sub_(vec##n r, vec##n const b) \
{ \
    for (int i = 0; i < n; ++i) \
        r[i] -= b[i]; \
} \
static inline void vec##n##_scale(vec##n r, vec##n const v, float const s) \
{ \
    for (int i = 0; i < n; ++i) \
        r[i] = v[i] * s; \
} \
static inline void vec##n##_scale_(vec##n r, float const s) \
{ \
    for (int i = 0; i < n; ++i) \
        r[i] *= s; \
} \
static inline float vec##n##_dot(vec##n const a, vec##n const b) \
{ \
    float p = 0.; \
    for(int i = 0; i < n; ++i) \
        p += b[i] * a[i]; \
    return p; \
} \
static inline float vec##n##_len(vec##n const v) \
{ \
    return sqrtf(vec##n##_dot(v,v)); \
} \
static inline void vec##n##_norm(vec##n r, vec##n const v) \
{ \
    float k = 1.0f / vec##n##_len(v); \
    vec##n##_scale(r, v, k); \
} \
static inline void vec##n##_norm_(vec##n r) \
{ \
    float k = 1.0f / vec##n##_len(r); \
    vec##n##_scale_(r, k); \
} \
static inline void vec##n##_min(vec##n r, const vec##n a, const vec##n b) \
{ \
    for (int i = 0; i < n; ++i) \
        r[i] = a[i]<b[i] ? a[i] : b[i]; \
} \
static inline void vec##n##_max(vec##n r, const vec##n a, const vec##n b) \
{ \
    for (int i = 0; i < n; ++i) \
        r[i] = a[i]>b[i] ? a[i] : b[i]; \
} \
static inline void vec##n##_lerp(vec##n r, const vec##n a, const vec##n b, float x) \
{ \
    vec##n ta, tb; \
    vec##n##_scale(ta, a, 1 - x); \
    vec##n##_scale(tb, b, x); \
    vec##n##_add(r, ta, tb); \
}

LINMATH_H_DEFINE_VEC(2)
LINMATH_H_DEFINE_VEC(3)
LINMATH_H_DEFINE_VEC(4)

static inline float vec2_dist(const vec2 a, const vec2 b)
{
    vec2 tmp;
    vec2_sub(tmp, a, b);
    return vec2_len(tmp);
}

static inline void vec2_swap(vec2 a, vec2 b)
{
    vec2 tmp;
    vec2_dup(tmp, a);
    vec2_dup(a, b);
    vec2_dup(b, tmp);
}

static inline void vec3_proj(vec3 r, vec3 a, vec3 b)
{
    vec3_scale(r, a, vec3_dot(a, b) / vec3_dot(a, a));
}

static inline void vec2_set(vec2 r, float x, float y)
{
    r[0] = x;
    r[1] = y;
}

static inline void vec3_set(vec3 r, float x, float y, float z)
{
    r[0] = x;
    r[1] = y;
    r[2] = z;
}

static inline void vec3_cross(vec3 r, vec3 const a, vec3 const b)
{
    r[0] = a[1]*b[2] - a[2]*b[1];
    r[1] = a[2]*b[0] - a[0]*b[2];
    r[2] = a[0]*b[1] - a[1]*b[0];
}

static inline void vec3_reflect(vec3 r, vec3 const v, vec3 const n)
{
    float p = 2.f * vec3_dot(v, n);
    for (int i = 0; i < 3; ++i)
        r[i] = v[i] - p*n[i];
}

static inline void vec4_cross(vec4 r, vec4 a, vec4 b)
{
    r[0] = a[1]*b[2] - a[2]*b[1];
    r[1] = a[2]*b[0] - a[0]*b[2];
    r[2] = a[0]*b[1] - a[1]*b[0];
    r[3] = 1.f;
}

static inline void vec4_reflect(vec4 r, vec4 v, vec4 n)
{
    float p = 2.f * vec4_dot(v, n);
    for (int i = 0; i < 4; ++i)
        r[i] = v[i] - p*n[i];
}

typedef vec4 mat4x4[4];
static inline void mat4x4_set(mat4x4 m,
    float m00, float m10, float m20, float m30,
    float m01, float m11, float m21, float m31,
    float m02, float m12, float m22, float m32,
    float m03, float m13, float m23, float m33)
{
    m[0][0] = m00;
    m[0][1] = m01;
    m[0][2] = m02;
    m[0][3] = m03;

    m[1][0] = m10;
    m[1][1] = m11;
    m[1][2] = m12;
    m[1][3] = m13;

    m[2][0] = m20;
    m[2][1] = m21;
    m[2][2] = m22;
    m[2][3] = m23;

    m[3][0] = m30;
    m[3][1] = m31;
    m[3][2] = m32;
    m[3][3] = m33;
}
/*
static inline void mat4x4_set_vec4(mat4x4 m, vec4 m0, vec4 m1, vec4 m2, vec4 m3)
{
    vec4_dup(m[0], m0);
    vec4_dup(m[1], m1);
    vec4_dup(m[2], m2);
    vec4_dup(m[3], m3);
}

static inline void mat4x4_set_vec3(mat4x4 m, vec3 m0, vec3 m1, vec3 m2)
{
    vec4_dup(m[0], m0);
    vec4_dup(m[1], m1);
    vec4_dup(m[2], m2);
    vec4_dup(m[3], m3);
}
*/
static inline void mat4x4_identity(mat4x4 M)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            M[i][j] = i == j ? 1.f : 0.f;
}
static inline void mat4x4_dup(mat4x4 M, mat4x4 N)
{
    for (int i = 0; i < 4; ++i)
        vec4_dup(M[i], N[i]);
}
static inline void mat4x4_row(vec4 r, mat4x4 M, int i)
{
    for(int k = 0; k < 4; ++k)
        r[k] = M[k][i];
}
static inline void mat4x4_col(vec4 r, mat4x4 M, int i)
{
    vec4_dup(r, M[i]);
}
static inline void mat4x4_transpose(mat4x4 M, mat4x4 N)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            M[i][j] = N[j][i];
}
static inline void mat4x4_add(mat4x4 M, mat4x4 a, mat4x4 b)
{
    for (int i = 0; i < 4; ++i)
        vec4_add(M[i], a[i], b[i]);
}
static inline void mat4x4_sub(mat4x4 M, mat4x4 a, mat4x4 b)
{
    for (int i = 0; i < 4; ++i)
        vec4_sub(M[i], a[i], b[i]);
}
static inline void mat4x4_scale(mat4x4 M, mat4x4 a, float k)
{
    for (int i = 0; i < 4; ++i)
        vec4_scale(M[i], a[i], k);
}
static inline void mat4x4_scale_aniso(mat4x4 M, mat4x4 a, float x, float y, float z)
{
    vec4_scale(M[0], a[0], x);
    vec4_scale(M[1], a[1], y);
    vec4_scale(M[2], a[2], z);

    for (int i = 0; i < 4; ++i)
        M[3][i] = a[3][i];
}
static inline void mat4x4_mul(mat4x4 M, mat4x4 a, mat4x4 b)
{
    mat4x4 temp;
    int k, r, c;
    for(c=0; c<4; ++c) for(r=0; r<4; ++r) {
        temp[c][r] = 0.f;
        for(k=0; k<4; ++k)
            temp[c][r] += a[k][r] * b[c][k];
    }
    mat4x4_dup(M, temp);
}
static inline void mat4x4_mul_vec4(vec4 r, mat4x4 M, vec4 v)
{
    vec4 temp;

    int i, j;
    for(j=0; j<4; ++j) {
        temp[j] = 0.f;
        for(i=0; i<4; ++i)
            temp[j] += M[i][j] * v[i];
    }

    vec4_dup(r, temp);
}
static inline void mat4x4_translation(mat4x4 T, float x, float y, float z)
{
    mat4x4_identity(T);
    T[3][0] = x;
    T[3][1] = y;
    T[3][2] = z;
}
static inline void mat4x4_translation_vec(mat4x4 T, vec3 v)
{
    mat4x4_translation(T, v[0], v[1], v[2]);
}
static inline void mat4x4_from_vec3_mul_outer(mat4x4 M, vec3 a, vec3 b)
{
    int i, j;
    for(i=0; i<4; ++i) for(j=0; j<4; ++j)
        M[i][j] = i<3 && j<3 ? a[i] * b[j] : 0.f;
}
static inline void mat4x4_rotate(mat4x4 R, mat4x4 M, float x, float y, float z, float angle)
{
    float s = sinf(angle);
    float c = cosf(angle);
    vec3 u = {x, y, z};

    if(vec3_len(u) > 1e-4) {
        vec3_norm(u, u);
        mat4x4 T;
        mat4x4_from_vec3_mul_outer(T, u, u);

        mat4x4 S = {
            {    0,  u[2], -u[1], 0},
            {-u[2],     0,  u[0], 0},
            { u[1], -u[0],     0, 0},
            {    0,     0,     0, 0}
        };
        mat4x4_scale(S, S, s);

        mat4x4 C;
        mat4x4_identity(C);
        mat4x4_sub(C, C, T);

        mat4x4_scale(C, C, c);

        mat4x4_add(T, T, C);
        mat4x4_add(T, T, S);

        T[3][3] = 1.;
        mat4x4_mul(R, M, T);
    } else {
        mat4x4_dup(R, M);
    }
}
static inline void mat4x4_rotate_X(mat4x4 Q, mat4x4 M, float angle)
{
    float s = sinf(angle);
    float c = cosf(angle);
    mat4x4 R = {
        {1.f, 0.f, 0.f, 0.f},
        {0.f,   c,   s, 0.f},
        {0.f,  -s,   c, 0.f},
        {0.f, 0.f, 0.f, 1.f}
    };
    mat4x4_mul(Q, M, R);
}
static inline void mat4x4_rotate_Y(mat4x4 Q, mat4x4 M, float angle)
{
    float s = sinf(angle);
    float c = cosf(angle);
    mat4x4 R = {
        {   c, 0.f,   s, 0.f},
        { 0.f, 1.f, 0.f, 0.f},
        {  -s, 0.f,   c, 0.f},
        { 0.f, 0.f, 0.f, 1.f}
    };
    mat4x4_mul(Q, M, R);
}
static inline void mat4x4_rotate_Z(mat4x4 Q, mat4x4 M, float angle)
{
    float s = sinf(angle);
    float c = cosf(angle);
    mat4x4 R = {
        {   c,   s, 0.f, 0.f},
        {  -s,   c, 0.f, 0.f},
        { 0.f, 0.f, 1.f, 0.f},
        { 0.f, 0.f, 0.f, 1.f}
    };
    mat4x4_mul(Q, M, R);
}
static inline void mat4x4_rotatation_z(mat4x4 Q, float angle)
{
    float s = sinf(angle);
    float c = cosf(angle);
    mat4x4 R = {
        { c,   s, 0.f, 0.f },
        { -s,   c, 0.f, 0.f },
        { 0.f, 0.f, 1.f, 0.f },
        { 0.f, 0.f, 0.f, 1.f }
    };
    mat4x4_dup(Q, R);
}
static inline void mat4x4_invert(mat4x4 T, mat4x4 M)
{
    float s[6];
    float c[6];
    s[0] = M[0][0]*M[1][1] - M[1][0]*M[0][1];
    s[1] = M[0][0]*M[1][2] - M[1][0]*M[0][2];
    s[2] = M[0][0]*M[1][3] - M[1][0]*M[0][3];
    s[3] = M[0][1]*M[1][2] - M[1][1]*M[0][2];
    s[4] = M[0][1]*M[1][3] - M[1][1]*M[0][3];
    s[5] = M[0][2]*M[1][3] - M[1][2]*M[0][3];

    c[0] = M[2][0]*M[3][1] - M[3][0]*M[2][1];
    c[1] = M[2][0]*M[3][2] - M[3][0]*M[2][2];
    c[2] = M[2][0]*M[3][3] - M[3][0]*M[2][3];
    c[3] = M[2][1]*M[3][2] - M[3][1]*M[2][2];
    c[4] = M[2][1]*M[3][3] - M[3][1]*M[2][3];
    c[5] = M[2][2]*M[3][3] - M[3][2]*M[2][3];

    /* Assumes it is invertible */
    float idet = 1.0f/( s[0]*c[5]-s[1]*c[4]+s[2]*c[3]+s[3]*c[2]-s[4]*c[1]+s[5]*c[0] );

    T[0][0] = ( M[1][1] * c[5] - M[1][2] * c[4] + M[1][3] * c[3]) * idet;
    T[0][1] = (-M[0][1] * c[5] + M[0][2] * c[4] - M[0][3] * c[3]) * idet;
    T[0][2] = ( M[3][1] * s[5] - M[3][2] * s[4] + M[3][3] * s[3]) * idet;
    T[0][3] = (-M[2][1] * s[5] + M[2][2] * s[4] - M[2][3] * s[3]) * idet;

    T[1][0] = (-M[1][0] * c[5] + M[1][2] * c[2] - M[1][3] * c[1]) * idet;
    T[1][1] = ( M[0][0] * c[5] - M[0][2] * c[2] + M[0][3] * c[1]) * idet;
    T[1][2] = (-M[3][0] * s[5] + M[3][2] * s[2] - M[3][3] * s[1]) * idet;
    T[1][3] = ( M[2][0] * s[5] - M[2][2] * s[2] + M[2][3] * s[1]) * idet;

    T[2][0] = ( M[1][0] * c[4] - M[1][1] * c[2] + M[1][3] * c[0]) * idet;
    T[2][1] = (-M[0][0] * c[4] + M[0][1] * c[2] - M[0][3] * c[0]) * idet;
    T[2][2] = ( M[3][0] * s[4] - M[3][1] * s[2] + M[3][3] * s[0]) * idet;
    T[2][3] = (-M[2][0] * s[4] + M[2][1] * s[2] - M[2][3] * s[0]) * idet;

    T[3][0] = (-M[1][0] * c[3] + M[1][1] * c[1] - M[1][2] * c[0]) * idet;
    T[3][1] = ( M[0][0] * c[3] - M[0][1] * c[1] + M[0][2] * c[0]) * idet;
    T[3][2] = (-M[3][0] * s[3] + M[3][1] * s[1] - M[3][2] * s[0]) * idet;
    T[3][3] = ( M[2][0] * s[3] - M[2][1] * s[1] + M[2][2] * s[0]) * idet;
}
static inline void mat4x4_orthonormalize(mat4x4 R, mat4x4 M)
{
    mat4x4_dup(R, M);
    float s = 1.;
    vec3 h;

    vec3_norm(R[2], R[2]);

    s = vec3_dot(R[1], R[2]);
    vec3_scale(h, R[2], s);
    vec3_sub(R[1], R[1], h);
    vec3_norm(R[2], R[2]);

    s = vec3_dot(R[1], R[2]);
    vec3_scale(h, R[2], s);
    vec3_sub(R[1], R[1], h);
    vec3_norm(R[1], R[1]);

    s = vec3_dot(R[0], R[1]);
    vec3_scale(h, R[1], s);
    vec3_sub(R[0], R[0], h);
    vec3_norm(R[0], R[0]);
}
/*
it doesn't use depth 0 .. 1
static inline void mat4x4_frustum(mat4x4 M, float l, float r, float b, float t, float n, float f)
{
    M[0][0] = 2.f*n/(r-l);
    M[0][1] = M[0][2] = M[0][3] = 0.f;

    M[1][1] = 2.0f*n/(t-b);
    M[1][0] = M[1][2] = M[1][3] = 0.f;

    M[2][0] = (r+l)/(r-l);
    M[2][1] = (t+b)/(t-b);
    M[2][2] = -(f+n)/(f-n);
    M[2][3] = -1.f;

    M[3][2] = -2.f*(f*n)/(f-n);
    M[3][0] = M[3][1] = M[3][3] = 0.f;
}
*/
static inline void mat4x4_ortho(mat4x4 m, float w, float h, float n, float f)
{
    mat4x4_set(m,
        2.0f / w,      0.0f,          0.0f,       0.0f,
            0.0f, -2.0f / h,          0.0f,       0.0f,
            0.0f,      0.0f, -1.0f / (f-n), -n / (f-n),
            //0, 0, 0, 0.5f,
            0.0f,      0.0f,          0.0f,       1.0f);
}
static inline void mat4x4_perspective(mat4x4 m, float y_fov, float aspect, float n, float f)
{
    /* NOTE: Degrees are an unhandy unit to work with.
     * linmath.h uses radians for everything! */
    float const a = 1.f / tanf(y_fov / 2.f);

    mat4x4_set(m,
        a / aspect, 0.0f,      0.0f,           0.0f,
        0.0f,         -a,      0.0f,           0.0f,
        0.0f,       0.0f, f / (n-f), -(f*n) / (f-n),
        0.0f,       0.0f,     -1.0f,           0.0f);
}
static inline void mat4x4_look_at(mat4x4 m, vec3 eye, vec3 center, vec3 up)
{
    vec3 z;
    vec3_sub(z, eye, center);
    vec3_norm_(z);

    vec3 x;
    vec3_cross(x, up, z);
    vec3_norm_(x);

    vec3 y;
    vec3_cross(y, z, x);

    mat4x4_set(m,
        x[0], x[1], x[2], -vec3_dot(eye, x),
        y[0], y[1], y[2], -vec3_dot(eye, y),
        z[0], z[1], z[2], -vec3_dot(eye, z),
        0.0f, 0.0f, 0.0f, 1.0f);
}

typedef float quat[4];
static inline void quat_identity(quat q)
{
    q[0] = q[1] = q[2] = 0.f;
    q[3] = 1.f;
}
static inline void quat_add(quat r, quat a, quat b)
{
    for (int i = 0; i < 4; ++i)
        r[i] = a[i] + b[i];
}
static inline void quat_sub(quat r, quat a, quat b)
{
    for (int i = 0; i < 4; ++i)
        r[i] = a[i] - b[i];
}
static inline void quat_mul(quat r, quat p, quat q)
{
    vec3 w;
    vec3_cross(r, p, q);
    vec3_scale(w, p, q[3]);
    vec3_add(r, r, w);
    vec3_scale(w, q, p[3]);
    vec3_add(r, r, w);
    r[3] = p[3]*q[3] - vec3_dot(p, q);
}
static inline void quat_scale(quat r, quat v, float s)
{
    for (int i = 0; i < 4; ++i)
        r[i] = v[i] * s;
}
static inline float quat_inner_product(quat a, quat b)
{
    float p = 0.f;
    for (int i = 0; i < 4; ++i)
        p += b[i]*a[i];
    return p;
}
static inline void quat_conj(quat r, quat q)
{
    for (int i = 0; i < 3; ++i)
        r[i] = -q[i];

    r[3] = q[3];
}
static inline void quat_rotate(quat r, float angle, vec3 axis)
{
    vec3 v;
    vec3_scale(v, axis, sinf(angle / 2));
    for (int i = 0; i < 3; ++i)
        r[i] = v[i];
    r[3] = cosf(angle / 2);
}
static inline void quat_lerp(quat r, quat q1, quat q2, float x)
{
    vec4_lerp(r, q1, q2, x);
    vec4_norm_(r);
}
#define quat_norm vec4_norm
static inline void quat_mul_vec3(vec3 r, quat q, vec3 v)
{
/*
 * Method by Fabian 'ryg' Giessen (of Farbrausch)
t = 2 * cross(q.xyz, v)
v' = v + q.w * t + cross(q.xyz, t)
 */
    vec3 t;
    vec3 q_xyz = {q[0], q[1], q[2]};
    vec3 u = {q[0], q[1], q[2]};

    vec3_cross(t, q_xyz, v);
    vec3_scale(t, t, 2);

    vec3_cross(u, q_xyz, t);
    vec3_scale(t, t, q[3]);

    vec3_add(r, v, t);
    vec3_add(r, r, u);
}
static inline void mat4x4_from_quat(mat4x4 M, quat q)
{
    float a = q[3];
    float b = q[0];
    float c = q[1];
    float d = q[2];
    float a2 = a*a;
    float b2 = b*b;
    float c2 = c*c;
    float d2 = d*d;

    M[0][0] = a2 + b2 - c2 - d2;
    M[0][1] = 2.f*(b*c + a*d);
    M[0][2] = 2.f*(b*d - a*c);
    M[0][3] = 0.f;

    M[1][0] = 2*(b*c - a*d);
    M[1][1] = a2 - b2 + c2 - d2;
    M[1][2] = 2.f*(c*d + a*b);
    M[1][3] = 0.f;

    M[2][0] = 2.f*(b*d + a*c);
    M[2][1] = 2.f*(c*d - a*b);
    M[2][2] = a2 - b2 - c2 + d2;
    M[2][3] = 0.f;

    M[3][0] = M[3][1] = M[3][2] = 0.f;
    M[3][3] = 1.f;
}

static inline void mat4x4o_mul_quat(mat4x4 R, mat4x4 M, quat q)
{
/*  XXX: The way this is written only works for othogonal matrices. */
/* TODO: Take care of non-orthogonal case. */
    quat_mul_vec3(R[0], q, M[0]);
    quat_mul_vec3(R[1], q, M[1]);
    quat_mul_vec3(R[2], q, M[2]);

    R[3][0] = R[3][1] = R[3][2] = 0.f;
    R[3][3] = 1.f;
}
static inline void quat_from_mat4x4(quat q, mat4x4 M)
{
    float r=0.f;
    int i;

    int perm[] = { 0, 1, 2, 0, 1 };
    int *p = perm;

    for(i = 0; i<3; i++) {
        float m = M[i][i];
        if( m < r )
            continue;
        m = r;
        p = &perm[i];
    }

    r = sqrtf(1.f + M[p[0]][p[0]] - M[p[1]][p[1]] - M[p[2]][p[2]] );

    if(r < 1e-6) {
        q[0] = 1.f;
        q[1] = q[2] = q[3] = 0.f;
        return;
    }

    q[0] = r/2.f;
    q[1] = (M[p[0]][p[1]] - M[p[1]][p[0]])/(2.f*r);
    q[2] = (M[p[2]][p[0]] - M[p[0]][p[2]])/(2.f*r);
    q[3] = (M[p[2]][p[1]] - M[p[1]][p[2]])/(2.f*r);
}

static inline void mat4x4_lerp(mat4x4 R, mat4x4 m1, mat4x4 m2, float x)
{
    mat4x4 tmp;
    mat4x4_scale(R, m1, 1 - x);
    mat4x4_scale(tmp, m2, x);
    mat4x4_add(R, R, tmp);
}

static inline void mat4x4_interpolate(mat4x4 R, mat4x4 m1, mat4x4 m2, float x)
{
    vec3 s1 = { vec3_len(m1[0]), vec3_len(m1[1]), vec3_len(m1[2]) };
    vec3 s2 = { vec3_len(m2[0]), vec3_len(m2[1]), vec3_len(m2[2]) };

    quat r1, r2;
    quat_from_mat4x4(r1, m1);
    quat_from_mat4x4(r2, m2);

    vec3 t, s;
    quat r;

    vec3_lerp(t, m1[3], m2[3], x);
    vec3_lerp(s, s1, s2, x);
    quat_lerp(r, r1, r2, x);

    mat4x4_from_quat(R, r);
    mat4x4_scale_aniso(R, R, s[0], s[1], s[2]);
    vec3_dup(R[3], t);
}
#endif

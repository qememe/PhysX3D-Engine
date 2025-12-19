#pragma once

#include <cmath>
#include <iostream>
#include <array>
#include <immintrin.h>

namespace physx3d {

struct alignas(16) Vec3 {
    float x, y, z;
    float _pad;
    
    constexpr Vec3() : x(0), y(0), z(0), _pad(0) {}
    constexpr Vec3(float x, float y, float z) : x(x), y(y), z(z), _pad(0) {}
    
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(float s) const { return Vec3(x * s, y * s, z * s); }
    Vec3 operator/(float s) const { return Vec3(x / s, y / s, z / s); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    
    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    Vec3& operator*=(float s) { x *= s; y *= s; z *= s; return *this; }
    
    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    
    Vec3 cross(const Vec3& v) const {
        return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    
    float lengthSq() const { return x * x + y * y + z * z; }
    float length() const { return std::sqrt(lengthSq()); }
    
    Vec3 normalized() const {
        float len = length();
        return len > 1e-6f ? (*this / len) : Vec3(0, 0, 0);
    }
    
    void normalize() {
        float len = length();
        if (len > 1e-6f) { x /= len; y /= len; z /= len; }
    }
    
    float& operator[](size_t i) { return (&x)[i]; }
    const float& operator[](size_t i) const { return (&x)[i]; }
    
    bool isZero() const { return lengthSq() < 1e-8f; }
};

inline Vec3 operator*(float s, const Vec3& v) { return v * s; }

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

struct Mat3 {
    std::array<float, 9> m;
    
    Mat3() : m{1,0,0, 0,1,0, 0,0,1} {}
    
    Mat3(float m00, float m01, float m02,
         float m10, float m11, float m12,
         float m20, float m21, float m22)
        : m{m00, m01, m02, m10, m11, m12, m20, m21, m22} {}
    
    static Mat3 identity() { return Mat3(); }
    static Mat3 scale(float s) { return Mat3(s, 0, 0, 0, s, 0, 0, 0, s); }
    static Mat3 diagonal(float x, float y, float z) {
        return Mat3(x, 0, 0, 0, y, 0, 0, 0, z);
    }
    
    Vec3 operator*(const Vec3& v) const {
        return Vec3(
            m[0] * v.x + m[1] * v.y + m[2] * v.z,
            m[3] * v.x + m[4] * v.y + m[5] * v.z,
            m[6] * v.x + m[7] * v.y + m[8] * v.z
        );
    }
    
    Mat3 operator*(const Mat3& other) const {
        Mat3 result;
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                result.m[row * 3 + col] = 
                    m[row * 3 + 0] * other.m[0 * 3 + col] +
                    m[row * 3 + 1] * other.m[1 * 3 + col] +
                    m[row * 3 + 2] * other.m[2 * 3 + col];
            }
        }
        return result;
    }
    
    Mat3 transpose() const {
        return Mat3(m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]);
    }
    
    float determinant() const {
        return m[0] * (m[4] * m[8] - m[5] * m[7]) -
               m[1] * (m[3] * m[8] - m[5] * m[6]) +
               m[2] * (m[3] * m[7] - m[4] * m[6]);
    }
    
    Mat3 inverse() const {
        float det = determinant();
        if (std::abs(det) < 1e-8f) return Mat3::identity();
        
        float invDet = 1.0f / det;
        return Mat3(
            (m[4] * m[8] - m[5] * m[7]) * invDet,
            (m[2] * m[7] - m[1] * m[8]) * invDet,
            (m[1] * m[5] - m[2] * m[4]) * invDet,
            (m[5] * m[6] - m[3] * m[8]) * invDet,
            (m[0] * m[8] - m[2] * m[6]) * invDet,
            (m[2] * m[3] - m[0] * m[5]) * invDet,
            (m[3] * m[7] - m[4] * m[6]) * invDet,
            (m[1] * m[6] - m[0] * m[7]) * invDet,
            (m[0] * m[4] - m[1] * m[3]) * invDet
        );
    }
    
    float& operator()(int row, int col) { return m[row * 3 + col]; }
    const float& operator()(int row, int col) const { return m[row * 3 + col]; }
};

struct Quaternion {
    float w, x, y, z;
    
    Quaternion() : w(1), x(0), y(0), z(0) {}
    Quaternion(float w, float x, float y, float z) : w(w), x(x), y(y), z(z) {}
    
    static Quaternion identity() { return Quaternion(); }
    
    static Quaternion fromAxisAngle(const Vec3& axis, float angle) {
        float halfAngle = angle * 0.5f;
        float s = std::sin(halfAngle);
        Vec3 normalizedAxis = axis.normalized();
        return Quaternion(
            std::cos(halfAngle),
            normalizedAxis.x * s,
            normalizedAxis.y * s,
            normalizedAxis.z * s
        );
    }
    
    Quaternion operator*(const Quaternion& q) const {
        return Quaternion(
            w * q.w - x * q.x - y * q.y - z * q.z,
            w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x,
            w * q.z + x * q.y - y * q.x + z * q.w
        );
    }
    
    Quaternion operator*(float s) const {
        return Quaternion(w * s, x * s, y * s, z * s);
    }
    
    Quaternion& operator*=(const Quaternion& q) {
        *this = *this * q;
        return *this;
    }
    
    Quaternion operator+(const Quaternion& q) const {
        return Quaternion(w + q.w, x + q.x, y + q.y, z + q.z);
    }
    
    Quaternion& operator+=(const Quaternion& q) {
        w += q.w; x += q.x; y += q.y; z += q.z;
        return *this;
    }
    
    float lengthSq() const { return w * w + x * x + y * y + z * z; }
    float length() const { return std::sqrt(lengthSq()); }
    
    void normalize() {
        float len = length();
        if (len > 1e-6f) {
            w /= len; x /= len; y /= len; z /= len;
        }
    }
    
    Quaternion normalized() const {
        Quaternion q = *this;
        q.normalize();
        return q;
    }
    
    Quaternion conjugate() const {
        return Quaternion(w, -x, -y, -z);
    }
    
    Vec3 rotate(const Vec3& v) const {
        Quaternion qv(0, v.x, v.y, v.z);
        Quaternion result = (*this) * qv * conjugate();
        return Vec3(result.x, result.y, result.z);
    }
    
    Mat3 toMat3() const {
        float xx = x * x, yy = y * y, zz = z * z;
        float xy = x * y, xz = x * z, yz = y * z;
        float wx = w * x, wy = w * y, wz = w * z;
        
        return Mat3(
            1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy),
            2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx),
            2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy)
        );
    }
};

inline float clamp(float value, float min, float max) {
    return value < min ? min : (value > max ? max : value);
}

inline Vec3 clamp(const Vec3& v, float min, float max) {
    return Vec3(clamp(v.x, min, max), clamp(v.y, min, max), clamp(v.z, min, max));
}

}

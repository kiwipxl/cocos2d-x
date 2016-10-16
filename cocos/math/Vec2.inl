/**
 Copyright 2013 BlackBerry Inc.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 Original file from GamePlay3D: http://gameplay3d.org

 This file was modified to fit the cocos2d-x project
 */

#include "math/Vec2.h"

NS_CC_MATH_BEGIN

template <class T> inline bool Vector2<T>::isZero() const
{
    return x == 0.0f && y == 0.0f;
}

template <class T> bool Vector2<T>::isOne() const
{
    return x == 1.0f && y == 1.0f;
}

template <class T> inline void Vector2<T>::add(const Vector2<T>& v)
{
    x += v.x;
    y += v.y;
}

template <class T> inline T Vector2<T>::distanceSquared(const Vector2<T>& v) const
{
    T dx = v.x - x;
    T dy = v.y - y;
    return (dx * dx + dy * dy);
}

template <class T> inline T Vector2<T>::dot(const Vector2<T>& v) const
{
    return (x * v.x + y * v.y);
}

template <class T> inline T Vector2<T>::lengthSquared() const
{
    return (x * x + y * y);
}

template <class T> inline void Vector2<T>::negate()
{
    x = -x;
    y = -y;
}

template <class T> inline void Vector2<T>::scale(T scalar)
{
    x *= scalar;
    y *= scalar;
}

template <class T> inline void Vector2<T>::scale(const Vector2<T>& scale)
{
    x *= scale.x;
    y *= scale.y;
}

template <class T> inline void Vector2<T>::set(T xx, T yy)
{
    this->x = xx;
    this->y = yy;
}

template <class T> inline void Vector2<T>::set(const Vector2<T>& v)
{
    this->x = v.x;
    this->y = v.y;
}

template <class T> inline void Vector2<T>::set(const Vector2<T>& p1, const Vector2<T>& p2)
{
    x = p2.x - p1.x;
    y = p2.y - p1.y;
}

template <class T> void Vector2<T>::setZero()
{
    x = y = 0.0f;
}

template <class T> inline void Vector2<T>::subtract(const Vector2<T>& v)
{
    x -= v.x;
    y -= v.y;
}

template <class T> inline void Vector2<T>::smooth(const Vector2<T>& target, float elapsedTime, float responseTime)
{
    if (elapsedTime > 0)
    {
        *this += (target - *this) * (elapsedTime / (elapsedTime + responseTime));
    }
}

template <class T> inline const Vector2<T> Vector2<T>::operator+(const Vector2<T>& v) const
{
    Vector2<T> result(*this);
    result.add(v);
    return result;
}

template <class T> inline Vector2<T>& Vector2<T>::operator+=(const Vector2<T>& v)
{
    add(v);
    return *this;
}

template <class T> inline const Vector2<T> Vector2<T>::operator-(const Vector2<T>& v) const
{
    Vector2<T> result(*this);
    result.subtract(v);
    return result;
}

template <class T> inline Vector2<T>& Vector2<T>::operator-=(const Vector2<T>& v)
{
    subtract(v);
    return *this;
}

template <class T> inline const Vector2<T> Vector2<T>::operator-() const
{
    Vector2<T> result(*this);
    result.negate();
    return result;
}

template <class T> inline const Vector2<T> Vector2<T>::operator*(T s) const
{
    Vector2<T> result(*this);
    result.scale(s);
    return result;
}

template <class T> inline Vector2<T>& Vector2<T>::operator*=(T s)
{
    scale(s);
    return *this;
}

template <class T> inline const Vector2<T> Vector2<T>::operator/(const T s) const
{
    return Vector2<T>(this->x / s, this->y / s);
}

template <class T> inline bool Vector2<T>::operator<(const Vector2<T>& v) const
{
    if (x == v.x)
    {
        return y < v.y;
    }
    return x < v.x;
}

template <class T> inline bool Vector2<T>::operator>(const Vector2<T>& v) const
{
    if (x == v.x)
    {
        return y > v.y;
    }
    return x > v.x;
}

template <class T> inline bool Vector2<T>::operator==(const Vector2<T>& v) const
{
    return x==v.x && y==v.y;
}

template <class T> inline bool Vector2<T>::operator!=(const Vector2<T>& v) const
{
    return x!=v.x || y!=v.y;
}

template <class T> inline const Vector2<T> operator*(T x, const Vector2<T>& v)
{
    Vector2<T> result(v);
    result.scale(x);
    return result;
}

template <class T> void Vector2<T>::setPoint(T xx, T yy)
{
    this->x = xx;
    this->y = yy;
}

// below is previously from Vec2.cpp


// returns true if segment A-B intersects with segment C-D. S->E is the overlap part
template<class T> bool isOneDimensionSegmentOverlap(T A, T B, T C, T D, T *S, T * E)
{
    T ABmin = std::min(A, B);
    T ABmax = std::max(A, B);
    T CDmin = std::min(C, D);
    T CDmax = std::max(C, D);

    if (ABmax < CDmin || CDmax < ABmin)
    {
        // ABmin->ABmax->CDmin->CDmax or CDmin->CDmax->ABmin->ABmax
        return false;
    }
    else
    {
        if (ABmin >= CDmin && ABmin <= CDmax)
        {
            // CDmin->ABmin->CDmax->ABmax or CDmin->ABmin->ABmax->CDmax
            if (S != nullptr) *S = ABmin;
            if (E != nullptr) *E = CDmax < ABmax ? CDmax : ABmax;
        }
        else if (ABmax >= CDmin && ABmax <= CDmax)
        {
            // ABmin->CDmin->ABmax->CDmax
            if (S != nullptr) *S = CDmin;
            if (E != nullptr) *E = ABmax;
        }
        else
        {
            // ABmin->CDmin->CDmax->ABmax
            if (S != nullptr) *S = CDmin;
            if (E != nullptr) *E = CDmax;
        }
        return true;
    }
}

// cross product of 2 vector. A->B X C->D
template<class T> T crossProduct2Vector(const Vector2<T>& A, const Vector2<T>& B, const Vector2<T>& C, const Vector2<T>& D)
{
    return (D.y - C.y) * (B.x - A.x) - (D.x - C.x) * (B.y - A.y);
}

template<class T> T Vector2<T>::angle(const Vector2<T>& v1, const Vector2<T>& v2)
{
    T dz = v1.x * v2.y - v1.y * v2.x;
    return atan2f(fabsf(dz) + MATH_FLOAT_SMALL, dot(v1, v2));
}

template<class T> void Vector2<T>::add(const Vector2<T>& v1, const Vector2<T>& v2, Vector2<T>* dst)
{
    GP_ASSERT(dst);

    dst->x = v1.x + v2.x;
    dst->y = v1.y + v2.y;
}

template<class T> void Vector2<T>::clamp(const Vector2<T>& min, const Vector2<T>& max)
{
    GP_ASSERT(!(min.x > max.x || min.y > max.y));

    // Clamp the x value.
    if (x < min.x)
        x = min.x;
    if (x > max.x)
        x = max.x;

    // Clamp the y value.
    if (y < min.y)
        y = min.y;
    if (y > max.y)
        y = max.y;
}

template<class T> void Vector2<T>::clamp(const Vector2<T>& v, const Vector2<T>& min, const Vector2<T>& max, Vector2<T>* dst)
{
    GP_ASSERT(dst);
    GP_ASSERT(!(min.x > max.x || min.y > max.y));

    // Clamp the x value.
    dst->x = v.x;
    if (dst->x < min.x)
        dst->x = min.x;
    if (dst->x > max.x)
        dst->x = max.x;

    // Clamp the y value.
    dst->y = v.y;
    if (dst->y < min.y)
        dst->y = min.y;
    if (dst->y > max.y)
        dst->y = max.y;
}

template<class T> T Vector2<T>::distance(const Vector2<T>& v) const
{
    T dx = v.x - x;
    T dy = v.y - y;

    return std::sqrt(dx * dx + dy * dy);
}

template<class T> T Vector2<T>::dot(const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x * v2.x + v1.y * v2.y);
}

template<class T> T Vector2<T>::length() const
{
    return std::sqrt(x * x + y * y);
}

template<class T> void Vector2<T>::normalize()
{
    T n = x * x + y * y;
    // Already normalized.
    if (n == 1.0f)
        return;

    n = std::sqrt(n);
    // Too close to zero.
    if (n < MATH_TOLERANCE)
        return;

    n = 1.0f / n;
    x *= n;
    y *= n;
}

template<class T> Vector2<T> Vector2<T>::getNormalized() const
{
    Vector2<T> v(*this);
    v.normalize();
    return v;
}

template<class T> void Vector2<T>::rotate(const Vector2<T>& point, T angle)
{
    T sinAngle = std::sin(angle);
    T cosAngle = std::cos(angle);

    if (point.isZero())
    {
        T tempX = x * cosAngle - y * sinAngle;
        y = y * cosAngle + x * sinAngle;
        x = tempX;
    }
    else
    {
        T tempX = x - point.x;
        T tempY = y - point.y;

        x = tempX * cosAngle - tempY * sinAngle + point.x;
        y = tempY * cosAngle + tempX * sinAngle + point.y;
    }
}

template<class T> void Vector2<T>::set(const T* array)
{
    GP_ASSERT(array);

    x = array[0];
    y = array[1];
}

template<class T> void Vector2<T>::subtract(const Vector2<T>& v1, const Vector2<T>& v2, Vector2<T>* dst)
{
    GP_ASSERT(dst);

    dst->x = v1.x - v2.x;
    dst->y = v1.y - v2.y;
}

template<class T> bool Vector2<T>::equals(const Vector2<T>& target) const
{
    return (std::abs(this->x - target.x) < FLT_EPSILON)
        && (std::abs(this->y - target.y) < FLT_EPSILON);
}

template<class T> bool Vector2<T>::fuzzyEquals(const Vector2<T>& b, T var) const
{
    if (x - var <= b.x && b.x <= x + var)
        if (y - var <= b.y && b.y <= y + var)
            return true;
    return false;
}

template<class T> T Vector2<T>::getAngle(const Vector2<T>& other) const
{
    Vector2<T> a2 = getNormalized();
    Vector2<T> b2 = other.getNormalized();
    T angle = atan2f(a2.cross(b2), a2.dot(b2));
    if (std::abs(angle) < FLT_EPSILON) return 0.f;
    return angle;
}

template<class T> Vector2<T> Vector2<T>::rotateByAngle(const Vector2<T>& pivot, T angle) const
{
    return pivot + (*this - pivot).rotate(Vector2<T>::forAngle(angle));
}

template<class T> bool Vector2<T>::isLineIntersect(const Vector2<T>& A, const Vector2<T>& B,
    const Vector2<T>& C, const Vector2<T>& D,
    T *s, T *t)
{
    // FAIL: Line undefined
    if ((A.x == B.x && A.y == B.y) || (C.x == D.x && C.y == D.y))
    {
        return false;
    }

    const T denom = crossProduct2Vector(A, B, C, D);

    if (denom == 0)
    {
        // Lines parallel or overlap
        return false;
    }

    if (s != nullptr) *s = crossProduct2Vector(C, D, C, A) / denom;
    if (t != nullptr) *t = crossProduct2Vector(A, B, C, A) / denom;

    return true;
}

template<class T> bool Vector2<T>::isLineParallel(const Vector2<T>& A, const Vector2<T>& B,
    const Vector2<T>& C, const Vector2<T>& D)
{
    // FAIL: Line undefined
    if ((A.x == B.x && A.y == B.y) || (C.x == D.x && C.y == D.y))
    {
        return false;
    }

    if (crossProduct2Vector(A, B, C, D) == 0)
    {
        // line overlap
        if (crossProduct2Vector(C, D, C, A) == 0 || crossProduct2Vector(A, B, C, A) == 0)
        {
            return false;
        }

        return true;
    }

    return false;
}

template<class T> bool Vector2<T>::isLineOverlap(const Vector2<T>& A, const Vector2<T>& B,
    const Vector2<T>& C, const Vector2<T>& D)
{
    // FAIL: Line undefined
    if ((A.x == B.x && A.y == B.y) || (C.x == D.x && C.y == D.y))
    {
        return false;
    }

    if (crossProduct2Vector(A, B, C, D) == 0 &&
        (crossProduct2Vector(C, D, C, A) == 0 || crossProduct2Vector(A, B, C, A) == 0))
    {
        return true;
    }

    return false;
}

template<class T> bool Vector2<T>::isSegmentOverlap(const Vector2<T>& A, const Vector2<T>& B, const Vector2<T>& C, const Vector2<T>& D, Vector2<T>* S, Vector2<T>* E)
{

    if (isLineOverlap(A, B, C, D))
    {
        return isOneDimensionSegmentOverlap(A.x, B.x, C.x, D.x, &S->x, &E->x) &&
            isOneDimensionSegmentOverlap(A.y, B.y, C.y, D.y, &S->y, &E->y);
    }

    return false;
}

template<class T> bool Vector2<T>::isSegmentIntersect(const Vector2<T>& A, const Vector2<T>& B, const Vector2<T>& C, const Vector2<T>& D)
{
    T s, t;

    if (isLineIntersect(A, B, C, D, &s, &t) &&
        (s >= 0.0f && s <= 1.0f && t >= 0.0f && t <= 1.0f))
    {
        return true;
    }

    return false;
}

template<class T> Vector2<T> Vector2<T>::getIntersectPoint(const Vector2<T>& A, const Vector2<T>& B, const Vector2<T>& C, const Vector2<T>& D)
{
    T s, t;

    if (isLineIntersect(A, B, C, D, &s, &t))
    {
        // Vector2<T> of intersection
        Vector2<T> P;
        P.x = A.x + s * (B.x - A.x);
        P.y = A.y + s * (B.y - A.y);
        return P;
    }

    return Vector2<T>::ZERO;
}

template<class T> const Vector2<T> Vector2<T>::ZERO = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::ONE = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::UNIT_X = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::UNIT_Y = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::ANCHOR_MIDDLE = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::ANCHOR_BOTTOM_LEFT = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::ANCHOR_TOP_LEFT = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::ANCHOR_BOTTOM_RIGHT = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::ANCHOR_TOP_RIGHT = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::ANCHOR_MIDDLE_RIGHT = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::ANCHOR_MIDDLE_LEFT = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::ANCHOR_MIDDLE_TOP = Vector2<T>();
template<class T> const Vector2<T> Vector2<T>::ANCHOR_MIDDLE_BOTTOM = Vector2<T>();

template<> const Vector2<float> Vector2<float>::ZERO = Vector2<float>(0.0f, 0.0f);
template<> const Vector2<float> Vector2<float>::ONE = Vector2<float>(1.0f, 1.0f);
template<> const Vector2<float> Vector2<float>::UNIT_X = Vector2<float>(1.0f, 0.0f);
template<> const Vector2<float> Vector2<float>::UNIT_Y = Vector2<float>(0.0f, 1.0f);
template<> const Vector2<float> Vector2<float>::ANCHOR_MIDDLE = Vector2<float>(0.5f, 0.5f);
template<> const Vector2<float> Vector2<float>::ANCHOR_BOTTOM_LEFT = Vector2<float>(0.0f, 0.0f);
template<> const Vector2<float> Vector2<float>::ANCHOR_TOP_LEFT = Vector2<float>(0.0f, 1.0f);
template<> const Vector2<float> Vector2<float>::ANCHOR_BOTTOM_RIGHT = Vector2<float>(1.0f, 0.0f);
template<> const Vector2<float> Vector2<float>::ANCHOR_TOP_RIGHT = Vector2<float>(1.0f, 1.0f);
template<> const Vector2<float> Vector2<float>::ANCHOR_MIDDLE_RIGHT = Vector2<float>(1.0f, 0.5f);
template<> const Vector2<float> Vector2<float>::ANCHOR_MIDDLE_LEFT = Vector2<float>(0.0f, 0.5f);
template<> const Vector2<float> Vector2<float>::ANCHOR_MIDDLE_TOP = Vector2<float>(0.5f, 1.0f);
template<> const Vector2<float> Vector2<float>::ANCHOR_MIDDLE_BOTTOM = Vector2<float>(0.5f, 0.0f);

NS_CC_MATH_END

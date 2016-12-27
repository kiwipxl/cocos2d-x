/**
 Copyright 2013 BlackBerry Inc.
 Copyright (c) 2014-2015 Chukong Technologies

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

#ifndef MATH_VEC2_H
#define MATH_VEC2_H

#include <algorithm>
#include <functional>
#include <cmath>
#include "math/CCMathBase.h"

/**
 * @addtogroup base
 * @{
 */

NS_CC_MATH_BEGIN

/** Clamp a value between from and to.
 */

inline float clampf(float value, float min_inclusive, float max_inclusive)
{
    if (min_inclusive > max_inclusive) {
        std::swap(min_inclusive, max_inclusive);
    }
    return value < min_inclusive ? min_inclusive : value < max_inclusive? value : max_inclusive;
}

class Mat4;

/**
 * Defines a 2-element floating point vector.
 */
template <class T>
class Vector2
{
public:

    /**
     * The x coordinate.
     */
    T x;

    /**
     * The y coordinate.
     */
    T y;

    /**
     * Constructs a new vector initialized to all zeros.
     */
    Vector2() : x(0.0f), y(0.0f) { }

    /**
     * Constructs a new vector initialized to the specified values.
     *
     * @param xx The x coordinate.
     * @param yy The y coordinate.
     */
    Vector2(T xx, T yy) : x(xx), y(yy) { }

    /**
     * Constructs a new vector from the values in the specified array.
     *
     * @param array An array containing the elements of the vector in the order x, y.
     */
    Vector2(const T* array) {
        set(array);
    }

    /**
     * Constructs a vector that describes the direction between the specified points.
     *
     * @param p1 The first point.
     * @param p2 The second point.
     */
    Vector2(const Vector2<T>& p1, const Vector2<T>& p2) {
        set(p1, p2);
    }

    /**
     * Constructs a new vector that is a copy of the specified vector.
     *
     * @param copy The vector to copy.
     */
    Vector2(const Vector2<T>& copy) {
        set(copy);
    }

    /**
     * Destructor.
     */
    ~Vector2() { }

    /**
     * Indicates whether this vector contains all zeros.
     *
     * @return true if this vector contains all zeros, false otherwise.
     */
    inline bool isZero() const;

    /**
     * Indicates whether this vector contains all ones.
     *
     * @return true if this vector contains all ones, false otherwise.
     */
    inline bool isOne() const;

    /**
     * Returns the angle (in radians) between the specified vectors.
     *
     * @param v1 The first vector.
     * @param v2 The second vector.
     * 
     * @return The angle between the two vectors (in radians).
     */
    static T angle(const Vector2<T>& v1, const Vector2<T>& v2);

    /**
     * Adds the elements of the specified vector to this one.
     *
     * @param v The vector to add.
     */
    inline void add(const Vector2<T>& v);

    /**
     * Adds the specified vectors and stores the result in dst.
     *
     * @param v1 The first vector.
     * @param v2 The second vector.
     * @param dst A vector to store the result in.
     */
    static void add(const Vector2<T>& v1, const Vector2<T>& v2, Vector2<T>* dst);

    /**
     * Clamps this vector within the specified range.
     *
     * @param min The minimum value.
     * @param max The maximum value.
     */
    void clamp(const Vector2<T>& min, const Vector2<T>& max);

    /**
     * Clamps the specified vector within the specified range and returns it in dst.
     *
     * @param v The vector to clamp.
     * @param min The minimum value.
     * @param max The maximum value.
     * @param dst A vector to store the result in.
     */
    static void clamp(const Vector2<T>& v, const Vector2<T>& min, const Vector2<T>& max, Vector2<T>* dst);

    /**
     * Returns the distance between this vector and v.
     *
     * @param v The other vector.
     * 
     * @return The distance between this vector and v.
     * 
     * @see distanceSquared
     */
    T distance(const Vector2<T>& v) const;

    /**
     * Returns the squared distance between this vector and v.
     *
     * When it is not necessary to get the exact distance between
     * two vectors (for example, when simply comparing the
     * distance between different vectors), it is advised to use
     * this method instead of distance.
     *
     * @param v The other vector.
     * 
     * @return The squared distance between this vector and v.
     * 
     * @see distance
     */
    inline T distanceSquared(const Vector2<T>& v) const;

    /**
     * Returns the dot product of this vector and the specified vector.
     *
     * @param v The vector to compute the dot product with.
     * 
     * @return The dot product.
     */
    inline T dot(const Vector2<T>& v) const;

    /**
     * Returns the dot product between the specified vectors.
     *
     * @param v1 The first vector.
     * @param v2 The second vector.
     * 
     * @return The dot product between the vectors.
     */
    static T dot(const Vector2<T>& v1, const Vector2<T>& v2);

    /**
     * Computes the length of this vector.
     *
     * @return The length of the vector.
     * 
     * @see lengthSquared
     */
    T length() const;

    /**
     * Returns the squared length of this vector.
     *
     * When it is not necessary to get the exact length of a
     * vector (for example, when simply comparing the lengths of
     * different vectors), it is advised to use this method
     * instead of length.
     *
     * @return The squared length of the vector.
     * 
     * @see length
     */
    inline T lengthSquared() const;

    /**
     * Negates this vector.
     */
    inline void negate();

    /**
     * Normalizes this vector.
     *
     * This method normalizes this Vec2 so that it is of
     * unit length (in other words, the length of the vector
     * after calling this method will be 1.0f). If the vector
     * already has unit length or if the length of the vector
     * is zero, this method does nothing.
     * 
     * @return This vector, after the normalization occurs.
     */
    void normalize();

    /**
     Get the normalized vector.
     */
    Vector2<T> getNormalized() const;

    /**
     * Scales all elements of this vector by the specified value.
     *
     * @param scalar The scalar value.
     */
    inline void scale(T scalar);

    /**
     * Scales each element of this vector by the matching component of scale.
     *
     * @param scale The vector to scale by.
     */
    inline void scale(const Vector2<T>& scale);

    /**
     * Rotates this vector by angle (specified in radians) around the given point.
     *
     * @param point The point to rotate around.
     * @param angle The angle to rotate by (in radians).
     */
    void rotate(const Vector2<T>& point, T angle);

    /**
     * Sets the elements of this vector to the specified values.
     *
     * @param xx The new x coordinate.
     * @param yy The new y coordinate.
     */
    inline void set(T xx, T yy);

    /**
     * Sets the elements of this vector from the values in the specified array.
     *
     * @param array An array containing the elements of the vector in the order x, y.
     */
    void set(const T* array);

    /**
     * Sets the elements of this vector to those in the specified vector.
     *
     * @param v The vector to copy.
     */
    inline void set(const Vector2<T>& v);

    /**
     * Sets this vector to the directional vector between the specified points.
     * 
     * @param p1 The first point.
     * @param p2 The second point.
     */
    inline void set(const Vector2<T>& p1, const Vector2<T>& p2);

    /**
    * Sets the elements of this vector to zero.
    */
    inline void setZero();

    /**
     * Subtracts this vector and the specified vector as (this - v)
     * and stores the result in this vector.
     *
     * @param v The vector to subtract.
     */
    inline void subtract(const Vector2<T>& v);

    /**
     * Subtracts the specified vectors and stores the result in dst.
     * The resulting vector is computed as (v1 - v2).
     *
     * @param v1 The first vector.
     * @param v2 The second vector.
     * @param dst The destination vector.
     */
    static void subtract(const Vector2<T>& v1, const Vector2<T>& v2, Vector2<T>* dst);

    /**
     * Updates this vector towards the given target using a smoothing function.
     * The given response time determines the amount of smoothing (lag). A longer
     * response time yields a smoother result and more lag. To force this vector to
     * follow the target closely, provide a response time that is very small relative
     * to the given elapsed time.
     *
     * @param target target value.
     * @param elapsedTime elapsed time between calls.
     * @param responseTime response time (in the same units as elapsedTime).
     */
    inline void smooth(const Vector2<T>& target, float elapsedTime, float responseTime);

    /**
     * Calculates the sum of this vector with the given vector.
     * 
     * Note: this does not modify this vector.
     * 
     * @param v The vector to add.
     * @return The vector sum.
     */
    inline Vector2<T> operator+(const Vector2<T>& v) const;

    /**
     * Adds the given vector to this vector.
     * 
     * @param v The vector to add.
     * @return This vector, after the addition occurs.
     */
    inline Vector2<T>& operator+=(const Vector2<T>& v);

    /**
     * Calculates the sum of this vector with the given vector.
     * 
     * Note: this does not modify this vector.
     * 
     * @param v The vector to add.
     * @return The vector sum.
     */
    inline Vector2<T> operator-(const Vector2<T>& v) const;

    /**
     * Subtracts the given vector from this vector.
     * 
     * @param v The vector to subtract.
     * @return This vector, after the subtraction occurs.
     */
    inline Vector2<T>& operator-=(const Vector2<T>& v);

    /**
     * Calculates the negation of this vector.
     * 
     * Note: this does not modify this vector.
     * 
     * @return The negation of this vector.
     */
    inline Vector2<T> operator-() const;

    /**
     * Calculates the scalar product of this vector with the given value.
     * 
     * Note: this does not modify this vector.
     * 
     * @param s The value to scale by.
     * @return The scaled vector.
     */
    inline Vector2<T> operator*(T s) const;

    /**
     * Scales this vector by the given value.
     * 
     * @param s The value to scale by.
     * @return This vector, after the scale occurs.
     */
    inline Vector2<T>& operator*=(T s);
    
    /**
     * Returns the components of this vector divided by the given constant
     *
     * Note: this does not modify this vector.
     *
     * @param s the constant to divide this vector with
     * @return a smaller vector
     */
    inline Vector2<T> operator/(T s) const;

    /**
     * Determines if this vector is less than the given vector.
     * 
     * @param v The vector to compare against.
     * 
     * @return True if this vector is less than the given vector, false otherwise.
     */
    inline bool operator<(const Vector2<T>& v) const;
    
    /**
     * Determines if this vector is greater than the given vector.
     *
     * @param v The vector to compare against.
     *
     * @return True if this vector is greater than the given vector, false otherwise.
     */
    inline bool operator>(const Vector2<T>& v) const;

    /**
     * Determines if this vector is equal to the given vector.
     * 
     * @param v The vector to compare against.
     * 
     * @return True if this vector is equal to the given vector, false otherwise.
     */
    inline bool operator==(const Vector2<T>& v) const;

    /**
     * Determines if this vector is not equal to the given vector.
     * 
     * @param v The vector to compare against.
     * 
     * @return True if this vector is not equal to the given vector, false otherwise.
     */
    inline bool operator!=(const Vector2<T>& v) const;

    //code added compatible for Point
public:
      /**
     * @js NA
     * @lua NA
     */
    inline void setPoint(T xx, T yy);
    /**
     * @js NA
     */
    bool equals(const Vector2<T>& target) const;
    
    /** @returns if points have fuzzy equality which means equal with some degree of variance.
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    bool fuzzyEquals(const Vector2<T>& target, T variance) const;

    /** Calculates distance between point an origin
     @return float
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline T getLength() const {
        return sqrtf(x*x + y*y);
    }

    /** Calculates the square length of a Vec2 (not calling sqrt() )
     @return float
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline T getLengthSq() const {
        return dot(*this); //x*x + y*y;
    }

    /** Calculates the square distance between two points (not calling sqrt() )
     @return float
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline T getDistanceSq(const Vector2<T>& other) const {
        return (*this - other).getLengthSq();
    }

    /** Calculates the distance between two points
     @return float
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline T getDistance(const Vector2<T>& other) const {
        return (*this - other).getLength();
    }

    /** @returns the angle in radians between this vector and the x axis
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline T getAngle() const {
        return atan2f(y, x);
    }

    /** @returns the angle in radians between two vector directions
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    T getAngle(const Vector2<T>& other) const;

    /** Calculates cross product of two points.
     @return float
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline T cross(const Vector2<T>& other) const {
        return x*other.y - y*other.x;
    }

    /** Calculates perpendicular of v, rotated 90 degrees counter-clockwise -- cross(v, perp(v)) >= 0
     @return Vec2
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline Vector2<T> getPerp() const {
        return Vector2<T>(-y, x);
    }
    
    /** Calculates midpoint between two points.
     @return Vec2
     @since v3.0
     * @js NA
     * @lua NA
     */
    inline Vector2<T> getMidpoint(const Vector2<T>& other) const
    {
        return Vector2<T>((x + other.x) / 2.0f, (y + other.y) / 2.0f);
    }
    
    /** Clamp a point between from and to.
     @since v3.0
     * @js NA
     * @lua NA
     */
    inline Vector2<T> getClampPoint(const Vector2<T>& min_inclusive, const Vector2<T>& max_inclusive) const
    {
        return Vector2<T>(clampf(x,min_inclusive.x,max_inclusive.x), clampf(y, min_inclusive.y, max_inclusive.y));
    }
    
    /** Run a math operation function on each point component
     * absf, floorf, ceilf, roundf
     * any function that has the signature: float func(float);
     * For example: let's try to take the floor of x,y
     * p.compOp(floorf);
     @since v3.0
     * @js NA
     * @lua NA
     */
    inline Vector2<T> compOp(std::function<T(T)> function) const
    {
        return Vector2<T>(function(x), function(y));
    }

    /** Calculates perpendicular of v, rotated 90 degrees clockwise -- cross(v, rperp(v)) <= 0
     @return Vec2
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline Vector2<T> getRPerp() const {
        return Vector2<T>(y, -x);
    }

    /** Calculates the projection of this over other.
     @return Vec2
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline Vector2<T> project(const Vector2<T>& other) const {
        return other * (dot(other)/other.dot(other));
    }

    /** Complex multiplication of two points ("rotates" two points).
     @return Vec2 vector with an angle of this.getAngle() + other.getAngle(),
     and a length of this.getLength() * other.getLength().
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline Vector2<T> rotate(const Vector2<T>& other) const {
        return Vector2<T>(x*other.x - y*other.y, x*other.y + y*other.x);
    }

    /** Unrotates two points.
     @return Vec2 vector with an angle of this.getAngle() - other.getAngle(),
     and a length of this.getLength() * other.getLength().
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline Vector2<T> unrotate(const Vector2<T>& other) const {
        return Vector2<T>(x*other.x + y*other.y, y*other.x - x*other.y);
    }

    /** Linear Interpolation between two points a and b
     @returns
        alpha == 0 ? a
        alpha == 1 ? b
        otherwise a value between a..b
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    inline Vector2<T> lerp(const Vector2<T>& other, float alpha) const {
        return *this * (1.f - alpha) + other * alpha;
    }

    /** Rotates a point counter clockwise by the angle around a pivot
     @param pivot is the pivot, naturally
     @param angle is the angle of rotation ccw in radians
     @returns the rotated point
     @since v2.1.4
     * @js NA
     * @lua NA
     */
    Vector2<T> rotateByAngle(const Vector2<T>& pivot, T angle) const;

    /**
     * @js NA
     * @lua NA
     */
    static inline Vector2<T> forAngle(const T a)
    {
        return Vector2<T>(cosf(a), sinf(a));
    }
    
    /** A general line-line intersection test
     @param A   the startpoint for the first line L1 = (A - B)
     @param B   the endpoint for the first line L1 = (A - B)
     @param C   the startpoint for the second line L2 = (C - D)
     @param D   the endpoint for the second line L2 = (C - D)
     @param S   the range for a hitpoint in L1 (p = A + S*(B - A))
     @param T   the range for a hitpoint in L2 (p = C + T*(D - C))
     @return    whether these two lines intersects.

     Note that to truly test intersection for segments we have to make
     sure that S & T lie within [0..1] and for rays, make sure S & T > 0
     the hit point is        C + T * (D - C);
     the hit point also is   A + S * (B - A);
     @since 3.0
     * @js NA
     * @lua NA
     */
    static bool isLineIntersect(const Vector2<T>& A, const Vector2<T>& B,
                                 const Vector2<T>& C, const Vector2<T>& D,
                                 T *s = nullptr, T *t = nullptr);
    
    /**
     returns true if Line A-B overlap with segment C-D
     @since v3.0
     * @js NA
     * @lua NA
     */
    static bool isLineOverlap(const Vector2<T>& A, const Vector2<T>& B,
                                const Vector2<T>& C, const Vector2<T>& D);
    
    /**
     returns true if Line A-B parallel with segment C-D
     @since v3.0
     * @js NA
     * @lua NA
     */
    static bool isLineParallel(const Vector2<T>& A, const Vector2<T>& B,
                   const Vector2<T>& C, const Vector2<T>& D);
    
    /**
     returns true if Segment A-B overlap with segment C-D
     @since v3.0
     * @js NA
     * @lua NA
     */
    static bool isSegmentOverlap(const Vector2<T>& A, const Vector2<T>& B,
                                 const Vector2<T>& C, const Vector2<T>& D,
                                 Vector2<T>* S = nullptr, Vector2<T>* E = nullptr);
    
    /**
     returns true if Segment A-B intersects with segment C-D
     @since v3.0
     * @js NA
     * @lua NA
     */
    static bool isSegmentIntersect(const Vector2<T>& A, const Vector2<T>& B, const Vector2<T>& C, const Vector2<T>& D);
    
    /**
     returns the intersection point of line A-B, C-D
     @since v3.0
     * @js NA
     * @lua NA
     */
    static Vector2<T> getIntersectPoint(const Vector2<T>& A, const Vector2<T>& B, const Vector2<T>& C, const Vector2<T>& D);
    
    /** equals to Vec2(0,0) */
    static const Vector2<T> ZERO;
    /** equals to Vec2(1,1) */
    static const Vector2<T> ONE;
    /** equals to Vec2(1,0) */
    static const Vector2<T> UNIT_X;
    /** equals to Vec2(0,1) */
    static const Vector2<T> UNIT_Y;
    /** equals to Vec2(0.5, 0.5) */
    static const Vector2<T> ANCHOR_MIDDLE;
    /** equals to Vec2(0, 0) */
    static const Vector2<T> ANCHOR_BOTTOM_LEFT;
    /** equals to Vec2(0, 1) */
    static const Vector2<T> ANCHOR_TOP_LEFT;
    /** equals to Vec2(1, 0) */
    static const Vector2<T> ANCHOR_BOTTOM_RIGHT;
    /** equals to Vec2(1, 1) */
    static const Vector2<T> ANCHOR_TOP_RIGHT;
    /** equals to Vec2(1, 0.5) */
    static const Vector2<T> ANCHOR_MIDDLE_RIGHT;
    /** equals to Vec2(0, 0.5) */
    static const Vector2<T> ANCHOR_MIDDLE_LEFT;
    /** equals to Vec2(0.5, 1) */
    static const Vector2<T> ANCHOR_MIDDLE_TOP;
    /** equals to Vec2(0.5, 0) */
    static const Vector2<T> ANCHOR_MIDDLE_BOTTOM;
};

/**
 * Calculates the scalar product of the given vector with the given value.
 * 
 * @param x The value to scale by.
 * @param v The vector to scale.
 * @return The scaled vector.
 */
template <class T> inline Vector2<T> operator*(T x, const Vector2<T>& v);

typedef Vector2<float> Vec2;
typedef Vector2<float> Point;

NS_CC_MATH_END

/**
 end of base group
 @}
 */

#include "math/Vec2.inl"

#endif // MATH_VEC2_H

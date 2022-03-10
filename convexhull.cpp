/**
 */
#include "convexhull.h"

#ifdef _DEBUG
#    include <stdio.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <random>

#include "halfedge.h"

namespace convexhull
{
namespace
{
    static constexpr u32 Invalid = static_cast<u32>(-1);
    static constexpr f32 ConvexEpsilon = 1.0e-6f;
    static constexpr f32 MinusConvexEpsilon = -1.0e-6f;
    static constexpr u32 MaxBlockSize = 256;

    f32 absolute(f32 x)
    {
        return fabsf(x);
    }

    f32 maximum(f32 x0, f32 x1)
    {
        return x0 < x1 ? x1 : x0;
    }

    f32 minimum(f32 x0, f32 x1)
    {
        return x0 < x1 ? x0 : x1;
    }

    bool sameSign(f32 x0, f32 x1)
    {
        return (zero < x0) == (zero < x1);
    }

    bool isEqual(f32 x0, f32 x1)
    {
        return std::abs(x1 - x0) < ConvexEpsilon;
    }

} // namespace

//---------------------------------------------
Vector2 Vector2::operator-() const
{
    return {-x_, -y_};
}

f32 Vector2::lengthSqr() const
{
    return x_ * x_ + y_ * y_;
}

Vector2 operator+(const Vector2& x0, const Vector2& x1)
{
    return {x0.x_ + x1.x_, x0.y_ + x1.y_};
}

Vector2 operator-(const Vector2& x0, const Vector2& x1)
{
    return {x0.x_ - x1.x_, x0.y_ - x1.y_};
}

Vector2 perpendicular(const Vector2& x)
{
    return {x.y_, -x.x_};
}

f32 dot(const Vector2& x0, const Vector2& x1)
{
    return x0.x_ * x1.x_ + x0.y_ * x1.y_;
}

Vector2 tripleProduct(const Vector2& x0, const Vector2& x1, const Vector2& x2)
{
    f32 d02 = dot(x0, x2);
    f32 d12 = dot(x1, x2);
    return {x1.x_ * d02 - x0.x_ * d12, x1.y_ * d02 - x0.y_ * d12};
}

//---------------------------------------------
Vector3 Vector3::operator-() const
{
    return {-x_, -y_, -z_};
}

f32 Vector3::lengthSqr() const
{
    return x_ * x_ + y_ * y_ + z_ * z_;
}

Vector3& Vector3::operator+=(const Vector3& x)
{
    x_ += x.x_;
    y_ += x.y_;
    z_ += x.z_;
    return *this;
}

Vector3 operator+(const Vector3& x0, const Vector3& x1)
{
    return {x0.x_ + x1.x_, x0.y_ + x1.y_, x0.z_ + x1.z_};
}

Vector3 operator-(const Vector3& x0, const Vector3& x1)
{
    return {x0.x_ - x1.x_, x0.y_ - x1.y_, x0.z_ - x1.z_};
}

Vector3 operator*(f32 x0, const Vector3& x1)
{
    return {x0 * x1.x_, x0 * x1.y_, x0 * x1.z_};
}

Vector3 operator*(const Vector3& x0, f32 x1)
{
    return {x0.x_ * x1, x0.y_ * x1, x0.z_ * x1};
}

f32 dot(const Vector3& x0, const Vector3& x1)
{
    return x0.x_ * x1.x_ + x0.y_ * x1.y_ + x0.z_ * x1.z_;
}

Vector3 cross(const Vector3& x0, const Vector3& x1)
{
    f32 x = x0.y_ * x1.z_ - x0.z_ * x1.y_;
    f32 y = x0.z_ * x1.x_ - x0.x_ * x1.z_;
    f32 z = x0.x_ * x1.y_ - x0.y_ * x1.x_;
    return {x, y, z};
}

Vector3 normalize(const Vector3& x)
{
    f32 inv = 1.0f / ::sqrtf(x.lengthSqr());
    return {x.x_ * inv, x.y_ * inv, x.z_ * inv};
}

Vector3 normalizeSafe(const Vector3& x)
{
    f32 l = ::sqrtf(x.lengthSqr());
    if(l < 1.0e-6f) {
        return {0.0f, 0.0f, 0.0f};
    }
    f32 inv = 1.0f / l;
    return {x.x_ * inv, x.y_ * inv, x.z_ * inv};
}

bool isEqual(const Vector3& x0, const Vector3& x1)
{
    return isEqual(x0.x_, x1.x_) && isEqual(x0.y_, x1.y_) && isEqual(x0.z_, x1.z_);
}

void orthonormalBasis(Vector3& binormal0, Vector3& binormal1, const Vector3& normal)
{
    if(normal.z_ < -0.9999999f) {
        binormal0 = {0.0f, -1.0f, 0.0f};
        binormal1 = {-1.0f, 0.0f, 0.0f};
        return;
    }

    const f32 a = 1.0f / (1.0f + normal.z_);
    const f32 b = -normal.x_ * normal.y_ * a;
    binormal0 = {1.0f - normal.x_ * normal.x_ * a, b, -normal.x_};
    binormal1 = {b, 1.0f - normal.y_ * normal.y_ * a, -normal.y_};
}

namespace
{
    //--------------------------------------------------------------
    struct Flags
    {
        Flags();
        ~Flags();
        void clear(u32 capacity);
        void set(u32 index);
        bool check(u32 index) const;

        u32 capacity_;
        u8* flags_;
    };

    Flags::Flags()
        : capacity_(0)
        , flags_(nullptr)
    {
    }

    Flags::~Flags()
    {
        ::free(flags_);
    }

    void Flags::clear(u32 capacity)
    {
        capacity = (capacity + 7U) >> 3U;
        if(capacity_ < capacity) {
            while(capacity_ < capacity) {
                capacity_ += 128U;
            }
            ::free(flags_);
            flags_ = reinterpret_cast<u8*>(::malloc(sizeof(u8) * capacity_));
        }
        ::memset(flags_, 0, sizeof(u8) * capacity_);
    }

    void Flags::set(u32 index)
    {
        u32 i = index >> 3U;
        u8 bit = static_cast<u8>(0x01U << (index & 0x07U));
        flags_[i] |= bit;
    }

    bool Flags::check(u32 index) const
    {
        u32 i = index >> 3U;
        u8 bit = static_cast<u8>(0x01U << (index & 0x07U));
        return bit == (flags_[i] & bit);
    }

    //--------------------------------------------------------------
    struct Plane
    {
        Vector3 normal_;
        f32 d_;
    };

    //--------------------------------------------------------------
    f32 angle(const Vector2& p0, const Vector2& p1, const Vector2& p2)
    {
        f32 dx0 = (p1.x_ - p0.x_);
        f32 dy0 = (p1.y_ - p0.y_);
        f32 dx1 = (p2.x_ - p1.x_);
        f32 dy1 = (p2.y_ - p1.y_);
        return dy0 * dx1 - dx0 * dy1;
    }

    bool distance(const Vector2& p0, const Vector2& p1, const Vector2& p2)
    {
        Vector2 d0 = p1 - p0;
        Vector2 d1 = p2 - p0;
        return (dot(d0, d0) < dot(d1, d1)) ? false : true;
    }

    bool less(const Vector2& p0, const Vector2& p1, const Vector2& p2)
    {
        f32 dx0 = (p1.x_ - p0.x_);
        f32 dy0 = (p1.y_ - p0.y_);
        f32 dx1 = (p2.x_ - p0.x_);
        f32 dy1 = (p2.y_ - p0.y_);
        return dy0 * dx1 < dy1 * dx0;
    }

    struct Predicate2
    {
        bool operator()(const Vector2& p1, const Vector2& p2) const noexcept
        {
            return less(p0_, p1, p2);
        }
        Vector2 p0_;
    };

    u32 orientation(const Vector2& p0, const Vector2& p1, const Vector2& p2)
    {
        f32 d = angle(p0, p1, p2);
        if(d < MinusConvexEpsilon) {
            return 2;
        }
        if(ConvexEpsilon < d) {
            return 1;
        }
        return 0;
    }

    f32 angle(const Vector3& p0, const Vector3& p1, const Vector3& p2)
    {
        Vector3 d0 = p1 - p0;
        Vector3 d1 = p2 - p1;
        return dot(d0, d1);
    }

    bool less(const Vector3& p0, const Vector3& p1, const Vector3& p2)
    {
        Vector3 d0 = p1 - p0;
        Vector3 d1 = p2 - p1;
        f32 d = dot(d0, d1);
        return 0.0f <= d;
    }

    u32 GrahamScan(Vector2* result, u32 begin, u32 end, Vector2* source)
    {
        u32 size = end - begin;
        assert(3 <= size);
        source += begin;
        u32 i0 = 0;
        for(u32 i = 1; i < size; ++i) {
            const Vector2& p0 = source[i0];
            const Vector2& p1 = source[i];
            if((p1.y_ < p0.y_) || ((p1.y_ - p0.y_) < ConvexEpsilon && p1.x_ < p0.x_)) {
                i0 = i;
            }
        }
        std::swap(source[0], source[i0]);
        // Sort in counter-clockwise
        Predicate2 predicate;
        predicate.p0_ = source[0];
        std::sort(source + 1, source + size, predicate);
#if 0
            u32 m = 1;
            for (u32 i = 1; i < size; ++i) {
                const Vector2& p0 = source[begin];
                while (i < (size - 1) && orient(p0, source[i], source[i + 1]) <= 0) {
                    ++i;
                }
                source[m] = source[i];
                ++m;
            }
            size = m;
#endif

        u32 count = 3;
        result[0] = source[0];
        result[1] = source[1];
        result[2] = source[2];
        for(u32 i = 3; i < size; ++i) {
            const Vector2& p0 = source[i];
            while(1 < count && orientation(result[count - 2], result[count - 1], p0) < 2) {
                --count;
            }
            result[count] = p0;
            ++count;
        }
        return count;
    }

    std::optional<u32> JarvisMarch(Vector2* result, u32 size, const Vector2* source, u32 maxCount = 0xFFFFFFFFU)
    {
        u32 i0 = 0;
        for(u32 i = 1; i < size; ++i) {
            const Vector2& p0 = source[i0];
            const Vector2& p1 = source[i];
            if(p1.x_ < p0.x_) {
                i0 = i;
            }
        }

        u32 count = 0;
        u32 p = i0;
        do {
            if(maxCount < count) {
                return std::nullopt;
            }
            result[count] = source[p];
            ++count;
            u32 q = (p + 1) % size;
            for(u32 i = 0; i < size; ++i) {
                if(2 == orientation(source[p], source[i], source[q])) {
                    q = i;
                }
            }
            p = q;
        } while(p != i0);
        return count;
    }

    f32 distanceSqr(const Vector3& origin, const Vector3& direction, const Vector3& point)
    {
        Vector3 diff = origin - point;
        f32 d = dot(diff, direction);
        Vector3 p = diff - d * direction;
        return p.lengthSqr();
    }

    f32 distance(f32 distance, const Vector3& normal, const Vector3& point)
    {
        return dot(normal, point) + distance;
    }

    f32 squaredArea(const Vector3& p0, const Vector3& p1, const Vector3& p2)
    {
        Vector3 d0 = p1 - p0;
        Vector3 d1 = p2 - p0;
        return cross(d0, d1).lengthSqr();
    }

    f32 signedVolume(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3)
    {
        Vector3 d0 = p1 - p0;
        Vector3 d1 = p2 - p0;
        Vector3 d2 = p3 - p0;
        return dot(cross(d0, d1), d2);
    }

    //---------------------------------------------
    class Quickhull
    {
    public:
        /**
         * @brief 
         * @param indices 
         * @param size 
         * @param source 
         * @return convexhull indices 
        */
        Result build(Array<u32>& indices, u32 size, const Vector3* source);

    private:
        static constexpr f32 StaticTolerance = ConvexEpsilon;

        struct Edge
        {
            u32 p0_;
            u32 p1_;

            bool operator==(const Edge& other) const
            {
                return p0_ == other.p0_ && p1_ == other.p1_;
            }

            bool operator<(const Edge& other) const
            {
                return p0_<other.p0_;
            }
        };

        /**
         * @brief Find max extent
         * @return min and max point's indices
        */
        std::tuple<u32, u32> findMinMax();

        /**
         * @brief Find a simplex
         * @param indices [in] ...
         * @param simplex [out] ... 
         * @return 
        */
        Result findInitialSimplex(Array<u32>& indices, u32 simplex[4]);
        bool checkSimplex(const u32 simplex[4]) const;

        u32 getFurthestPointOnFace(u32 p0, u32 p1, u32 p2);
        void findVisibleFaces(Array<Edge>& visibles, Flags& visibility, u32 point);
        void findContour(Array<Edge>& edges, const Array<Edge>& visibles, const Flags& visibility) const;
        void removeFaces(Array<Edge>& faces, Array<Edge>& visibles);
        void addFaces(Array<Edge>& faces, const Array<Edge>& edges, u32 point);
        bool isConvex(const u32 triangle0[3], const u32 triangle1[3]) const;
        bool isConvex(const Edge& e0, const Edge& e1, u32 point) const;
        bool mergeNonConvexFaces(Array<Edge>& faces, u32 point);
        Plane facePlane(u32 p0, u32 p1, u32 p2) const;
        bool isOnFace(const Plane& plane, const Vector3& point) const;
        bool checkFullContour(const Array<Edge>& edges) const;
        void checkInsideSimplex(u32 p0, u32 p1, u32 p2, u32 p3);

        static Vector2 project(const Vector3& p, const Plane& b0, const Plane& b1);
        Result build2d(Array<u32>& indices, const Plane& b0, const Plane& b1);

        u32 size_;
        const Vector3* source_;

        Vector3 normal_;
        Vector3 binormal0_;
        Vector3 binormal1_;
        f32 tolerance_;
        Flags visited_;
        halfedge::HalfedgeMesh mesh_;
    };

    Result Quickhull::build(Array<u32>& indices, u32 size, const Vector3* source)
    {
        size_ = size;
        source_ = source;
        tolerance_ = StaticTolerance;

        u32 simplex[4];
        Result result = findInitialSimplex(indices, simplex);
        if(Result::Convexhull != result) {
            return result;
        }
        assert(checkSimplex(simplex));
        visited_.clear(size_);
        for(u32 i = 0; i < 4; ++i) {
            visited_.set(simplex[i]);
        }
        checkInsideSimplex(simplex[0], simplex[1], simplex[2], simplex[3]);

        {
            mesh_.addTriangle(simplex[0], simplex[1], simplex[2]);
            mesh_.addTriangle(simplex[1], simplex[0], simplex[3]);
            mesh_.addTriangle(simplex[2], simplex[1], simplex[3]);
            mesh_.addTriangle(simplex[0], simplex[2], simplex[3]);
            assert(mesh_.validate());
        }

        Array<Edge> faces; //Should be sorted
        faces.push({simplex[0], simplex[1]});
        faces.push({simplex[0], simplex[2]});
        faces.push({simplex[1], simplex[0]});
        faces.push({simplex[2], simplex[1]});

        Array<Edge> edges;
        Array<Edge> visibles;
        Flags visibility;
        while(0 < faces.size()) {
            const Edge face = faces.pop();
            if(Invalid == mesh_.findHalfedge(face.p0_, face.p1_)){
                continue;
            }
            u32 facePoints[3];
            mesh_.getFacePoints(facePoints, face.p0_, face.p1_);
            u32 furthest = getFurthestPointOnFace(facePoints[0], facePoints[1], facePoints[2]);
            if(Invalid == furthest) {
                continue;
            }
            visited_.set(furthest);
            findVisibleFaces(visibles, visibility, furthest);
            findContour(edges, visibles, visibility);
            if(!mergeNonConvexFaces(edges, furthest)) {
                continue;
            }
            if(edges.size() < 3) {
                continue;
            }
            if(!checkFullContour(edges)){
                continue;
            }
            assert(0 == (visibles.size()%3));
            removeFaces(faces, visibles);
            assert(mesh_.validate(false));
            addFaces(faces, edges, furthest);
            checkInsideSimplex(facePoints[0], facePoints[2], facePoints[1], furthest);
            assert(mesh_.validate());
        }

        indices.clear();
        indices.reserve(3 * mesh_.getNumHalfedges());
        for(u32 i = mesh_.beginHalfedge(); i != mesh_.endHalfedge(); i = mesh_.nextHalfedge(i)) {
            u32 points[3];
            mesh_.getFacePoints(points, mesh_.getHalfedge(i).edge_.p0_, mesh_.getHalfedge(i).edge_.p1_);
            indices.push(points[0]);
            indices.push(points[1]);
            indices.push(points[2]);
        }
        return Result::Convexhull;
    }

    std::tuple<u32, u32> Quickhull::findMinMax()
    {
        u32 minX[3];
        u32 maxX[3];
        f32 minV[3];
        f32 maxV[3];
        for(u32 i = 0; i < 3; ++i) {
            minX[i] = 0;
            maxX[i] = 0;
            minV[i] = source_[0][i];
            maxV[i] = source_[0][i];
            for(u32 j = 1; j < size_; ++j) {
                if(source_[j][i] < minV[i]) {
                    minX[i] = j;
                    minV[i] = source_[j][i];
                }
                if(maxV[i] < source_[j][i]) {
                    maxX[i] = j;
                    maxV[i] = source_[j][i];
                }
            }
        }
        f32 maxD = maxV[0] - minV[0];
        u32 result = 0;
        for(u32 i = 1; i < 3; ++i) {
            f32 d = maxV[i] - minV[i];
            if(maxD < d) {
                maxD = d;
                result = i;
            }
        }

        f32 t = maximum(absolute(minV[0]), absolute(maxV[0]))
                + maximum(absolute(minV[1]), absolute(maxV[1]))
                + maximum(absolute(minV[2]), absolute(maxV[2]));
        tolerance_ = std::numeric_limits<f32>::epsilon() * 3.0f * t;

        return tolerance_ < maxD ? std::make_tuple(minX[result], maxX[result]) : std::make_tuple(0U, 0U);
    }

    Result Quickhull::findInitialSimplex(Array<u32>& indices, u32 simplex[4])
    {
        auto [p0, p1] = findMinMax();
        if(p0 == p1) {
            indices.push(p0);
            indices.push(p0);
            indices.push(p0);
            return Result::Point;
        }
        Vector3 origin = source_[p0];
        Vector3 direction = normalize(source_[p1] - source_[p0]);

        f32 maxDistance = distanceSqr(origin, direction, source_[0]);
        u32 p2 = 0;
        for(u32 i = 1; i < size_; ++i) {
            f32 d = distanceSqr(origin, direction, source_[i]);
            if(maxDistance < d) {
                maxDistance = d;
                p2 = i;
            }
        }
        if(maxDistance < (tolerance_ * 10.0f)) {
            indices.push(p0);
            indices.push(p0);
            indices.push(p1);
            return Result::Colinear;
        }
        Vector3 d0 = source_[p1] - source_[p0];
        Vector3 d1 = source_[p2] - source_[p0];
        Vector3 normal = normalizeSafe(cross(d0, d1));
        f32 pd = dot(normal, source_[p0]);

        maxDistance = 0.0f;
        u32 p3 = 0;
        for(u32 i = 0; i < size_; ++i) {
            f32 d = absolute(distance(pd, normal, source_[i]));
            if(maxDistance < d) {
                maxDistance = d;
                p3 = i;
            }
        }
        if(maxDistance < (tolerance_ * 10.0f)) {
            Plane b0, b1;
            orthonormalBasis(b0.normal_, b1.normal_, normal);
            b0.d_ = -dot(b0.normal_, source_[p0]);
            b1.d_ = -dot(b1.normal_, source_[p0]);
            return build2d(indices, b0, b1);
        }

        f32 d = dot(source_[p3], normal) - pd;
        if(d < 0.0f) {
            simplex[0] = p0;
            simplex[1] = p1;
            simplex[2] = p2;
            simplex[3] = p3;
        } else {
            simplex[0] = p0;
            simplex[1] = p2;
            simplex[2] = p1;
            simplex[3] = p3;
        }
        return Result::Convexhull;
    }

    bool Quickhull::checkSimplex(const u32 simplex[4]) const
    {
        u32 triangle0[3] = {simplex[0], simplex[1], simplex[2]};
        u32 triangle1[3] = {simplex[1], simplex[0], simplex[3]};
        u32 triangle2[3] = {simplex[2], simplex[1], simplex[3]};
        u32 triangle3[3] = {simplex[0], simplex[2], simplex[3]};

        if(!isConvex(triangle0, triangle1)) {
            return false;
        }

        if(!isConvex(triangle0, triangle2)) {
            return false;
        }

        if(!isConvex(triangle0, triangle3)) {
            return false;
        }

        if(!isConvex(triangle1, triangle2)) {
            return false;
        }

        if(!isConvex(triangle1, triangle3)) {
            return false;
        }

        if(!isConvex(triangle2, triangle3)) {
            return false;
        }
        return true;
    }

    u32 Quickhull::getFurthestPointOnFace(u32 i0, u32 i1, u32 i2)
    {
        const Vector3& p0 = source_[i0];
        const Vector3& p1 = source_[i1];
        const Vector3& p2 = source_[i2];
        normal_ = normalizeSafe(cross(p1 - p0, p2 - p0));
        orthonormalBasis(binormal0_, binormal1_, normal_);
        f32 d = dot(normal_, p0);
        f32 distance = tolerance_;
        u32 point = Invalid;
        for(u32 i = 0; i < size_; ++i) {
            if(visited_.check(i)) {
                continue;
            }
            f32 t = dot(normal_, source_[i]) - d;
            if(distance < t) {
                point = i;
                distance = t;
            }
        }
        return point;
    }

    void Quickhull::findVisibleFaces(Array<Edge>& visibles, Flags& visibility, u32 p)
    {
        const Vector3& point = source_[p];
        visibles.clear();
        visibility.clear(mesh_.getCapacityHalfedges());
        u32 facePoints[3];
        for(u32 h = mesh_.beginHalfedge(); h != mesh_.endHalfedge(); h = mesh_.nextHalfedge(h)) {
            if(visibility.check(h)) {
                continue;
            }
            const halfedge::Halfedge& halfedge = mesh_.getHalfedge(h);
            mesh_.getFacePoints(facePoints, halfedge.edge_.p0_, halfedge.edge_.p1_);
            const Vector3& p0 = source_[facePoints[0]];
            const Vector3& p1 = source_[facePoints[1]];
            const Vector3& p2 = source_[facePoints[2]];
            Vector3 n = normalizeSafe(cross(p1 - p0, p2 - p0));
            f32 distance = dot(n, point) - dot(n, p0);
            if(distance <= tolerance_) {
                continue;
            }
            visibles.push({halfedge.edge_.p0_, halfedge.edge_.p1_});
            visibility.set(h);
            u32 next = mesh_.findHalfedge(halfedge.next_.p0_, halfedge.next_.p1_);
            do {
                const halfedge::Halfedge& tmp = mesh_.getHalfedge(next);
                visibles.push({tmp.edge_.p0_, tmp.edge_.p1_});
                visibility.set(next);
                next = mesh_.findHalfedge(tmp.next_.p0_, tmp.next_.p1_);
            } while(h != next);
        }
    }

    void Quickhull::findContour(Array<Edge>& edges, const Array<Edge>& visibles, const Flags& visibility) const
    {
        edges.clear();
#if 1
        for(u32 i = 0; i < visibles.size(); ++i) {
            u32 h = mesh_.findHalfedge(visibles[i].p0_, visibles[i].p1_);
            assert(Invalid != h);
            assert(visibility.check(h));
            const halfedge::Halfedge& halfedge = mesh_.getHalfedge(h);
            assert(halfedge.edge_.p0_ != halfedge.edge_.p1_);

            u32 opposite = mesh_.findHalfedge(halfedge.edge_.p1_, halfedge.edge_.p0_);
            assert(Invalid != opposite);
            if(!visibility.check(opposite)) {
                edges.push({halfedge.edge_.p0_, halfedge.edge_.p1_});
            }
        }
#endif
    }

    void Quickhull::removeFaces(Array<Edge>& faces, Array<Edge>& visibles)
    {
        u32 points[3];
        for(u32 i = 0; i < visibles.size(); ++i) {
            u32 p0 = visibles[i].p0_;
            u32 p1 = visibles[i].p1_;

            u32 index = faces.binarySearch({p0, p1});
            if(Array<Edge>::Invalid != index){
                faces.removeAt(index);
            }

            u32 h = mesh_.findHalfedge(p0, p1);
            if(Invalid == h) {
                continue;
            }
            mesh_.getFacePoints(points, p0, p1);
            mesh_.removeTriangle(points);
        }
    }

    void Quickhull::addFaces(Array<Edge>& faces, const Array<Edge>& edges, u32 point)
    {
        for(u32 i = 0; i < edges.size(); ++i) {
            Edge edge = edges[i];
            mesh_.addTriangle(edge.p0_, edge.p1_, point);
            u32 index = faces.binarySearch(edge);
            if(Array<Edge>::Invalid == index){
                faces.insert(edge);
            }
        }
    }

    bool Quickhull::isConvex(const u32 triangle0[3], const u32 triangle1[3]) const
    {
        Vector3 d0 = source_[triangle0[1]] - source_[triangle0[0]];
        Vector3 d1 = source_[triangle0[2]] - source_[triangle0[0]];
        Vector3 n0 = normalizeSafe(cross(d0, d1));
        Vector3 c0 = (source_[triangle0[0]] + source_[triangle0[1]] + source_[triangle0[2]]) * (1.0f / 3.0f);

        Vector3 d2 = source_[triangle1[1]] - source_[triangle1[0]];
        Vector3 d3 = source_[triangle1[2]] - source_[triangle1[0]];
        Vector3 n1 = normalizeSafe(cross(d2, d3));
        Vector3 c1 = (source_[triangle1[0]] + source_[triangle1[1]] + source_[triangle1[2]]) * (1.0f / 3.0f);

        f32 t0 = dot(n0, c1) - dot(n0, c0);
        if(tolerance_ < t0) {
            return false;
        }
        f32 t1 = dot(n1, c0) - dot(n1, c1);
        if(tolerance_ < t1) {
            return false;
        }
        return true;
    }

    bool Quickhull::isConvex(const Edge& e0, const Edge& e1, u32 point) const
    {
        Vector3 d0 = source_[e0.p1_] - source_[e0.p0_];
        Vector3 d1 = source_[point] - source_[e0.p0_];
        Vector3 n0 = normalizeSafe(cross(d0, d1));
        Vector3 c0 = (source_[e0.p0_] + source_[e0.p1_] + source_[point]) * (1.0f / 3.0f);

        Vector3 d2 = source_[e1.p1_] - source_[e1.p0_];
        Vector3 d3 = source_[point] - source_[e1.p0_];
        Vector3 n1 = normalizeSafe(cross(d2, d3));
        Vector3 c1 = (source_[e1.p0_] + source_[e1.p1_] + source_[point]) * (1.0f / 3.0f);

        f32 t0 = dot(n0, c1) - dot(n0, c0);
        if(tolerance_ < t0) {
            return false;
        }
        f32 t1 = dot(n1, c0) - dot(n1, c1);
        if(tolerance_ < t1) {
            return false;
        }
        return true;
    }

    bool Quickhull::mergeNonConvexFaces(Array<Edge>& faces, u32 point)
    {
        for(u32 i = 0; i < faces.size();) {
            u32 next = Invalid;
            for(u32 j = i + 1; j < faces.size(); ++j) {
                if(faces[i].p1_ == faces[j].p0_) {
                    next = j;
                    break;
                }
            }
            if(Invalid == next) {
                ++i;
                continue;
            }
            if(isConvex(faces[i], faces[next], point)) {
                ++i;
            } else {
                faces[i].p1_ = faces[next].p1_;
                faces.removeAt(next);
            }
        }
        return true;
    }

    Plane Quickhull::facePlane(u32 p0, u32 p1, u32 p2) const
    {
        Vector3 normal = normalizeSafe(cross(source_[p1] - source_[p0], source_[p2] - source_[p0]));
        return {normal, -dot(normal, source_[p0])};
    }

    bool Quickhull::isOnFace(const Plane& plane, const Vector3& point) const
    {
        return tolerance_ < (dot(plane.normal_, point) + plane.d_);
    }

    bool Quickhull::checkFullContour(const Array<Edge>& edges) const
    {
        for(u32 i=0; i<edges.size(); ++i){
            u32 next = Invalid;
            for(u32 j=0; j<edges.size(); ++j){
                if(edges[i].p1_ == edges[j].p0_){
                    assert(i!=j);
                    next = j;
                    break;
                }
            }
            if(Invalid == next){
                return false;
            }
            for(u32 j=i+1; j<edges.size(); ++j){
                if(edges[i].p1_ == edges[j].p1_){
                    return false;
                }
            }
        }
        return true;
    }

    void Quickhull::checkInsideSimplex(u32 p0, u32 p1, u32 p2, u32 p3)
    {
        Plane plane0 = facePlane(p0, p1, p2);
        Plane plane1 = facePlane(p1, p0, p3);
        Plane plane2 = facePlane(p2, p1, p3);
        Plane plane3 = facePlane(p0, p2, p3);
        for(u32 i = 0; i < size_; ++i) {
            if(visited_.check(i)) {
                continue;
            }
            if(isOnFace(plane0, source_[i])) {
                continue;
            }
            if(isOnFace(plane1, source_[i])) {
                continue;
            }
            if(isOnFace(plane2, source_[i])) {
                continue;
            }
            if(isOnFace(plane3, source_[i])) {
                continue;
            }
            visited_.set(i);
        }
    }

    Vector2 Quickhull::project(const Vector3& p, const Plane& b0, const Plane& b1)
    {
        return {dot(p, b0.normal_)+b0.d_, dot(p, b1.normal_)+b1.d_};
    }

    Result Quickhull::build2d(Array<u32>& indices, const Plane& b0, const Plane& b1)
    {
        u32 i0 = 0;
        for(u32 i = 1; i < size_; ++i) {
            Vector2 p0 = project(source_[i0], b0, b1);
            Vector2 p1 = project(source_[i], b0, b1);
            if(p1.x_ < p0.x_) {
                i0 = i;
            }
        }

        Array<u32> points;
        points.reserve(size_);
        u32 count = 0;
        u32 p = i0;
        do {
            points.push(p);
            ++count;
            u32 q = (p + 1) % size_;

            Vector2 p0 = project(source_[p], b0, b1);

            for(u32 i = 0; i < size_; ++i) {
                Vector2 p1 = project(source_[i], b0, b1);
                Vector2 p2 = project(source_[q], b0, b1);
                if(2 == orientation(p0, p1, p2)) {
                    q = i;
                }
            }
            p = q;
        } while(p != i0);
        switch(points.size()){
        case 0:
            return Result::Error;
        case 1:
            indices.push(points[0]);
            indices.push(points[0]);
            indices.push(points[0]);
            return Result::Point;
        case 2:
            indices.push(points[0]);
            indices.push(points[0]);
            indices.push(points[1]);
            return Result::Colinear;
        }
        for(u32 i=2; i<points.size(); ++i){
            indices.push(points[0]);
            indices.push(points[i-1]);
            indices.push(points[i]);
        }
        return Result::Coplanar;
    }
} // namespace

void convexhull2(Array<Vector2>& result, u32 size, Vector2* points)
{
    result.clear();
    if(size <= 3) {
        for(u32 i = 0; i < size; ++i) {
            result.push(points[i]);
        }
        return;
    }
    result.resize(size);
    Array<Vector2> buffer;
    buffer.resize(size);
    u32 m;
    u32 t = 2;
    do {
        m = 1U << (1U << t);
        ++t;
        u32 count = 0;
        for(u32 i = 0; i < size; i += m) {
            u32 end = i + m;
            end = end <= size ? end : size;
            count += GrahamScan(&buffer[count], i, end, points);
            assert(count <= size);
        }
        std::optional<u32> num = JarvisMarch(&result[0], count, &buffer[0], m);
        if(num.has_value()) {
            result.resize(num.value());
            return;
        }
    } while(m < MaxBlockSize);
    u32 count = GrahamScan(&result[0], 0, size, points);
    result.resize(count);
}

Result convexhull3(Array<u32>& indices, u32 size, const Vector3* points)
{
    indices.clear();
    Quickhull quickhull;
    return quickhull.build(indices, size, points);
}

bool validate(const Array<Vector2>& points, const char* name)
{
    if(points.size() <= 3) {
        return true;
    }
    FILE* file = fopen(name, "wb");
    if(nullptr == file){
        return true;
    }
    for(u32 i = 0; i < points.size(); ++i) {
        fprintf(file, "%f, %f\n", points[i].x_, points[i].y_);
    }
    fclose(file);
    for(u32 i = 2; i < points.size(); ++i) {
        const Vector2& p0 = points[i - 2];
        const Vector2& p1 = points[i - 1];
        const Vector2& p2 = points[i];
        u32 d = orientation(p0, p1, p2);
        if(d == 1) {
            return false;
        }
    }
    return true;
}

Validation validate(const Array<u32>& indices, u32 size, const Vector3* points, const char* name)
{
#if 1
    u32 triangles = indices.size() / 3;
    FILE* file = fopen(name, "wb");
    if(nullptr == file){
        return {0, 0};
    }
    fprintf(file, "ply\nformat ascii 1.0\n");
    fprintf(file, "element vertex %d\n", size);
    fprintf(file, "property float x\nproperty float y\nproperty float z\n");
    fprintf(file, "element face %d\n", triangles);
    fprintf(file, "property list uchar int vertex_index\n");
    fprintf(file, "end_header\n");
    for(u32 i = 0; i < size; ++i) {
        fprintf(file, "%f %f %f\n", points[i].x_, points[i].y_, points[i].z_);
    }
    for(u32 i = 0; i < triangles; ++i) {
        u32 n = i * 3;
        fprintf(file, "3 %d %d %d\n", indices[n + 0], indices[n + 1], indices[n + 2]);
    }
    fclose(file);

    Flags flags;
    flags.clear(size);
    f32 maxDistance = -std::numeric_limits<f32>::infinity();
    for(u32 i = 0; i < triangles; ++i) {
        u32 n = i * 3;
        u32 i0 = indices[n + 0];
        u32 i1 = indices[n + 1];
        u32 i2 = indices[n + 2];
        Vector3 normal = normalizeSafe(cross(points[i1]-points[i0], points[i2]-points[i0]));
        f32 d = dot(normal, points[i0]);
        for(u32 j=0; j<size; ++j){
            f32 distance = dot(points[j], normal) - d;
            if(maxDistance<distance){
                maxDistance = distance;
            }
            if(0.0f<distance){
                flags.set(j);
            }
        }
    }
    u32 count = 0;
    for(u32 i=0; i<size; ++i){
        if(flags.check(i)){
            ++count;
        }
    }
    return {maxDistance, count};
#endif
}
} // namespace convexhull

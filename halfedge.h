#ifndef INC_HALFEDGE_H_
#define INC_HALFEDGE_H_
/**
 */
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <utility>

namespace halfedge
{
using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

using f32 = float;
using f64 = double;

static constexpr u32 Invalid = static_cast<u32>(-1);

//--------------------------------------------------------------
u32 wyhash(u32 x0, u32 x1);

//--------------------------------------------------------------
struct Edge
{
    u32 p0_;
    u32 p1_;

    bool operator==(const Edge& other) const
    {
        return p0_ == other.p0_ && p1_ == other.p1_;
    }
};

struct Halfedge
{
    Edge edge_;
    Edge next_;

    static u32 calcHash(const Halfedge& edge)
    {
        return wyhash(edge.edge_.p0_, edge.edge_.p1_);
    }

    bool operator==(const Halfedge& other) const
    {
        return edge_ == other.edge_;
    }
};

//--------------------------------------------------------------
class HashMap
{
public:
    static constexpr u32 Expand = 128U;

    HashMap(u32 capacity=Expand);
    ~HashMap();

    u32 size() const;
    u32 capacity() const;
    void reserve(u32 capacity);

    u32 exists(u32 p0, u32 p1) const;
    u32 add(u32 p0, u32 p1, u32 p2, u32 p3);
    void removeAt(u32 x);
    void remove(u32 p0, u32 p1);

    u32 begin() const;
    u32 end() const;
    u32 next(u32 x) const;

    const Halfedge& operator[](u32 x) const;
    Halfedge& operator[](u32 x);

private:
    HashMap(const HashMap&) = delete;
    HashMap& operator=(const HashMap&) = delete;

    struct Entry
    {
        u32 hash_;
        Halfedge value_;
    };

    void expand();
    u32 add_(u32 capacity, Entry* entries, const Halfedge& x, u32 hash);

    u32 capacity_;
    u32 size_;
    Entry* entries_;
};

//--------------------------------------------------------------
class HalfedgeMesh
{
public:

    HalfedgeMesh();
    ~HalfedgeMesh();

    bool validate(bool checkOpposite=true) const;

    void addTriangle(u32 p0, u32 p1, u32 p2);
    void removeTriangle(const u32 points[3]);

    u32 findHalfedge(u32 p0, u32 p1) const;
    u32 getNumHalfedges() const;
    u32 getCapacityHalfedges() const;
    u32 beginHalfedge() const;
    u32 endHalfedge() const;
    u32 nextHalfedge(u32 h) const;

    const Halfedge& getHalfedge(u32 halfedge) const;
    Halfedge& getHalfedge(u32 halfedge);
    void getFacePoints(u32 points[3], u32 p0, u32 p1) const;
private:
    HalfedgeMesh(const HalfedgeMesh&) = delete;
    HalfedgeMesh& operator=(const HalfedgeMesh&) = delete;

    HashMap halfedges_;
};
} // namespace halfedge
#endif // INC_HALFEDGE_H_

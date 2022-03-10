/**
 */
#include "halfedge.h"
#ifdef _DEBUG
#    include <stdio.h>
#endif

namespace halfedge
{
namespace
{
    static constexpr u64 wyhash_p0 = 0xa0761d65ull;
    static constexpr u64 wyhash_p1 = 0xe7037ed1ull;
    static constexpr u64 wyhash_p2 = 0x8ebc6af1ull;
    static constexpr u64 wyhash_p3 = 0x589965cdull;
    static constexpr u64 wyhash_p4 = 0x1d8e4e27ull;
    static constexpr u64 wyhash_p5 = 0xeb44accbull;

    static inline u64 wyhash_mum(const u64 A, const uint64_t B)
    {
        u64 r = A * B;
        return r - (r >> 32U);
    }
} // namespace

//--------------------------------------------------------------
u32 wyhash(u32 x0, u32 x1)
{
    u64 seed = 2685821657736338717ULL + wyhash_p5;
    seed = wyhash_mum(static_cast<u64>(x0) ^ seed, static_cast<u64>(x1) ^ wyhash_p0);
    seed = (seed ^ seed << 16) * (8U ^ wyhash_p0);
	return static_cast<u32>(seed - (seed >> 32U));
}

//--------------------------------------------------------------
HashMap::HashMap(u32 capacity)
    : capacity_(capacity)
    , size_(0)
    , entries_(nullptr)
{
    entries_ = reinterpret_cast<Entry*>(::malloc(sizeof(Entry) * capacity_));
    ::memset(entries_, 0, sizeof(Entry) * capacity_);
}

HashMap::~HashMap()
{
    ::free(entries_);
}

u32 HashMap::size() const
{
    return size_;
}

u32 HashMap::capacity() const
{
    return capacity_;
}

void HashMap::reserve(u32 capacity)
{
    //capacity <<= 1;
    if(capacity_<capacity){
        expand();
    }
}

u32 HashMap::exists(u32 p0, u32 p1) const
{
    u32 hash = wyhash(p0, p1) & 0x7FFFFFFFU;
    u32 start = hash % capacity_;
    u32 index = start;
    hash |= 0x80000000U;
    do {
        u32 h = entries_[index].hash_;
        if(hash == h) {
            if(p0 == entries_[index].value_.edge_.p0_ && p1 == entries_[index].value_.edge_.p1_) {
                return index;
            }
        }
        index = (index + 1) & (capacity_-1);
    } while(index != start);
    return Invalid;
}

u32 HashMap::add(u32 p0, u32 p1, u32 p2, u32 p3)
{
    assert(size_<capacity_);
    u32 index = exists(p0, p1);
    if(Invalid != index) {
        entries_[index].value_ = {{p0,p1},{p2,p3}};
        return index;
    }
    assert(size_<capacity_);
    ++size_;
    return add_(capacity_, entries_, {{p0,p1},{p2,p3}}, wyhash(p0, p1));
}

void HashMap::expand()
{
    u32 capacity = capacity_ << 1U;
    Entry* entries = reinterpret_cast<Entry*>(::malloc(sizeof(Entry) * capacity));
    ::memset(entries, 0, sizeof(Entry) * capacity);
    for(u32 i = begin(); i != end(); i = next(i)) {
        add_(capacity, entries, entries_[i].value_, entries_[i].hash_);
    }
    ::free(entries_);
    capacity_ = capacity;
    entries_ = entries;
}

u32 HashMap::add_(u32 capacity, Entry* entries, const Halfedge& x, u32 hash)
{
    hash &= 0x7FFFFFFFU;
    u32 start = hash % capacity;
    u32 index = start;
    do {
        u32 h = entries[index].hash_;
        if(0 == (h & 0x80000000U)) {
            assert(Invalid != x.edge_.p0_);
            assert(Invalid != x.edge_.p1_);
            assert(Invalid != x.next_.p0_);
            assert(Invalid != x.next_.p1_);
            entries[index].hash_ = hash | 0x80000000U;
            entries[index].value_ = x;
            return index;
        }
        index = (index + 1) & (capacity_-1);
    } while(index != start);
    return Invalid;
}

void HashMap::removeAt(u32 index)
{
    assert(index < capacity_);
    assert(0 != (entries_[index].hash_ & 0x80000000U));
    if(0 != (entries_[index].hash_ & 0x80000000U)) {
        --size_;
        entries_[index].hash_ = 0;
        entries_[index].value_ = {{Invalid,Invalid}, {Invalid,Invalid}};
    }
}

void HashMap::remove(u32 p0, u32 p1)
{
    u32 index = exists(p0, p1);
    if(Invalid == index){
        return;
    }
    removeAt(index);
}

u32 HashMap::begin() const
{
    for(u32 i = 0; i < capacity_; ++i) {
        if(0 != (entries_[i].hash_ & 0x80000000U)) {
            return i;
        }
    }
    return capacity_;
}

u32 HashMap::end() const
{
    return capacity_;
}

u32 HashMap::next(u32 index) const
{
    for(u32 i = index + 1; i < capacity_; ++i) {
        if(0 != (entries_[i].hash_ & 0x80000000U)) {
            return i;
        }
    }
    return capacity_;
}


const Halfedge& HashMap::operator[](u32 index) const
{
    assert(index < capacity_);
    return entries_[index].value_;
}


Halfedge& HashMap::operator[](u32 index)
{
    assert(index < capacity_);
    return entries_[index].value_;
}

//--------------------------------------------------------------
HalfedgeMesh::HalfedgeMesh()
{
}

HalfedgeMesh::~HalfedgeMesh()
{
}

bool HalfedgeMesh::validate(bool checkOpposite) const
{
    for(u32 i = beginHalfedge(); i != endHalfedge(); i = nextHalfedge(i)) {
        const Halfedge& halfedge = getHalfedge(i);
        if(Invalid == halfedge.edge_.p0_) {
            assert(false);
            return false;
        }
        if(Invalid == halfedge.edge_.p1_) {
            assert(false);
            return false;
        }
        if(checkOpposite) {
            u32 o = findHalfedge(halfedge.edge_.p1_, halfedge.edge_.p0_);
            assert(o != Invalid);
            if(Invalid == o){
                return false;
            }
            const Halfedge& opposite = getHalfedge(o);
            if(Invalid == opposite.edge_.p0_ || halfedge.edge_.p0_ != opposite.edge_.p1_) {
                assert(false);
                return false;
            }
            if(Invalid == opposite.edge_.p1_ || halfedge.edge_.p1_ != opposite.edge_.p0_) {
                assert(false);
                return false;
            }
        }
        u32 next = findHalfedge(halfedge.next_.p0_, halfedge.next_.p1_);
        if(Invalid == next){
            assert(false);
            return false;
        }
        u32 count = 1;
        u32 ids[4];
        Halfedge logs[4];
        ids[0] = i;
        logs[0] = halfedge;
        do {
            const Halfedge& tmp = getHalfedge(next);
            ids[count] = next;
            logs[count] = tmp;
            ++count;
            assert(count<=3);
            next = findHalfedge(tmp.next_.p0_, tmp.next_.p1_);
            assert(Invalid != next);
        } while(i != next);
        if(3 != count) {
            assert(false);
            return false;
        }
    }
    return true;
}

void HalfedgeMesh::addTriangle(u32 p0, u32 p1, u32 p2)
{
    assert(Invalid != p0);
    assert(Invalid != p1);
    assert(Invalid != p2);
    //assert(Invalid == halfedges_.exists(p0, p1));
    //assert(Invalid == halfedges_.exists(p1, p2));
    //assert(Invalid == halfedges_.exists(p2, p0));

    halfedges_.reserve(halfedges_.size()+3);
    u32 h0 = halfedges_.add(p0, p1, p1, p2);
    u32 h1 = halfedges_.add(p1, p2, p2, p0);
    u32 h2 = halfedges_.add(p2, p0, p0, p1);

    assert(Invalid != h0);
    assert(Invalid != h1);
    assert(Invalid != h2);

    assert(h0 != h1);
    assert(h1 != h2);
    assert(h2 != h0);

    assert(validate(false));
}

void HalfedgeMesh::removeTriangle(const u32 points[3])
{
    assert(Invalid != points[0]);
    assert(Invalid != points[1]);
    assert(Invalid != points[2]);

#ifdef _DEBUG
    u32 h0 = halfedges_.exists(points[0], points[1]);
    u32 h1 = halfedges_.exists(points[1], points[2]);
    u32 h2 = halfedges_.exists(points[2], points[0]);

    assert(Invalid != h0);
    assert(Invalid != h1);
    assert(Invalid != h2);

    assert(h0 != h1);
    assert(h1 != h2);
    assert(h2 != h0);

    assert(halfedges_[h0].next_.p0_ == halfedges_[h1].edge_.p0_);
    assert(halfedges_[h0].next_.p1_ == halfedges_[h1].edge_.p1_);
    assert(halfedges_[h1].next_.p0_ == halfedges_[h2].edge_.p0_);
    assert(halfedges_[h1].next_.p1_ == halfedges_[h2].edge_.p1_);
    assert(halfedges_[h2].next_.p0_ == halfedges_[h0].edge_.p0_);
    assert(halfedges_[h2].next_.p1_ == halfedges_[h0].edge_.p1_);

#endif

    halfedges_.remove(points[0], points[1]);
    halfedges_.remove(points[1], points[2]);
    halfedges_.remove(points[2], points[0]);
    assert(validate(false));
}

u32 HalfedgeMesh::findHalfedge(u32 p0, u32 p1) const
{
    return halfedges_.exists(p0, p1);
}

u32 HalfedgeMesh::getNumHalfedges() const
{
    return halfedges_.size();
}

u32 HalfedgeMesh::getCapacityHalfedges() const
{
    return halfedges_.capacity();
}

u32 HalfedgeMesh::beginHalfedge() const
{
    return halfedges_.begin();
}

u32 HalfedgeMesh::endHalfedge() const
{
    return halfedges_.end();
}

u32 HalfedgeMesh::nextHalfedge(u32 h) const
{
    return halfedges_.next(h);
}

const Halfedge& HalfedgeMesh::getHalfedge(u32 halfedge) const
{
    return halfedges_[halfedge];
}

Halfedge& HalfedgeMesh::getHalfedge(u32 halfedge)
{
    return halfedges_[halfedge];
}

void HalfedgeMesh::getFacePoints(u32 points[3], u32 p0, u32 p1) const
{
    u32 start,next;
    next = start = findHalfedge(p0, p1);
    u32 count = 0;
    do {
        const Halfedge& h = halfedges_[next];
        points[count] = h.edge_.p0_;
        ++count;
        next = findHalfedge(h.next_.p0_, h.next_.p1_);
    } while(next != start);
}
} // namespace halfedge

#ifndef INC_QUICKHULL_H_
#define INC_QUICKHULL_H_
/**
*/
#include <cstdint>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <utility>

namespace convexhull
{
    using s32 = int32_t;

    using u8 = uint8_t;
    using u16 = uint16_t;
    using u32 = uint32_t;
    using u64 = uint64_t;

    using f32 = float;
    using f64 = double;

    static constexpr f32 zero = 0.0f;
    static constexpr f32 one = 1.0f;

    enum class Result
    {
        Error,
        Point,
        Colinear,
        Coplanar,
        Convexhull,
    };

    //---------------------------------------------
    template<class T>
    class Array
    {
    public:
        static constexpr u32 Invalid = static_cast<u32>(-1);

        Array();
        ~Array();
        u32 capacity() const;
        u32 size() const;
        void clear();
        void reserve(u32 capacity);
        void resize(u32 size);
        void push(const T& x);
        T pop();

        void removeAt(u32 index);

        /**
         * @brief Add one item and sort items.
         * @param x 
         * @pre This is already sorted.
        */
        void insert(const T& x);

        /**
         * @brief Search the item
         * @param x 
         * @return The index if found.
         * @pre This is already sorted.
        */
        u32 binarySearch(const T& x) const;

        const T& operator[](u32 index) const;
        T& operator[](u32 index);
    private:
        Array(const Array&) = delete;
        Array& operator=(const Array&) = delete;
        inline static u32 align(u32 x)
        {
            return (x+63U)&~63U;
        }

        u32 capacity_;
        u32 size_;
        T* items_;
    };

    template<class T>
    Array<T>::Array()
        :capacity_(0)
        ,size_(0)
        ,items_(nullptr)
    {}

    template<class T>
    Array<T>::~Array()
    {
        ::free(items_);
    }

    template<class T>
    u32 Array<T>::capacity() const
    {
        return capacity_;
    }

    template<class T>
    u32 Array<T>::size() const
    {
        return size_;
    }

    template<class T>
    void Array<T>::clear()
    {
        size_ = 0;
    }

    template<class T>
    void Array<T>::reserve(u32 capacity)
    {
        capacity = align(capacity);
        if(capacity<=capacity_){
            return;
        }
        T* items = reinterpret_cast<T*>(::malloc(sizeof(T)*capacity));
        if(0<size_){
            ::memcpy(items, items_, sizeof(T)*size_);
        }
        ::free(items_);
        capacity_ = capacity;
        items_ = items;
    }

    template<class T>
    void Array<T>::resize(u32 size)
    {
        reserve(size);
        size_ = size;
    }

    template<class T>
    void Array<T>::push(const T& x)
    {
        if(capacity_<=size_){
            reserve(capacity_+64U);
        }
        items_[size_] = x;
        ++size_;
    }

    template<class T>
    T Array<T>::pop()
    {
        assert(0<size_);
        --size_;
        return items_[size_];
    }

    template<class T>
    void Array<T>::removeAt(u32 index)
    {
        assert(index<size_);
        for(u32 i=index+1; i<size_; ++i){
            items_[i-1] = items_[i];
        }
        --size_;
    }

    template<class T>
    void Array<T>::insert(const T& x)
    {
        push(x);
        assert(0<size_);
        for(u32 i=size_-1; 0<i; --i){
            if(items_[i] < items_[i-1]){
                std::swap(items_[i], items_[i-1]);
            }else{
                break;
            }
        }
    }

    template<class T>
    u32 Array<T>::binarySearch(const T& x) const
    {
        u32 first = 0;
        { // lower bound
            s32 count = static_cast<s32>(size_);
            while(0 < count) {
                u32 d = count / 2;
                u32 m = first + d;
                if(items_[m] < x){
                    first = m + 1;
                    count -= d + 1;
                } else {
                    count = d;
                }
            }
            if(size_<=first){
                return Invalid;
            }
        }
        assert(x.p0_<=items_[first].p0_);
        for(u32 i=first; i<size_; ++i){
            if(x==items_[i]){
                return i;
            }
            if(x<items_[i]){
                break;
            }
        }
        return Invalid;
    }

    template<class T>
    const T& Array<T>::operator[](u32 index) const
    {
        assert(index<size_);
        return items_[index];
    }

    template<class T>
    T& Array<T>::operator[](u32 index)
    {
        assert(index<size_);
        return items_[index];
    }

    //---------------------------------------------
    struct Vector2
    {
        f32 x_;
        f32 y_;
        Vector2 operator-() const;
        f32 operator[](u32 index) const
        {
            return reinterpret_cast<const f32*>(this)[index];
        }
        f32 lengthSqr() const;
    };

    Vector2 operator+(const Vector2& x0, const Vector2& x1);

    Vector2 operator-(const Vector2& x0, const Vector2& x1);

    Vector2 perpendicular(const Vector2& x);

    f32 dot(const Vector2& x0, const Vector2& x1);
    Vector2 tripleProduct(const Vector2& x0, const Vector2& x1, const Vector2& x2);

    //---------------------------------------------
    struct Vector3
    {
        f32 x_;
        f32 y_;
        f32 z_;

        Vector3 operator-() const;

        f32 operator[](u32 index) const
        {
            return reinterpret_cast<const f32*>(this)[index];
        }
        f32 lengthSqr() const;
        Vector3& operator+=(const Vector3& x);
    };

    Vector3 operator+(const Vector3& x0, const Vector3& x1);
    Vector3 operator-(const Vector3& x0, const Vector3& x1);
    Vector3 operator*(f32 x0, const Vector3& x1);
    Vector3 operator*(const Vector3& x0, f32 x1);

    f32 dot(const Vector3& x0, const Vector3& x1);
    Vector3 cross(const Vector3& x0, const Vector3& x1);
    Vector3 normalize(const Vector3& x);
    Vector3 normalizeSafe(const Vector3& x);
    bool isEqual(const Vector3& x0, const Vector3& x1);
    void orthonormalBasis(Vector3& binormal0, Vector3& binormal1, const Vector3& normal);

    //---------------------------------------------
    void convexhull2(Array<Vector2>& result, u32 size, Vector2* points);
    Result convexhull3(Array<u32>& indices, u32 size, const Vector3* points);

    //---------------------------------------------
    bool validate(const Array<Vector2>& points, const char* name);

    struct Validation
    {
        f32 maxDistance_;
        u32 countOuter_;
    };
    Validation validate(const Array<u32>& indices, u32 size, const Vector3* points, const char* name);
}
#endif //INC_QUICKHULL_H_


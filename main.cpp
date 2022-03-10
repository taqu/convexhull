#include "convexhull.h"
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <tuple>
#include <random>

using namespace convexhull;

static constexpr u64 Increment = 0x14057B7EF767814FULL;
static constexpr u64 Multiplier = 0x5851F42D4C957F2DULL;
inline f32 toF32_1(u32 x)
{
    static const u32 m0 = 0x3F800000U;
    static const u32 m1 = 0x007FFFFFU;
    x = m0 | (x & m1);
    return (*(f32*)&x) - 1.000000000f;
}

u64 state_ = 0xCAFEF00DD15EA5E5ULL;

namespace
{
inline u32 rotr32(u32 x, u32 r)
{
    return (x >> r) | (x << ((~r + 1) & 31U));
}
} // namespace

u32 urand()
{
    u64 x = state_;
    u32 count = static_cast<u32>(x >> 59);
    state_ = x * Multiplier + Increment;
    x ^= x >> 18;
    return rotr32(static_cast<u32>(x >> 27), count);
}

f32 frand()
{
    return toF32_1(urand());
}

void srandom(u64 seed)
{
    state_ = Increment + seed;
    urand();
}

float perturbation()
{
    return frand() * ((rand() % 2) ? 1.0f : -1.0f);
}

Vector2 jostle(Vector2 a)
{
    Vector2 b;
    b.x_ = a.x_ + perturbation();
    b.y_ = a.y_ + perturbation();
    return b;
}

std::chrono::nanoseconds proc2d(u32 numSamples, f32 scale)
{
    Vector2* points = reinterpret_cast<Vector2*>(malloc(sizeof(Vector2) * numSamples));
    for(u32 i = 0; i < numSamples; ++i) {
        points[i].x_ = frand() * scale;
        points[i].y_ = frand() * scale;
    }

    Array<Vector2> result;
    std::chrono::time_point start = std::chrono::high_resolution_clock::now();
    convexhull2(result, numSamples, points);
    std::chrono::duration duration = std::chrono::high_resolution_clock::now() - start;
    if(!validate(result, "convexhull2.csv")) {
        assert(false);
    }
    free(points);
    return std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
}

std::tuple<std::chrono::nanoseconds, f32, u32> proc3d(u32 numSamples, u32 count, f32 scale, bool skip)
{
    if(skip) {
        for(u32 i = 0; i < numSamples; ++i) {
            (void)urand();
            (void)urand();
            (void)urand();
        }
        return std::make_tuple(std::chrono::nanoseconds::zero(), 0.0f, 0);
    }
    Vector3* points = reinterpret_cast<Vector3*>(malloc(sizeof(Vector3) * numSamples));

    for(u32 i = 0; i < numSamples; ++i) {
        points[i].x_ = frand() * scale;
        points[i].y_ = frand() * scale;
        points[i].z_ = frand() * scale;
    }
    Array<u32> result;
    std::chrono::time_point start = std::chrono::high_resolution_clock::now();
    convexhull3(result, numSamples, points);
    std::chrono::duration duration = std::chrono::high_resolution_clock::now() - start;
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "./out/convexhull3_%04d.ply", count);
    auto validation = validate(result, numSamples, points, buffer);
    free(points);
    return std::make_tuple(std::chrono::duration_cast<std::chrono::nanoseconds>(duration), validation.maxDistance_, validation.countOuter_);
}

int main(void)
{
    static constexpr s32 Samples = 60000;
    //srandom(123456U);
    srandom(std::random_device()());
    static constexpr s32 Count = 1000;
    static constexpr f32 Scale = 1.0e3f;

    printf("---\nconvexhull 2d\n---\n");
    for(s32 i = 0; i < Count; ++i) {
        std::chrono::nanoseconds duration = proc2d(Samples, Scale);
        long long time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
        printf("[%d] %lld milliseconds\n", i, time);
    }

    printf("---\nconvexhull 3d\n---\n");
    std::chrono::nanoseconds maxDuration = {};
    std::chrono::nanoseconds totalDuration = {};
    f32 maxDistance = -std::numeric_limits<f32>::infinity();
    f64 totalDistance = 0.0;
    u32 maxCount = 0;
    u64 totalCount = 0;
    s32 maxDistanceIndex = -1;
    s32 maxCountIndex = -1;
    for(s32 i = 0; i < Count; ++i) {
        auto statistics = proc3d(Samples, i, Scale, i < 0);
        long long time = std::chrono::duration_cast<std::chrono::milliseconds>(std::get<0>(statistics)).count();
        printf("[%d] %lld milliseconds\n", i, time);
        if(maxDuration < std::get<0>(statistics)) {
            maxDuration = std::get<0>(statistics);
        }
        totalDuration += std::get<0>(statistics);
        if(maxDistance < std::get<1>(statistics)) {
            maxDistance = std::get<1>(statistics);
            maxDistanceIndex = i;
        }
        totalDistance += std::get<1>(statistics);
        if(maxCount < std::get<2>(statistics)) {
            maxCount = std::get<2>(statistics);
            maxCountIndex = i;
        }
        totalCount += std::get<2>(statistics);
    }

    {
        FILE* file = fopen("statistics.txt", "wb");
        fprintf(file, "cout: %d, points: %d\n", Count, Samples);
        std::chrono::nanoseconds avgDuration = totalDuration / Count;
        f64 avgDistance = totalDistance / Count;
        u64 avgCount = totalCount / Count;

        long long avgTime = std::chrono::duration_cast<std::chrono::milliseconds>(avgDuration).count();
        long long maxTime = std::chrono::duration_cast<std::chrono::milliseconds>(maxDuration).count();

        fprintf(file, "[time] avarage: %lld milliseconds, max: %lld milliseconds\n", avgTime, maxTime);
        fprintf(file, "[distance] avarage: %f, max: %f (%d)\n", avgDistance, maxDistance, maxDistanceIndex);
        fprintf(file, "[count] avarage: %zu, max: %d (%d)\n", avgCount, maxCount, maxCountIndex);
        fclose(file);
    }
    return 0;
}

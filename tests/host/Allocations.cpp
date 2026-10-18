#include "host/Allocations.hpp"

#include <atomic>
#include <cstdlib>
#include <new>

namespace {

    std::atomic<std::size_t> gCount{0}, gLive{0};

    // Each block carries its size in front, so a delete without one still knows what it frees.
    constexpr std::size_t kHeader = alignof(std::max_align_t);

    void *allocate(std::size_t size) {
        void *raw = std::malloc(size + kHeader);
        if (!raw) throw std::bad_alloc();
        *static_cast<std::size_t *>(raw) = size;
        gCount.fetch_add(1, std::memory_order_relaxed);
        gLive.fetch_add(size, std::memory_order_relaxed);
        return static_cast<char *>(raw) + kHeader;
    }

    void release(void *p) noexcept {
        if (!p) return;
        void *raw = static_cast<char *>(p) - kHeader;
        gLive.fetch_sub(*static_cast<std::size_t *>(raw), std::memory_order_relaxed);
        std::free(raw);
    }

} // namespace

namespace host {
    Allocations allocations() { return {gCount.load(), gLive.load()}; }
} // namespace host

void *operator new(std::size_t size) { return allocate(size); }
void *operator new[](std::size_t size) { return allocate(size); }
void *operator new(std::size_t size, const std::nothrow_t &) noexcept {
    try {
        return allocate(size);
    } catch (...) {
        return nullptr;
    }
}
void *operator new[](std::size_t size, const std::nothrow_t &t) noexcept { return operator new(size, t); }
void operator delete(void *p) noexcept { release(p); }
void operator delete[](void *p) noexcept { release(p); }
void operator delete(void *p, std::size_t) noexcept { release(p); }
void operator delete[](void *p, std::size_t) noexcept { release(p); }
void operator delete(void *p, const std::nothrow_t &) noexcept { release(p); }
void operator delete[](void *p, const std::nothrow_t &) noexcept { release(p); }

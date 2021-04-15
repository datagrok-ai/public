/**
    Datagrok memory cache WASM
    (c) 2021 Datagrok, Inc.

    Anatoliy Kuznetsov
*/

#include <chrono>
#include <vector>
#include <string>
#include <memory>
#include <cassert>


#ifdef __EMSCRIPTEN__

#include <emscripten/threading.h>
#include <emscripten/emscripten.h>
#include <emscripten/bind.h>
using namespace emscripten;

#endif

using namespace std;
using namespace std::chrono;

typedef int JS_heap_ptr_type; ///< JavaScript heap offset


namespace grok {

enum class CompressMethod {
    NONE,
    ZIP
};



// ------------------------------------------------------------------------

struct CacheRecord
{
    int          id;
    std::string  key;
    int          timestamp;
};


/**
    Return version string for the WASM component of the package
 */
static
std::string wasm_version()
{
    return string("0.0.1");
}

/**
    Timestamp in seconds since epoch
 */
static
int get_timestamp()
{
    system_clock::duration dtn = system_clock::now().time_since_epoch();
    return dtn.count() * system_clock::period::num / system_clock::period::den;
}

static
CacheRecord get_cache_record()
{
    CacheRecord cr;
    cr.id = 1; cr.key = "test_key"; cr.timestamp = get_timestamp();
    return cr;
}

/// ===========================================================================
///
/// Collection of memory buffers
/// (available from the JS side)
///
class ByteBufferStore
{
public:
    typedef unsigned char                                value_type;
    typedef std::vector<unsigned char>                   uchar8_vector_type;
    typedef std::vector<unique_ptr<uchar8_vector_type> > uchar8_vector_ptr_type;
public:

    /// Reset storage, delete all vectors
    void reset() { uchar8p_vectors_.resize(0); }

    /// Return total size of the stored objects
    size_t size() const noexcept { return uchar8p_vectors_.size(); }

    /// Add a new buffer, return JS legal pointer to the memory
    JS_heap_ptr_type add_buffer_js(size_t size)
        { return  (JS_heap_ptr_type) add_buffer(size);}

    /// Get buffer pointer (return 0 if not available)
    JS_heap_ptr_type get_ptr_js(size_t idx) const noexcept
    {
        value_type* p = nullptr;
        if (idx >= uchar8p_vectors_.size())
            return 0;
        p = uchar8p_vectors_[idx]->data();
        return (JS_heap_ptr_type) p;
    }

    /// Get buffer size (return 0 if not available)
    size_t get_size(size_t idx) const noexcept
    {
        if (idx >= uchar8p_vectors_.size())
            return 0;
        return uchar8p_vectors_[idx]->size();
    }


protected:

    /// Return buffer as (value_type*)
    value_type* add_buffer(size_t size)
    {
        uchar8_vector_type* buf_v;
        uchar8p_vectors_.emplace_back(buf_v = new uchar8_vector_type(size));
        return buf_v->data();
    }
protected:
    uchar8_vector_ptr_type   uchar8p_vectors_; ///< buffer vectors
};


/// ------------------------------------------------------------------------
/// Factory for ByteBufferStore
///
EMSCRIPTEN_KEEPALIVE
ByteBufferStore* create_ByteBufferStore()
{
    return new ByteBufferStore();
}


} // namespace grok


// ------------------------------------------------------------------------
// Bindings for JavaScript
//


EMSCRIPTEN_BINDINGS(GrokCacheWASM) {
    emscripten::function("wasm_version", &grok::wasm_version);
    emscripten::function("create_ByteBufferStore",
                         &grok::create_ByteBufferStore, allow_raw_pointers());

    emscripten::class_<grok::ByteBufferStore>("GrokByteBufferStore")
        .function("size",       &grok::ByteBufferStore::size)
        .function("reset",      &grok::ByteBufferStore::reset)
        .function("add_buffer", &grok::ByteBufferStore::add_buffer_js)
        .function("get_ptr",    &grok::ByteBufferStore::get_ptr_js)
        .function("get_size",   &grok::ByteBufferStore::get_size)
    ;

    enum_<grok::CompressMethod>("CompressMethod")
        .value("NONE", grok::CompressMethod::NONE)
        .value("ZIP",  grok::CompressMethod::ZIP)
        ;

    emscripten::value_object<grok::CacheRecord>("CacheRecord")
            .field("id",   &grok::CacheRecord::id)
            .field("key",  &grok::CacheRecord::key)
            .field("time", &grok::CacheRecord::timestamp)
            ;

    emscripten::function("getCacheRecord", &grok::get_cache_record);



} // EMSCRIPTEN_BINDINGS



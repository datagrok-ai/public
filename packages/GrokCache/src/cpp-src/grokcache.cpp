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

} // namespace grok

// ------------------------------------------------------------------------
// Bindings for JavaScript
//


EMSCRIPTEN_BINDINGS(GrokCacheWASM) {
    emscripten::function("wasm_version", &grok::wasm_version);

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



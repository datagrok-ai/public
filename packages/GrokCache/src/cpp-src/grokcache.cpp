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

namespace grok {

// ------------------------------------------------------------------------

/**
    Return version string for the WASM component of the package
 */
static
std::string wasm_version()
{
    return string("0.0.1");
}

} // namespace grok

// ------------------------------------------------------------------------
// Bindings for JavaScript
//


EMSCRIPTEN_BINDINGS(GrokCacheWASM) {
    emscripten::function("wasm_version", &grok::wasm_version);
}



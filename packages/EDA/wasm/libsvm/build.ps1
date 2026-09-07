# Builds SVMAPI.js + .wasm for the Datagrok EDA package from the sibling
# LIBSVM clone.
#
# The clone is never modified: LIBSVM is a single self-contained svm.cpp, so
# no worktree/patch is needed (unlike the xgboost build). Artifacts are copied
# into the EDA package only if smoke/api-test.mjs passes.
#
# Usage (from wasm/libsvm/):
#   powershell -File build.ps1 [-LibsvmDir <path>] [-EmsdkDir <path>] [-EdaWasmDir <path>]
param(
  # Defaults assume the layout Datagrok/{public/packages/EDA, ml-cpp/{libsvm,emsdk}}
  [string]$LibsvmDir  = "$PSScriptRoot\..\..\..\..\..\..\ml-cpp\libsvm",
  [string]$EmsdkDir   = "$PSScriptRoot\..\..\..\..\..\..\ml-cpp\emsdk",
  [string]$EdaWasmDir = "$PSScriptRoot\..",
  [string]$Tag        = "v3.37"
)
$ErrorActionPreference = 'Stop'

$LibsvmDir  = (Resolve-Path $LibsvmDir).Path
$EmsdkDir   = (Resolve-Path $EmsdkDir).Path
$EdaWasmDir = (Resolve-Path $EdaWasmDir).Path
$build      = Join-Path $env:TEMP "svm-wasm-build"
New-Item -ItemType Directory -Force -Path $build | Out-Null

Write-Host "libsvm clone : $LibsvmDir"
Write-Host "emsdk        : $EmsdkDir"
Write-Host "EDA wasm dir : $EdaWasmDir"
Write-Host "Expected tag : $Tag (LIBSVM_VERSION in svm.h below)"

# 1. Emscripten environment
& "$EmsdkDir\emsdk_env.ps1" | Out-Null

# 2. Common emcc flags (see README.md for the rationale of each).
$exported = '_svmTrain,_svmModelSize,_svmModelCopy,_svmLoadModel,_svmPredict,_svmFreeModel,_malloc,_free'
$common = @(
  '-O3',
  "$LibsvmDir\svm.cpp",
  "$PSScriptRoot\svm-api.cpp",
  '-I', $LibsvmDir,
  '-sMODULARIZE=1', '-sEXPORT_NAME=SVM',
  '-sALLOW_MEMORY_GROWTH=1', '-sGROWABLE_ARRAYBUFFERS=0', '-sFORCE_FILESYSTEM=1',
  "-sEXPORTED_FUNCTIONS=$exported",
  '-sEXPORTED_RUNTIME_METHODS=HEAPF32,HEAPU8'
)

# 3. Two link targets from the same sources:
#  - SVMAPI       (web,worker)      - the shipped artifact. Must NOT include the
#    node environment: its loader would contain require("node:fs") calls that
#    break webpack bundling.
#  - SVMAPI-node  (web,worker,node) - used only by smoke/api-test.mjs below;
#    never copied to wasm/.
emcc @common '-sENVIRONMENT=web,worker' -o "$build\SVMAPI.js"
if ($LASTEXITCODE -ne 0) { throw "build (web,worker) failed" }

emcc @common '-sENVIRONMENT=web,worker,node' -o "$build\SVMAPI-node.js"
if ($LASTEXITCODE -ne 0) { throw "build (node) failed" }

# 4. MANDATORY: api test against the fresh module; artifacts are copied only if
#    it passes.
node "$PSScriptRoot\smoke\api-test.mjs" "$build\SVMAPI-node.js"
if ($LASTEXITCODE -ne 0) { throw "api-test FAILED - artifacts NOT copied" }

# 5. Artifacts -> EDA package
Copy-Item "$build\SVMAPI.js", "$build\SVMAPI.wasm" $EdaWasmDir -Force
Write-Host "Artifacts copied to ${EdaWasmDir}:"
Get-Item "$EdaWasmDir\SVMAPI.js", "$EdaWasmDir\SVMAPI.wasm" |
  Format-Table Name, Length, LastWriteTime

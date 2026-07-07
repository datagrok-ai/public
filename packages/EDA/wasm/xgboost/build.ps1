# Builds XGBoostAPI.js + .wasm for the Datagrok EDA package from the sibling
# XGBoost clone.
#
# The clone itself is never modified: a temporary git worktree is created at
# the pinned tag, the minimal-build patch is applied there, and the worktree
# is removed afterwards. Artifacts are copied into the EDA package only if
# smoke/api-test.mjs passes.
#
# Usage (from wasm/xgboost/):
#   powershell -File build.ps1 [-XgboostDir <path>] [-EmsdkDir <path>] [-EdaWasmDir <path>]
param(
  # Defaults assume the layout Datagrok/{public-copy/public/packages/EDA,
  # ml-cpp/{xgboost,emsdk}}
  [string]$XgboostDir = "$PSScriptRoot\..\..\..\..\..\..\ml-cpp\xgboost",
  [string]$EmsdkDir   = "$PSScriptRoot\..\..\..\..\..\..\ml-cpp\emsdk",
  [string]$EdaWasmDir = "$PSScriptRoot\..",
  [string]$Tag        = "v3.3.0",
  [switch]$KeepWorktree
)
$ErrorActionPreference = 'Stop'

$XgboostDir = (Resolve-Path $XgboostDir).Path
$EmsdkDir   = (Resolve-Path $EmsdkDir).Path
$EdaWasmDir = (Resolve-Path $EdaWasmDir).Path
$patch      = (Resolve-Path "$PSScriptRoot\patches\xgboost-wasm-minimal-$Tag.patch").Path
$worktree   = Join-Path $env:TEMP "xgb-wasm-$Tag"
$buildDir   = Join-Path $env:TEMP "xgb-wasm-build-$Tag"

Write-Host "XGBoost clone : $XgboostDir"
Write-Host "emsdk         : $EmsdkDir"
Write-Host "EDA wasm dir  : $EdaWasmDir"
Write-Host "Tag / patch   : $Tag / $patch"

# 1. Emscripten environment
& "$EmsdkDir\emsdk_env.ps1" | Out-Null

# 2. Temporary worktree at the pinned tag + patch
if (Test-Path $worktree) {
  git -C $XgboostDir worktree remove --force $worktree
}
git -C $XgboostDir worktree add $worktree $Tag
try {
  git -C $worktree submodule update --init dmlc-core
  git -C $worktree apply $patch

  # 3. Configure + build (emcmake picks up cmake/ninja from emsdk)
  emcmake cmake -B $buildDir -S $PSScriptRoot "-DXGBOOST_SOURCE_DIR=$($worktree -replace '\\','/')"
  if ($LASTEXITCODE -ne 0) { throw "configure failed" }
  cmake --build $buildDir --parallel
  if ($LASTEXITCODE -ne 0) { throw "build failed" }

  # 4. MANDATORY: api test against the fresh module; artifacts are copied
  #    only if it passes. The node-enabled variant is used for the test; the
  #    shipped XGBoostAPI.js is web/worker-only (same wasm code, different
  #    loader environments).
  node "$PSScriptRoot\smoke\api-test.mjs" "$buildDir\XGBoostAPI-node.js"
  if ($LASTEXITCODE -ne 0) { throw "api-test FAILED - artifacts NOT copied" }

  # 5. Artifacts -> EDA package
  Copy-Item "$buildDir\XGBoostAPI.js"   $EdaWasmDir -Force
  Copy-Item "$buildDir\XGBoostAPI.wasm" $EdaWasmDir -Force
  Write-Host "Artifacts copied to ${EdaWasmDir}:"
  Get-Item "$EdaWasmDir\XGBoostAPI.js", "$EdaWasmDir\XGBoostAPI.wasm" |
    Format-Table Name, Length, LastWriteTime
} finally {
  if (-not $KeepWorktree) {
    git -C $XgboostDir worktree remove --force $worktree
  }
}

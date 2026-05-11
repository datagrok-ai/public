# Bootstrap the .kg/ Python venv and verify the graph works.
# Idempotent: safe to re-run.

$ErrorActionPreference = "Stop"

$KgDir = Split-Path -Parent $PSScriptRoot
$Venv  = Join-Path $KgDir ".venv"
$VenvPy = Join-Path $Venv "Scripts\python.exe"
$Reqs  = Join-Path $KgDir "requirements.txt"
$Db    = Join-Path $KgDir "kg.kuzu"

Write-Host "[bootstrap] kg dir: $KgDir"

# 1. Create venv if missing
if (-not (Test-Path $VenvPy)) {
    Write-Host "[bootstrap] creating venv at $Venv ..."
    py -m venv $Venv
} else {
    Write-Host "[bootstrap] venv already exists, skipping creation"
}

# 2. Install requirements (always — pip is fast on no-op)
Write-Host "[bootstrap] installing requirements ..."
& $VenvPy -m pip install --quiet --upgrade pip
& $VenvPy -m pip install --quiet -r $Reqs

# 3. Build DB if missing
if (-not (Test-Path $Db)) {
    Write-Host "[bootstrap] kg.kuzu/ missing — building from JSONL (~3 min) ..."
    & $VenvPy (Join-Path $KgDir "build.py")
} else {
    Write-Host "[bootstrap] kg.kuzu/ already exists, skipping build"
}

# 4. Smoke test
Write-Host "[bootstrap] smoke test ..."
& $VenvPy (Join-Path $KgDir "qq.py") "MATCH (n) RETURN count(n) AS nodes LIMIT 1"
& $VenvPy (Join-Path $KgDir "qq.py") --stop | Out-Null

Write-Host ""
Write-Host "[bootstrap] OK. Use:"
Write-Host "  $VenvPy $KgDir\qq.py `"<cypher>`""

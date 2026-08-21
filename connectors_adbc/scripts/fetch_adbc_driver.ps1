# Fetch pre-built Apache Arrow ADBC driver libraries into lib/.
#
#   powershell.exe -ExecutionPolicy Bypass -File scripts/fetch_adbc_driver.ps1 snowflake bigquery
#
# See fetch_adbc_driver.sh for why PyPI wheels are the distribution channel.
# `param` has to be the first statement in the file, before anything else.

param(
    [Parameter(Mandatory = $true, ValueFromRemainingArguments = $true)]
    [string[]] $Drivers
)

$ErrorActionPreference = 'Stop'

# Pin to a known-good Apache Arrow ADBC release.
$Version = '1.11.0'

$LibDir = if ($env:ADBC_LIB_DIR) { $env:ADBC_LIB_DIR } else { Join-Path $PSScriptRoot '..\lib' }
New-Item -ItemType Directory -Force -Path $LibDir | Out-Null

$tmpDir = Join-Path ([System.IO.Path]::GetTempPath()) ([System.Guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Force -Path $tmpDir | Out-Null
try {
    foreach ($name in $Drivers) {
        $driver = "adbc_driver_$name"
        $targetDll = Join-Path $LibDir "$driver.dll"
        if (Test-Path $targetDll) {
            Write-Host "already fetched: $targetDll"
            continue
        }

        $wheel = "$driver-$Version-py3-none-win_amd64.whl"
        $url = "https://files.pythonhosted.org/packages/py3/a/$driver/$wheel"
        $wheelPath = Join-Path $tmpDir 'wheel.zip'
        Write-Host "downloading $url"
        Invoke-WebRequest -Uri $url -OutFile $wheelPath
        Expand-Archive -LiteralPath $wheelPath -DestinationPath $tmpDir -Force

        # The wheel always names the file lib<driver>.so even on Windows — the file is
        # actually a PE32+ Windows DLL. Rename to <driver>.dll so
        # ManagedDriver::load_dynamic_from_name can find it via the OS loader.
        $src = Join-Path $tmpDir "$driver\lib$driver.so"
        if (-not (Test-Path $src)) {
            throw "wheel did not contain $driver/lib$driver.so"
        }
        Copy-Item -LiteralPath $src -Destination $targetDll -Force
        Write-Host "installed: $targetDll"
    }
}
finally {
    if (Test-Path $tmpDir) {
        Remove-Item -Recurse -Force $tmpDir
    }
}

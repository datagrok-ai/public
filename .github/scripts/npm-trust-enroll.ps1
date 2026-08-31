<#
.SYNOPSIS
Enroll Datagrok npm packages in npm trusted publishing (OIDC).

.DESCRIPTION
Maps every publishable package in this repo to the workflow that publishes it
and registers it with `npm trust github`. Idempotent: packages that already have
a trusted publisher, are not on the registry yet, or are not ours are skipped.

Must be run from an interactive terminal on a machine with a browser — every
`npm trust` call is gated on a one-time password, and npm waits for you to
approve a https://www.npmjs.com/auth/cli/... URL. On that first prompt, enable
"skip two-factor authentication for the next 5 minutes" on npmjs.com and the
rest of the run goes through unattended.

Requires npm >= 11.15.0 and `npm login` as a user with write access to the
@datagrok / @datagrok-libraries / @datagrok-misc orgs.
See .github/NPM_TRUSTED_PUBLISHING.md.

.EXAMPLE
.\.github\scripts\npm-trust-enroll.ps1 -DryRun

.EXAMPLE
.\.github\scripts\npm-trust-enroll.ps1 -Only '@datagrok/chem'
#>
[CmdletBinding()]
param(
    [switch] $DryRun,
    [string] $Only,
    [switch] $Refresh
)

$ErrorActionPreference = 'Stop'
# Windows PowerShell 5.1 still negotiates TLS 1.0 by default; the registry needs 1.2.
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

$Repo = 'datagrok-ai/public'
$Root = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path

$npmVersion = [version] (npm --version)
if ($npmVersion -lt [version] '11.15.0') {
    Write-Host "npm $npmVersion is too old; 'npm trust' needs >= 11.15.0. Run: npm install -g npm@latest" -ForegroundColor Red
    exit 1
}

$script:Enrolled = 0
$script:Skipped = 0
$script:Failed = 0
$script:Cached = 0

# Every npm trust call needs a fresh OTP, and the "skip 2FA" grace period is only
# five minutes - not enough for the whole fleet. Remember what is already trusted
# so a re-run resumes instead of re-verifying, and write after each change so an
# aborted run keeps its progress. Shared with the bash script.
$CacheFile = Join-Path ([Environment]::GetFolderPath('UserProfile')) '.datagrok-npm-trust.json'
$script:Cache = @{}
if ((Test-Path $CacheFile) -and -not $Refresh) {
    try {
        $loaded = Get-Content -Raw $CacheFile | ConvertFrom-Json
        if ($loaded.$Repo) {
            foreach ($p in $loaded.$Repo.PSObject.Properties) { $script:Cache[$p.Name] = $p.Value }
        }
    }
    catch { Write-Host "Ignoring unreadable cache $CacheFile" -ForegroundColor DarkGray }
}

function Save-Cache {
    $byRepo = @{}
    if (Test-Path $CacheFile) {
        try {
            $existing = Get-Content -Raw $CacheFile | ConvertFrom-Json
            foreach ($r in $existing.PSObject.Properties) {
                if ($r.Name -ne $Repo) { $byRepo[$r.Name] = $r.Value }
            }
        }
        catch { }
    }
    $byRepo[$Repo] = $script:Cache
    $byRepo | ConvertTo-Json -Depth 5 | Set-Content -Path $CacheFile -Encoding utf8
}

function Test-OnRegistry {
    param([string] $Name)
    $escaped = $Name.Replace('/', '%2f')
    try {
        $response = Invoke-WebRequest -Uri "https://registry.npmjs.org/$escaped" `
            -Method Head -UseBasicParsing -TimeoutSec 20
        return $response.StatusCode -eq 200
    }
    catch {
        return $false
    }
}

function Get-ExistingTrust {
    param([string] $Name)
    # PowerShell 5.1 wraps a redirected native command's stderr as a
    # NativeCommandError, which 'Stop' turns terminating the moment npm prints a
    # notice. Drop to 'Continue' so 2>$null can quietly swallow the probe's noise.
    $prevEap = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    $raw = (& npm trust list $Name --json 2>$null) -join "`n"
    $ok = $LASTEXITCODE -eq 0
    $ErrorActionPreference = $prevEap
    if (-not $ok -or -not $raw) { return $null }
    try { return $raw.Trim([char]0xFEFF, ' ', "`r", "`n", "`t") | ConvertFrom-Json }
    catch { return $null }
}

function Invoke-Enroll {
    param([string] $Dir, [string] $Workflow)

    $manifestPath = Join-Path $Dir 'package.json'
    if (-not (Test-Path $manifestPath)) { return }

    $manifest = Get-Content -Raw -Path $manifestPath | ConvertFrom-Json
    $name = $manifest.name
    if ([string]::IsNullOrWhiteSpace($name)) { return }
    if ($Only -and $Only -ne $name) { return }

    $owned = $name -like '@datagrok/*' -or $name -like '@datagrok-libraries/*' -or
             $name -like '@datagrok-misc/*' -or $name -eq 'datagrok-api' -or $name -eq 'datagrok-tools'
    if (-not $owned) {
        Write-Output "SKIP  $name - not a Datagrok-owned name ($Dir)"
        $script:Skipped++
        return
    }

    if ($manifest.private -eq $true) {
        Write-Output "SKIP  $name - private"
        $script:Skipped++
        return
    }

    # Cheapest check first - a cache hit costs no network and no OTP at all.
    if ($script:Cache[$name] -eq $Workflow) {
        Write-Output "CACHED $name - already trusted ($Workflow)"
        $script:Cached++
        return
    }

    if (-not (Test-OnRegistry $name)) {
        Write-Output "SKIP  $name - not on the registry yet (publish its first version by hand, then re-run)"
        $script:Skipped++
        return
    }

    $trustArgs = @('trust', 'github', $name, '--repo', $Repo, '--file', $Workflow,
                   '--allow-publish', '--yes', '--json')
    if ($DryRun) { $trustArgs += '--dry-run' }

    # Everything needed (error code, summary, approval URL) is on stdout in the
    # --json payload, so park stderr in a file rather than letting a wall of
    # expected 409 noise bury the results. It is replayed only on a real failure.
    # 'Continue' keeps the redirect from becoming a terminating NativeCommandError.
    $errFile = Join-Path ([IO.Path]::GetTempPath()) "npm-trust-$PID.err"
    $prevEap = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    $output = (& npm @trustArgs 2>$errFile) -join "`n"
    $exitCode = $LASTEXITCODE
    $ErrorActionPreference = $prevEap

    $errCode = $null
    $errText = $null
    $authUrl = $null
    if ($output) {
        try {
            # A BOM or stray whitespace makes ConvertFrom-Json throw "Invalid JSON primitive".
            $parsed = $output.Trim([char]0xFEFF, ' ', "`r", "`n", "`t") | ConvertFrom-Json
            if ($parsed.PSObject.Properties.Name -contains 'error') {
                $errCode = $parsed.error.code
                $errText = if ($parsed.error.summary) { $parsed.error.summary } else { $parsed.error.detail }
                $authUrl = $parsed.error.authUrl
            }
        }
        catch { }
        # The create path emits two concatenated JSON objects, which will never
        # parse as one, and npm masks the approval UUID in its human-readable
        # output but not in --json. Recover both by pattern instead.
        if (-not $authUrl -and $output -match 'https://www\.npmjs\.com/auth/cli/[0-9a-fA-F-]{36}') {
            $authUrl = $Matches[0]
        }
        if (-not $errCode -and $output -match '"code"\s*:\s*"([A-Z][A-Z0-9_]*)"') { $errCode = $Matches[1] }
        if (-not $errText -and $output -match '"summary"\s*:\s*"((?:[^"\\]|\\.)*)"') { $errText = $Matches[1] }
    }

    if ($exitCode -eq 0) {
        Write-Output "TRUST $name -> $Workflow"
        $script:Enrolled++
        if (-not $DryRun) {
            $script:Cache[$name] = $Workflow
            Save-Cache
        }
    }
    elseif ($errCode -eq 'EOTP') {
        Write-Host ''
        Write-Host 'npm requires two-factor approval before it will accept trust changes.' -ForegroundColor Yellow
        Write-Host 'Open this URL, tick "skip two-factor authentication for the next' -ForegroundColor Yellow
        Write-Host '5 minutes", approve, then re-run this script:' -ForegroundColor Yellow
        Write-Host ''
        Write-Host "    $authUrl" -ForegroundColor Cyan
        Write-Host ''
        Write-Host '(npm hides this UUID in its own error output; it is read from --json.)' -ForegroundColor DarkGray
        Write-Host 'The URL is single-use and expires in a few minutes.' -ForegroundColor DarkGray
        exit 3
    }
    elseif ($errCode -eq 'E409' -or ($errText -and $errText -match '(?i)already|exists|conflict')) {
        # Already trusted, but 409 does not say via which workflow - and a config
        # pointing at the wrong one silently breaks publishing. Confirm, then cache.
        $existing = Get-ExistingTrust $name
        if ($existing -and ($existing.file -ne $Workflow -or $existing.repository -ne $Repo)) {
            Write-Output ("FAIL  $name - trusted via {0}/{1}, expected {2}/{3}; fix with: npm trust revoke $name --id {4}" -f `
                $existing.repository, $existing.file, $Repo, $Workflow, $existing.id)
            $script:Failed++
        }
        else {
            Write-Output "SKIP  $name - already has a trusted publisher ($Workflow)"
            $script:Skipped++
            $script:Cache[$name] = $Workflow
            Save-Cache
        }
    }
    else {
        $reason = if ($errCode) { $errCode } else { "exit $exitCode" }
        $detail = if ($errText) { $errText } else { 'unexpected npm failure' }
        Write-Output "FAIL  $name - ${reason}: $detail"
        $script:Failed++
        if (Test-Path $errFile) { Get-Content $errFile | ForEach-Object { Write-Host "      $_" -ForegroundColor DarkGray } }
        # A failure before anything has succeeded is systemic (auth, network,
        # permissions) - bail rather than replay it across all 85 packages.
        if ($script:Enrolled -eq 0) {
            Write-Host ''
            Write-Host 'Stopping: first package failed and nothing has succeeded yet.' -ForegroundColor Red
            Write-Host 'If npm printed an auth URL above, open it, tick "skip two-factor' -ForegroundColor Red
            Write-Host 'authentication for the next 5 minutes", then re-run this script.' -ForegroundColor Red
            exit 3
        }
    }

    # Stay under the registry rate limit; pointless when nothing is being written.
    if (-not $DryRun) { Start-Sleep -Seconds 2 }
}

Invoke-Enroll (Join-Path $Root 'js-api') 'js-api.yml'
Invoke-Enroll (Join-Path $Root 'tools') 'tools.yml'
foreach ($d in Get-ChildItem -Directory (Join-Path $Root 'libraries')) { Invoke-Enroll $d.FullName 'libraries.yaml' }
foreach ($d in Get-ChildItem -Directory (Join-Path $Root 'misc')) { Invoke-Enroll $d.FullName 'misc.yaml' }
foreach ($d in Get-ChildItem -Directory (Join-Path $Root 'packages')) { Invoke-Enroll $d.FullName 'packages.yaml' }

Write-Output ''
Write-Output "enrolled=$script:Enrolled cached=$script:Cached skipped=$script:Skipped failed=$script:Failed"
if ($script:Cache.Count -gt 0) {
    Write-Host "$($script:Cache.Count) trusted package(s) remembered in $CacheFile (-Refresh to re-verify)." -ForegroundColor DarkGray
}

# PowerShell binds an unrecognised '--flag' as a positional value for -Only, which
# would otherwise match nothing and exit 0 as if there were no work to do.
if ($Only -and ($script:Enrolled + $script:Cached + $script:Skipped + $script:Failed) -eq 0) {
    Write-Host "No package matched -Only '$Only'." -ForegroundColor Red
    if ($Only.StartsWith('-')) {
        Write-Host "That looks like a flag. This is PowerShell: use -DryRun, not --dry-run." -ForegroundColor Red
    }
    exit 2
}
if ($script:Failed -gt 0) { exit 1 }

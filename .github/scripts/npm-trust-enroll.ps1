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
    [string] $Only
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

    if (-not (Test-OnRegistry $name)) {
        Write-Output "SKIP  $name - not on the registry yet (publish its first version by hand, then re-run)"
        $script:Skipped++
        return
    }

    # Ask before writing: the registry allows only one publisher per package and
    # answers a repeat create with 409. Checking also catches a config pointing at
    # the wrong workflow, which would silently break publishing.
    # PowerShell 5.1 wraps a redirected native command's stderr as a
    # NativeCommandError, which 'Stop' turns terminating the moment npm prints a
    # notice. Drop to 'Continue' so 2>$null can quietly swallow the probe's noise.
    $prevEap = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    $listRaw = (& npm trust list $name --json 2>$null) -join "`n"
    $listExit = $LASTEXITCODE
    $ErrorActionPreference = $prevEap

    if ($listExit -eq 0 -and $listRaw) {
        $existing = $null
        try { $existing = $listRaw.Trim([char]0xFEFF, ' ', "`r", "`n", "`t") | ConvertFrom-Json } catch { }
        if ($existing -and $existing.id) {
            if ($existing.file -eq $Workflow -and $existing.repository -eq $Repo) {
                Write-Output "SKIP  $name - already trusted ($($existing.repository)/$($existing.file))"
                $script:Skipped++
            }
            else {
                Write-Output ("FAIL  $name - trusted via {0}/{1}, expected {2}/{3}; fix with: npm trust revoke $name --id {4}" -f `
                    $existing.repository, $existing.file, $Repo, $Workflow, $existing.id)
                $script:Failed++
            }
            return
        }
    }

    $trustArgs = @('trust', 'github', $name, '--repo', $Repo, '--file', $Workflow,
                   '--allow-publish', '--yes', '--json')
    if ($DryRun) { $trustArgs += '--dry-run' }

    # stderr is deliberately NOT redirected: npm prints the OTP approval URL there
    # and waits, so swallowing it would look like a silent hang. Only stdout (the
    # --json payload) is captured.
    $output = (& npm @trustArgs) -join "`n"
    $exitCode = $LASTEXITCODE

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
        Write-Output "SKIP  $name - already has a trusted publisher"
        $script:Skipped++
    }
    else {
        $reason = if ($errCode) { $errCode } else { "exit $exitCode" }
        $detail = if ($errText) { $errText } else { 'see npm output above' }
        Write-Output "FAIL  $name - ${reason}: $detail"
        $script:Failed++
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
Write-Output "enrolled=$script:Enrolled skipped=$script:Skipped failed=$script:Failed"

# PowerShell binds an unrecognised '--flag' as a positional value for -Only, which
# would otherwise match nothing and exit 0 as if there were no work to do.
if ($Only -and ($script:Enrolled + $script:Skipped + $script:Failed) -eq 0) {
    Write-Host "No package matched -Only '$Only'." -ForegroundColor Red
    if ($Only.StartsWith('-')) {
        Write-Host "That looks like a flag. This is PowerShell: use -DryRun, not --dry-run." -ForegroundColor Red
    }
    exit 2
}
if ($script:Failed -gt 0) { exit 1 }

[CmdletBinding()]
param(
    [switch]$CheckOnly = $true,
    [switch]$InstallHost = $false,
    [switch]$InstallTui = $false,
    [switch]$InstallWsl = $false,
    [string]$EnvName = "nanodenovo-wsl"
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path

Write-Host "NanoDeNovo Dependency Setup ($RepoRoot)"
if ($InstallHost -or $InstallTui -or $InstallWsl) {
    $CheckOnly = $false
}

if ($CheckOnly) {
    py -3 -I (Join-Path $RepoRoot "nanodenovo.py") --doctor
    exit $LASTEXITCODE
}

if ($InstallHost) {
    Write-Host "Installing host Python lock dependencies..."
    py -3 -m pip install --upgrade -r (Join-Path $RepoRoot "requirements-lock.txt")
    if ($InstallTui) {
        py -3 -m pip install "textual>=0.85.2"
    }
}

if ($InstallWsl) {
    $wslYaml = Join-Path $RepoRoot "environment.wsl.yml"
    $wslPath = ($wslYaml -replace '\\', '/') -replace '^([A-Za-z]):', { "/mnt/$($args[0].Groups[1].Value.ToLower())" }
    Write-Host "Creating/updating WSL conda environment '$EnvName' from $wslPath..."
    wsl -e bash -lc "conda env update -n '$EnvName' -f '$wslPath' --prune || conda env create -n '$EnvName' -f '$wslPath'"
}

py -3 -I (Join-Path $RepoRoot "nanodenovo.py") --doctor
exit $LASTEXITCODE

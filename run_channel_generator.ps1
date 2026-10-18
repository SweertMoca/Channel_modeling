$ErrorActionPreference = 'Stop'
$bundledPython = 'C:\Users\M_laptop\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
if (Test-Path -LiteralPath $bundledPython) {
    & $bundledPython (Join-Path $PSScriptRoot 'Channel_model_generator.py') @args
} else {
    python (Join-Path $PSScriptRoot 'Channel_model_generator.py') @args
}

$ErrorActionPreference = 'Stop'
$paperDir = $PSScriptRoot
$texName = '0910_第一问_药材预热模型_v01_Codex'
Push-Location -LiteralPath $paperDir
try {
    New-Item -ItemType Directory -Path 'build' -Force | Out-Null
    1..2 | ForEach-Object {
        & xelatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build "$texName.tex"
        if ($LASTEXITCODE -ne 0) { throw "XeLaTeX compilation failed, pass $_." }
    }
    Copy-Item -LiteralPath "build/$texName.pdf" -Destination "$texName.pdf" -Force
    Write-Output "Compiled: $paperDir/$texName.pdf"
}
finally { Pop-Location }

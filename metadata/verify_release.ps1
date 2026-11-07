# Run from the repository root before staging a public release.
$ErrorActionPreference = 'Stop'

$oversize = Get-ChildItem -Recurse -File | Where-Object { $_.Length -ge 100MB }
if ($oversize) {
  $oversize | Select-Object FullName, Length | Format-Table -AutoSize
  throw 'At least one file is 100 MB or larger and cannot be pushed to standard GitHub without LFS.'
}

$credentialPatterns = @(
  'eyJ[a-zA-Z0-9_-]{20,}',
  '(?i)(api[_-]?key|password|secret)\s*(<-|=)\s*["'']'
)
$matches = Get-ChildItem -Recurse -File -Include *.R,*.py,*.ps1,*.md,*.yml,*.yaml,*.json |
  Select-String -Pattern $credentialPatterns
if ($matches) {
  $matches | Select-Object Path, LineNumber, Line | Format-Table -Wrap
  throw 'Potential embedded credential found. Remove it before release.'
}

Write-Host 'PASS: no tracked-sized file is >=100 MB and no credential pattern was detected.'

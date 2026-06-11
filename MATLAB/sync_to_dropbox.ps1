# sync_to_dropbox.ps1
# One-way adaptive-rate sync: Source (server) --> Destination (Dropbox)
# Safe to rerun: skips files that already match by size + modification time.
#
# Usage:
#   .\sync_to_dropbox.ps1 -Source "\\server\data" -Destination "C:\Dropbox\data"
#   .\sync_to_dropbox.ps1 -Source "\\server\data" -Destination "C:\Dropbox\data" -MaxRateMBps 2 -DryRun

param(
    [Parameter(Mandatory = $true)]
    [string]$Source,

    [Parameter(Mandatory = $true)]
    [string]$Destination,

    # Target max throughput in MB/s. The script sleeps between files to stay at or below this rate.
    [double]$MaxRateMBps = 5.0,

    # Print what would be copied without actually copying.
    [switch]$DryRun,

    # Path to log file. Default: sync_log_<timestamp>.txt next to this script.
    [string]$LogFile = ""
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
if ($LogFile -eq "") {
    $stamp = (Get-Date -Format "yyyyMMdd_HHmmss")
    $LogFile = Join-Path $PSScriptRoot "sync_log_$stamp.txt"
}

function Write-Log {
    param([string]$Message, [string]$Level = "INFO")
    $line = "[{0}] [{1}] {2}" -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), $Level, $Message
    Write-Host $line
    Add-Content -Path $LogFile -Value $line
}

# ---------------------------------------------------------------------------
# Validate paths
# ---------------------------------------------------------------------------
if (-not (Test-Path $Source)) {
    Write-Log "Source path not found: $Source" "ERROR"
    exit 1
}

if (-not $DryRun -and -not (Test-Path $Destination)) {
    New-Item -ItemType Directory -Path $Destination -Force | Out-Null
    Write-Log "Created destination directory: $Destination"
}

# ---------------------------------------------------------------------------
# Rate-limiter state (rolling window over the last ~10 files)
# ---------------------------------------------------------------------------
$windowBytes   = [System.Collections.Generic.Queue[long]]::new()
$windowSeconds = [System.Collections.Generic.Queue[double]]::new()
$windowSize    = 10   # files to average over
$maxBytesPerSec = $MaxRateMBps * 1MB

# Returns the current rolling average throughput in bytes/sec.
function Get-RollingThroughput {
    $totalBytes = 0; foreach ($b in $windowBytes)   { $totalBytes   += $b }
    $totalSecs  = 0; foreach ($s in $windowSeconds) { $totalSecs    += $s }
    if ($totalSecs -le 0) { return 0 }
    return $totalBytes / $totalSecs
}

# After copying a file, record its stats and sleep if we are too fast.
function Throttle {
    param([long]$FileBytes, [double]$ElapsedSec)

    # Update rolling window
    $windowBytes.Enqueue($FileBytes)
    $windowSeconds.Enqueue([Math]::Max($ElapsedSec, 0.001))
    while ($windowBytes.Count -gt $windowSize)   { $windowBytes.Dequeue()   | Out-Null }
    while ($windowSeconds.Count -gt $windowSize) { $windowSeconds.Dequeue() | Out-Null }

    $currentRate = Get-RollingThroughput
    if ($currentRate -le 0 -or $maxBytesPerSec -le 0) { return }

    if ($currentRate -gt $maxBytesPerSec) {
        # How long to sleep so that the window average drops to target rate?
        $totalBytes = 0; foreach ($b in $windowBytes) { $totalBytes += $b }
        $targetSecs = $totalBytes / $maxBytesPerSec
        $totalSecs  = 0; foreach ($s in $windowSeconds) { $totalSecs += $s }
        $sleepSec   = [Math]::Max($targetSecs - $totalSecs, 0)
        if ($sleepSec -gt 0.05) {
            Start-Sleep -Milliseconds ([int]($sleepSec * 1000))
        }
    }
}

# ---------------------------------------------------------------------------
# Pause/resume: checks if P was pressed; waits for another key to resume.
# ---------------------------------------------------------------------------
function Check-Pause {
    if ([Console]::KeyAvailable) {
        $key = [Console]::ReadKey($true)
        if ($key.Key -eq [ConsoleKey]::P) {
            Write-Host "`n[PAUSED] Press any key to resume, or Q to quit..." -ForegroundColor Yellow
            $resume = [Console]::ReadKey($true)
            if ($resume.Key -eq [ConsoleKey]::Q) {
                Write-Log "User quit during pause."
                exit 0
            }
            Write-Host "[RESUMED]`n" -ForegroundColor Green
        }
    }
}

# ---------------------------------------------------------------------------
# Determine whether a destination file is already up to date.
# Match criteria: same size AND modification time within 2 seconds.
# ---------------------------------------------------------------------------
function Is-UpToDate {
    param([System.IO.FileInfo]$SrcFile, [string]$DstPath)
    if (-not (Test-Path $DstPath -PathType Leaf)) { return $false }
    $dst = Get-Item $DstPath
    $sameSize = ($SrcFile.Length -eq $dst.Length)
    $timeDiff = [Math]::Abs(($SrcFile.LastWriteTimeUtc - $dst.LastWriteTimeUtc).TotalSeconds)
    return ($sameSize -and $timeDiff -le 2)
}

# ---------------------------------------------------------------------------
# Main sync loop
# ---------------------------------------------------------------------------
Write-Log ("Starting one-way sync" + $(if ($DryRun) { " [DRY RUN]" } else { "" }))
Write-Log "  Source      : $Source"
Write-Log "  Destination : $Destination"
Write-Log "  Max rate    : $MaxRateMBps MB/s"
Write-Log "  Log file    : $LogFile"

$stats = @{ Copied = 0; Skipped = 0; Failed = 0; BytesCopied = 0 }

$allFiles = Get-ChildItem -Path $Source -Recurse -File -ErrorAction SilentlyContinue

foreach ($srcFile in $allFiles) {
    # Build mirrored destination path
    $relativePath = $srcFile.FullName.Substring($Source.TrimEnd('\', '/').Length).TrimStart('\', '/')
    $dstPath      = Join-Path $Destination $relativePath

    if (Is-UpToDate $srcFile $dstPath) {
        Write-Log "SKIP  $relativePath" "SKIP"
        $stats.Skipped++
        continue
    }

    if ($DryRun) {
        Write-Log "WOULD COPY  $relativePath  ($([Math]::Round($srcFile.Length / 1MB, 2)) MB)" "DRY"
        $stats.Copied++
        continue
    }

    # Ensure destination directory exists
    $dstDir = Split-Path $dstPath -Parent
    if (-not (Test-Path $dstDir)) {
        New-Item -ItemType Directory -Path $dstDir -Force | Out-Null
    }

    # Copy with timing for rate control
    try {
        $sw = [System.Diagnostics.Stopwatch]::StartNew()
        Copy-Item -Path $srcFile.FullName -Destination $dstPath -Force
        $sw.Stop()

        # Preserve source modification time on the copy
        (Get-Item $dstPath).LastWriteTimeUtc = $srcFile.LastWriteTimeUtc

        $elapsed = $sw.Elapsed.TotalSeconds
        $mbCopied = [Math]::Round($srcFile.Length / 1MB, 2)
        $mbps     = if ($elapsed -gt 0) { [Math]::Round($srcFile.Length / 1MB / $elapsed, 2) } else { 0 }

        Write-Log "COPY  $relativePath  ($mbCopied MB at $mbps MB/s)"
        $stats.Copied++
        $stats.BytesCopied += $srcFile.Length

        Throttle -FileBytes $srcFile.Length -ElapsedSec $elapsed
    }
    catch {
        Write-Log "FAIL  $relativePath  -- $($_.Exception.Message)" "ERROR"
        $stats.Failed++
    }
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
$totalMB = [Math]::Round($stats.BytesCopied / 1MB, 2)
Write-Log "Done. Copied: $($stats.Copied) files ($totalMB MB)  |  Skipped: $($stats.Skipped)  |  Failed: $($stats.Failed)"

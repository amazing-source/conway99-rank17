#!/bin/bash
# Knuth walks with full leaf processing (build + Z4 + solver, no dedupe): estimates of Z4-survivor leaves and processing seconds
cd /c/home/researcher/Desktop/conway99-rank17
N=$1; shift
for p in "$@"; do
  powershell -NoProfile -ExecutionPolicy Bypass -File C:/home/researcher/.claude/tools/run_capped.ps1 -LimitMB 1500 -MinFreeMB 3000 -TimeoutSec 3600 -StatusSec 120 -WorkDir C:/home/researcher/Desktop/conway99-rank17/z17 -Log C:/home/researcher/Desktop/conway99-rank17/logs/knuthsolve_p$p.log -Exe C:/home/researcher/Desktop/conway99-rank17/z17/z17pipe21.exe -ArgList "out/pipe_in_p$p.txt out/knuthsolve_p$p knuth $N 5$p solve" > /dev/null
  echo "pair $p: $(tail -n 1 logs/knuthsolve_p$p.log.out)"
done

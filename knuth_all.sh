#!/bin/bash
# Knuth estimates of the stage-A tree for Z17 pairs (1 thread, capped)
cd /c/home/researcher/Desktop/conway99-rank17
for p in "$@"; do
  powershell -NoProfile -ExecutionPolicy Bypass -File C:/home/researcher/.claude/tools/run_capped.ps1 -LimitMB 1500 -MinFreeMB 3000 -TimeoutSec 1800 -StatusSec 60 -WorkDir C:/home/researcher/Desktop/conway99-rank17/z17 -Log C:/home/researcher/Desktop/conway99-rank17/logs/knuth_p$p.log -Exe C:/home/researcher/Desktop/conway99-rank17/z17/z17pipe20.exe -ArgList "out/pipe_in_p$p.txt out/knuth_p$p knuth 8192 $p" > /dev/null
  echo "pair $p: $(tail -n 1 logs/knuth_p$p.log.out)"
done

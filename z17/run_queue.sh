#!/bin/bash
# sequential PASS-2 queue: bash run_queue.sh PAIR [PAIR ...]   (uses z17pipe9.exe = lazy-compat build)
cd /c/home/researcher/Desktop/conway99-rank17/z17
for p in "$@"; do
  powershell -NoProfile -ExecutionPolicy Bypass -File "C:\Users\RESEARCHER\.claude\tools\run_capped.ps1" -LimitMB 1200 -TimeoutSec 36000 -StatusSec 1800 \
    -WorkDir "C:\Users\RESEARCHER\Desktop\conway99-rank17\z17" -Log "C:\Users\RESEARCHER\Desktop\conway99-rank17\logs\pipe4_p$p.log" \
    -Exe "C:\Users\RESEARCHER\Desktop\conway99-rank17\z17\z17pipe9.exe" -ArgList "out/pipe_in_p$p.txt out/pass2_p$p dedupe heavycap 50000000"
done

#!/bin/bash
# pass1_probe.sh PAIR...: 15-minute random-probe PASS-1 sample per pair (label PATTERN)
cd /c/home/researcher/Desktop/conway99-rank17/z17
for p in "$@"; do
  powershell -NoProfile -ExecutionPolicy Bypass -File "C:\Users\RESEARCHER\.claude\tools\run_capped.ps1" -LimitMB 1200 -TimeoutSec 900 -StatusSec 900 -WorkDir "C:\Users\RESEARCHER\Desktop\conway99-rank17\z17" -Log "C:\Users\RESEARCHER\Desktop\conway99-rank17\logs\probe_p$p.log" -Exe "C:\Users\RESEARCHER\Desktop\conway99-rank17\z17\z17pipe8.exe" -ArgList "out/pipe_in_p$p.txt out/pass1_p$p dedupe probe 0.04 $p"
done

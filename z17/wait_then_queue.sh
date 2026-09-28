#!/bin/bash
# wait_then_queue.sh WATCH_PAIR PAIR...: wait for the z17pipe run on pipe_in_pWATCH to exit, then run the queue
w=$1; shift
while powershell -NoProfile -Command "if (Get-CimInstance Win32_Process -Filter \"Name like 'z17pipe%'\" | Where-Object { \$_.CommandLine -like '*pipe_in_p$w.txt*' }) { exit 0 } else { exit 1 }"; do sleep 60; done
bash /c/home/researcher/Desktop/conway99-rank17/z17/run_queue.sh "$@"

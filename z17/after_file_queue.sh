#!/bin/bash
# after_file_queue.sh FILE PAIR...: wait until FILE exists (a finished run's .stats), then run the queue
f=$1; shift
while [ ! -f "$f" ]; do sleep 60; done
bash /c/home/researcher/Desktop/conway99-rank17/z17/run_queue.sh "$@"

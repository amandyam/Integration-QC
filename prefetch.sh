#!/usr/bin/env bash

function cmd {
local id=$1
prefetch.sh "$id"
}

source env_parallel.bash
# env_parallel cmd ::: $(cut -f 2 </storage/chentemp/atulyadm/cellrangertutorial | tail -n +2d)
env_parallel cmd ::: $(cut -f 2 </storage/chentemp/atulyadm/cellrangertutorial/fnameinfo.txt | tail -n +2)

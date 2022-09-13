#!/usr/bin/env bash
prog=$1
args=${@:2}
($1 -p 0 --localhost $args) >/dev/null & \
    ($1 -p 1 --localhost $args) >/dev/null & \
    ($1 -p 2 --localhost $args) >/dev/null & \
    ($1 -p 3 --localhost $args)

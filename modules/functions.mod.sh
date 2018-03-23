#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20180323
# Project: gpseq-centrality-estimate
# Description: this module contains common functions.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

function containsElement () {
  # from https://stackoverflow.com/a/8574392
  # Set $? to 1 if contained, 0 otherwise
  local e match="$1"; shift
  for e; do [[ "$e" == "$match" ]] && return 1; done
  return 0
}

function join_by () { local IFS="$1"; shift; echo "$*"; }

function split_by () {
  local IFS="$1"
  read -r -a out <<< "$2"
  for i in ${out[@]}; do
    echo -e "$i"
  done
}

# END ==========================================================================

################################################################################

#!/bin/bash
# Shared logging and utility functions for the JET pipeline.
# Source this file from every step script: source "${LIB_DIR}/common.sh"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

_ts() { date +"%Y-%m-%d %H:%M:%S"; }

log_info()  { echo -e "$(_ts) ${BLUE}[INFO]${NC}  $*"; }
log_step()  { echo -e "$(_ts) ${BOLD}${BLUE}[STEP]${NC}  $*"; }
log_ok()    { echo -e "$(_ts) ${GREEN}[OK]${NC}    $*"; }
log_warn()  { echo -e "$(_ts) ${YELLOW}[WARN]${NC}  $*"; }
log_error() { echo -e "$(_ts) ${RED}[ERROR]${NC} $*" >&2; }

# die <message> [exit_code]
die() {
    log_error "$1"
    exit "${2:-1}"
}

# require_file <path> <description>
require_file() {
    local path="$1" desc="$2"
    [ -f "$path" ] || die "${desc} not found: ${path}"
}

# require_dir <path> <description>
require_dir() {
    local path="$1" desc="$2"
    [ -d "$path" ] || die "${desc} not found: ${path}"
}

# run_cmd <description> <command...>
# Runs a command, timing it, logging start/end, and dying on failure.
run_cmd() {
    local desc="$1"; shift
    log_info "Running: ${desc}"
    local start=$(date +%s)
    if "$@"; then
        local end=$(date +%s)
        log_ok "${desc} (took $((end-start))s)"
    else
        local rc=$?
        die "${desc} FAILED (exit code ${rc})" "${rc}"
    fi
}

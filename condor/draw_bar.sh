# shellcheck shell=bash
_BAR_COLORS=(green blue red pink purple orange yellow cyan white)
_LAST_BAR_COLOR=""
_COLOR_QUEUE=()
BAR_COLOR=""
_BAR_START_TIME=""
_TERM_COLS=80

start_bar_timer() { _BAR_START_TIME=$(date +%s); }

# Real terminal width, falling back when stdout isn't a tty (e.g. Condor job
# output redirected to a log file). Sets _TERM_COLS (global, same convention
# as BAR_COLOR above) rather than a $(...)-captured return value -- capturing
# would redirect this function's own fd 1 to the capture pipe and break the
# -t 1 check below. Also avoids `local -n` (bash 4.3+), since macOS ships
# bash 3.2.
terminal_width() {
  _TERM_COLS=80
  if [[ -t 1 ]]; then
    local size cols
    size=$(stty size </dev/tty 2>/dev/null) || size=""
    cols="${size#* }"
    if [[ -n "$cols" && "$cols" -gt 0 ]]; then
      _TERM_COLS="$cols"
    fi
  fi
}

# Sets BAR_COLOR by drawing from a shuffled deck; refills when exhausted,
# guaranteeing every color appears once per cycle with no consecutive repeats.
pick_bar_color() {
  if ((${#_COLOR_QUEUE[@]} == 0)); then
    local shuffled=("${_BAR_COLORS[@]}")
    local n=${#shuffled[@]} i j tmp
    for ((i = n - 1; i > 0; i--)); do
      j=$((RANDOM % (i + 1)))
      tmp="${shuffled[i]}"
      shuffled[i]="${shuffled[j]}"
      shuffled[j]="${tmp}"
    done
    if [[ "${shuffled[0]}" == "$_LAST_BAR_COLOR" && n -gt 1 ]]; then
      tmp="${shuffled[0]}"
      shuffled[0]="${shuffled[1]}"
      shuffled[1]="${tmp}"
    fi
    _COLOR_QUEUE=("${shuffled[@]}")
  fi
  BAR_COLOR="${_COLOR_QUEUE[0]}"
  _LAST_BAR_COLOR="$BAR_COLOR"
  _COLOR_QUEUE=("${_COLOR_QUEUE[@]:1}")
}

draw_bar() {
  local color="$1" label="$2" current="$3" total="$4"
  local pct=$((current >= total ? 100 : current * 100 / total))
  local ansi_color
  case "$color" in
    green) ansi_color='\033[32m' ;;
    blue) ansi_color='\033[34m' ;;
    red) ansi_color='\033[31m' ;;
    pink) ansi_color='\033[35m' ;;
    purple) ansi_color='\033[95m' ;;
    orange) ansi_color='\033[33m' ;;
    yellow) ansi_color='\033[93m' ;;
    cyan) ansi_color='\033[96m' ;;
    white) ansi_color='\033[97m' ;;
    *) ansi_color='\033[0m' ;;
  esac
  local reset='\033[0m'
  local grey='\033[90m'
  local count_str="${current}/${total}"
  local rate_str="" time_str="  --:-- ETA"
  if ((current >= total)); then time_str=""; fi
  if [[ -n "$_BAR_START_TIME" && "$current" -gt 0 ]]; then
    local elapsed=$(($(date +%s) - _BAR_START_TIME))
    if ((elapsed >= 1)); then
      rate_str=$(awk "BEGIN{printf \"  %4.1f/s\", $current/$elapsed}")
      if ((current >= total)); then
        time_str=$(printf "  %02d:%02d" "$((elapsed / 60))" "$((elapsed % 60))")
      elif ((pct >= 5)); then
        local eta_secs=$((elapsed * (total - current) / current))
        time_str=$(printf "  %02d:%02d ETA" "$((eta_secs / 60))" "$((eta_secs % 60))")
      fi
    fi
  fi
  terminal_width

  # Progressively drop cosmetics as the terminal narrows -- rate first,
  # then ETA/elapsed, then the label -- so the bar itself is the last
  # thing to give up space. If even a bare minimum-width bar won't fit,
  # drop the bar too and fall back to plain "current/total pct%".
  local min_bar=10
  local use_label="$label" use_rate="$rate_str" use_time="$time_str" width
  local stage
  for stage in 0 1 2 3; do
    case "$stage" in
      1) use_rate="" ;;
      2) use_time="" ;;
      3) use_label="" ;;
    esac
    local overhead=$((2 + ${#use_label} + 2 + 2 + 1 + 2 + ${#count_str} + 2 + 4 + ${#use_rate} + ${#use_time} + 1))
    width=$((_TERM_COLS - overhead))
    ((width >= min_bar)) && break
  done

  if ((width < min_bar)); then
    printf "\r  %s  %3d%%%s%s\033[K" "$count_str" "$pct" "$use_rate" "$use_time"
    return
  fi

  local filled=$((current >= total ? width : current * width / total))
  local empty=$((width - filled))
  local filled_str="" empty_str=""
  ((filled > 0)) && filled_str="$(printf '%0.s█' $(seq 1 $filled))"
  ((empty > 0)) && empty_str="$(printf '%0.s░' $(seq 1 $empty))"
  printf "\r  %s  [${ansi_color}%s${reset}${grey}%s${reset}]  %s  %3d%%%s%s\033[K" \
    "$use_label" "$filled_str" "$empty_str" "$count_str" "$pct" "$use_rate" "$use_time"
}

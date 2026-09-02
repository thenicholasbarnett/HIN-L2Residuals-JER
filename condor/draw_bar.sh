# shellcheck shell=bash
_BAR_COLORS=(green blue red pink purple orange yellow cyan white)
_LAST_BAR_COLOR=""
_COLOR_QUEUE=()
BAR_COLOR=""
_BAR_START_TIME=""
_TERM_COLS=80
_BAR_HIST_CUR=()
_BAR_HIST_TIME=()
_LAST_RATE_STR=""
_LAST_ETA_STR=""
_LAST_ETA_TIME=""

trap 'printf "\r\033[2K"' WINCH

start_bar_timer() {
  _BAR_START_TIME=$(date +%s)
  _BAR_HIST_CUR=()
  _BAR_HIST_TIME=()
  _LAST_RATE_STR=""
  _LAST_ETA_STR=""
  _LAST_ETA_TIME=""
}

# scales with terminal width
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
  
  _TERM_COLS=$((_TERM_COLS - 1))
}

# color of progress bar is randomized
# every color appears once before any appear twice
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
  local rate_str="$_LAST_RATE_STR" time_str="  --:-- ETA"
  if ((current >= total)); then time_str=""; fi
  if [[ -n "$_BAR_START_TIME" && "$current" -gt 0 ]]; then
    local now
    now=$(date +%s)
    local elapsed=$((now - _BAR_START_TIME))

    # rate is a trailing window over the last ~5% of total, not a
    # since-start running average -- a real slowdown would get smoothed
    # away by a full-run average instead of showing up in the number
    _BAR_HIST_CUR+=("$current")
    _BAR_HIST_TIME+=("$now")
    local window=$((total * 5 / 100))
    ((window < 1)) && window=1
    local threshold=$((current - window))
    while ((${#_BAR_HIST_CUR[@]} > 1)) && ((_BAR_HIST_CUR[0] < threshold)); do
      _BAR_HIST_CUR=("${_BAR_HIST_CUR[@]:1}")
      _BAR_HIST_TIME=("${_BAR_HIST_TIME[@]:1}")
    done
    local win_elapsed=$((now - _BAR_HIST_TIME[0]))
    local win_delta=$((current - _BAR_HIST_CUR[0]))
    if ((win_elapsed >= 1 && win_delta > 0)); then
      _LAST_RATE_STR=$(awk "BEGIN{printf \"  %4.1f/s\", $win_delta/$win_elapsed}")
      rate_str="$_LAST_RATE_STR"
    fi

    if ((elapsed >= 1)); then
      if ((current >= total)); then
        time_str=$(printf "  %02d:%02d" "$((elapsed / 60))" "$((elapsed % 60))")
      elif ((pct >= 5)); then
        # ETA text only refreshes once a second -- the windowed rate above
        # still updates every call, this just stops the countdown from
        # jittering with it
        if [[ -z "$_LAST_ETA_TIME" || $((now - _LAST_ETA_TIME)) -ge 1 ]]; then
          local eta_secs
          if ((win_elapsed >= 1 && win_delta > 0)); then
            eta_secs=$((win_elapsed * (total - current) / win_delta))
          else
            eta_secs=$((elapsed * (total - current) / current))
          fi
          _LAST_ETA_STR=$(printf "  %02d:%02d ETA" "$((eta_secs / 60))" "$((eta_secs % 60))")
          _LAST_ETA_TIME="$now"
        fi
        time_str="$_LAST_ETA_STR"
      fi
    fi
  fi
  terminal_width

  # progressively drop cosmetics as terminal narrows
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

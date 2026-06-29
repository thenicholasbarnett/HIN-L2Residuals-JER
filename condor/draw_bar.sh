_BAR_COLORS=(green blue red pink purple orange yellow cyan white)
_LAST_BAR_COLOR=""
_COLOR_QUEUE=()
BAR_COLOR=""

# Sets BAR_COLOR by drawing from a shuffled deck; refills when exhausted,
# guaranteeing every color appears once per cycle with no consecutive repeats.
pick_bar_color() {
    if (( ${#_COLOR_QUEUE[@]} == 0 )); then
        local shuffled=("${_BAR_COLORS[@]}")
        local n=${#shuffled[@]} i j tmp
        for (( i = n-1; i > 0; i-- )); do
            j=$(( RANDOM % (i+1) ))
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
    local width=40
    local filled=$(( current >= total ? width : current * width / total ))
    local empty=$(( width - filled ))
    local pct=$(( current >= total ? 100 : current * 100 / total ))
    local ansi_color
    case "$color" in
        green)   ansi_color='\033[32m'  ;;
        blue)    ansi_color='\033[34m'  ;;
        red)     ansi_color='\033[31m'  ;;
        pink)    ansi_color='\033[35m'  ;;
        purple)  ansi_color='\033[95m'  ;;
        orange)  ansi_color='\033[33m'  ;;
        yellow)  ansi_color='\033[93m'  ;;
        cyan)    ansi_color='\033[96m'  ;;
        white)   ansi_color='\033[97m'  ;;
        *)       ansi_color='\033[0m'   ;;
    esac
    local reset='\033[0m'
    local grey='\033[90m'
    local filled_str="" empty_str=""
    (( filled > 0 )) && filled_str="$(printf '%0.s█' $(seq 1 $filled))"
    (( empty > 0 ))  && empty_str="$(printf '%0.s░' $(seq 1 $empty))"
    printf "\r  %-24s [${ansi_color}%s${reset}${grey}%s${reset}] %d/%d (%d%%)" \
        "$label" "$filled_str" "$empty_str" "$current" "$total" "$pct"
}

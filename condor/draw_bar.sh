_BAR_COLORS=(green blue red pink purple orange yellow cyan white)
_LAST_BAR_COLOR=""
BAR_COLOR=""

# Sets BAR_COLOR to a random color different from the previous call
pick_bar_color() {
    local color tries=0
    while (( tries++ < 20 )); do
        color="${_BAR_COLORS[$(( RANDOM % ${#_BAR_COLORS[@]} ))]}"
        [[ "$color" != "$_LAST_BAR_COLOR" ]] && break
    done
    _LAST_BAR_COLOR="$color"
    BAR_COLOR="$color"
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

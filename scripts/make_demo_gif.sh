#!/usr/bin/env bash
# Trim + encode the record_demo.cjs recording into docs/demo.gif. Needs ffmpeg + gifski.
# Usage: scripts/make_demo_gif.sh [INPUT_WEBM] [OUTPUT_GIF]   (defaults: newest webm -> docs/demo.gif)
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
IN="${1:-$(ls -t /tmp/demo_rec/*.webm 2>/dev/null | head -1)}"
OUT="${2:-$REPO/docs/demo.gif}"

# record_demo.cjs writes the exact top->bottom scroll window to markers.env beside the webm;
# source it so the GIF starts at the page top with no manual -ss/-t guessing (SS/DUR env still win).
MARKERS="$(dirname "${IN:-/tmp/demo_rec/x}")/markers.env"
[ -f "$MARKERS" ] && . "$MARKERS"
SS="${SS:-${GIF_SS:-2.4}}"
DUR="${DUR:-${GIF_DUR:-5.2}}"
FPS=16
# Crop at y=0 (not centered): the topbar is position:sticky, so the helix header must stay in frame.
VF="crop=1280:536:0:0,scale=1200:502:flags=lanczos"

if [ -z "${IN:-}" ] || [ ! -f "$IN" ]; then
  echo "No input webm found. Record first (scripts/record_demo.cjs) or pass a path." >&2
  exit 1
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "Input : $IN"
echo "Output: $OUT"
ffmpeg -loglevel error -ss "$SS" -t "$DUR" -i "$IN" -vf "$VF,fps=$FPS" "$TMP/f_%04d.png"
gifski --fps "$FPS" --width 1200 --height 502 --quality 84 -o "$OUT" "$TMP"/f_*.png

echo "Done."
ffprobe -v error -select_streams v:0 -show_entries stream=width,height,nb_frames -of default=noprint_wrappers=1 "$OUT"
du -h "$OUT"

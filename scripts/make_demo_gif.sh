#!/usr/bin/env bash
# Convert the Playwright recording (scripts/record_demo.cjs) into the README hero GIF:
# a cinematic 2.39:1 widescreen (1200x502) center-letterbox crop of the scrolling dashboard.
#
# Usage:
#   scripts/record_demo.cjs  ->  /tmp/demo_rec/<page>.webm     (record first; see that file's header)
#   scripts/make_demo_gif.sh [INPUT_WEBM] [OUTPUT_GIF]
# Defaults: newest *.webm in /tmp/demo_rec  ->  docs/demo.gif
#
# Requires: ffmpeg + gifski (both on PATH).
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
IN="${1:-$(ls -t /tmp/demo_rec/*.webm 2>/dev/null | head -1)}"
OUT="${2:-$REPO/docs/demo.gif}"

# Trim to the meaningful window (skip the blank load; keep structure -> mutations -> rankings).
SS=2.4      # start: the 3D structure has painted and key facts are on screen
DUR=5.2     # ends on the Target Priority rankings (avoids the empty-left frames at page bottom)
FPS=16      # ~16.7fps reference cadence; 5.2s x 16 ~= 83 frames
# Cinematic 2.39:1 band: from the 1280x860 capture take a centered 1280x536 strip, scale to 1200x502.
VF="crop=1280:536:0:162,scale=1200:502:flags=lanczos"

if [ -z "${IN:-}" ] || [ ! -f "$IN" ]; then
  echo "No input webm found. Record first (scripts/record_demo.cjs) or pass a path." >&2
  exit 1
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "Input : $IN"
echo "Output: $OUT"
ffmpeg -loglevel error -ss "$SS" -t "$DUR" -i "$IN" -vf "$VF,fps=$FPS" "$TMP/f_%04d.png"
gifski --fps "$FPS" --width 1200 --height 502 --quality 88 -o "$OUT" "$TMP"/f_*.png

echo "Done."
ffprobe -v error -select_streams v:0 -show_entries stream=width,height,nb_frames -of default=noprint_wrappers=1 "$OUT"
du -h "$OUT"

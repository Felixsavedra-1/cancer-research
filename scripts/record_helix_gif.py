#!/usr/bin/env python3
from __future__ import annotations

import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parent.parent
OUT_GIF = ROOT / "vr01-helix.gif"

CSS_W = 256
CSS_H = 160
SCALE = 2
FPS = 16
# drawHelix spins at spin=t*0.9; one full revolution is a seamless loop.
PERIOD = 2 * math.pi / 0.9
FRAMES = 96
OUT_W = CSS_W * SCALE
MAX_COLORS = 256

PAGE = """<!doctype html><html><head><meta charset="utf-8"><style>
  html,body{margin:0;background:#0A0A0B;}
  #stage{width:256px;height:160px;}
  canvas{
    width:256px;height:160px;display:block;
    border:1px solid #26262B;border-radius:2px;
    background:linear-gradient(180deg,#111113,#0A0A0B);
  }
</style></head><body>
<div id="stage"><canvas id="c" width="512" height="320"></canvas></div>
<script>
  var CRIM='178,58,58', GREY='140,140,147', BONE='#E8E3D6';

  function grid(ctx,W,H){
    ctx.strokeStyle='rgba(38,38,43,0.85)'; ctx.lineWidth=1;
    for(var x=0;x<=W;x+=W/8){ctx.beginPath();ctx.moveTo(x,0);ctx.lineTo(x,H);ctx.stroke();}
    for(var y=0;y<=H;y+=H/8){ctx.beginPath();ctx.moveTo(0,y);ctx.lineTo(W,y);ctx.stroke();}
  }
  function dot(ctx,x,y,r,fill){ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ctx.fillStyle=fill;ctx.fill();}

  function drawHelix(ctx,W,H,t){
    grid(ctx,W,H);
    var cx=W*0.5, top=H*0.12, bot=H*0.88, span=bot-top;
    var amp=W*0.20, turns=2.2, N=140, spin=t*0.9;
    ctx.strokeStyle='rgba('+GREY+',0.28)';ctx.lineWidth=1;
    ctx.beginPath();ctx.moveTo(cx,0);ctx.lineTo(cx,H);ctx.stroke();
    function strandPts(phase){
      var pts=[];
      for(var i=0;i<=N;i++){
        var s=i/N, y=top+s*span, th=s*turns*Math.PI*2+spin+phase;
        pts.push({x:cx+amp*Math.sin(th), y:y, z:Math.cos(th)});
      }
      return pts;
    }
    var A=strandPts(0), B=strandPts(Math.PI), segs=[];
    for(var i=0;i<N;i++){
      segs.push({k:'s',c:CRIM,p:A[i],q:A[i+1],z:(A[i].z+A[i+1].z)/2});
      segs.push({k:'s',c:GREY,p:B[i],q:B[i+1],z:(B[i].z+B[i+1].z)/2});
    }
    for(var r=0;r<=N;r+=7){ segs.push({k:'r',c:CRIM,p:A[r],q:B[r],z:(A[r].z+B[r].z)/2}); }
    segs.sort(function(m,n){return m.z-n.z;});
    for(var u=0;u<segs.length;u++){
      var sg=segs[u], depth=(sg.z+1)/2;
      if(sg.k==='r'){ ctx.strokeStyle='rgba('+sg.c+','+(0.14+0.30*depth).toFixed(3)+')'; ctx.lineWidth=1; }
      else{ ctx.strokeStyle='rgba('+sg.c+','+(0.22+0.62*depth).toFixed(3)+')'; ctx.lineWidth=1+1.3*depth; }
      ctx.beginPath();ctx.moveTo(sg.p.x,sg.p.y);ctx.lineTo(sg.q.x,sg.q.y);ctx.stroke();
    }
    for(var j=0;j<=N;j+=7){
      if(A[j].z>0.1){ dot(ctx,A[j].x,A[j].y,1.6+1.4*A[j].z,'rgba('+CRIM+','+(0.45+0.45*A[j].z).toFixed(3)+')'); }
      if(B[j].z>0.1){ dot(ctx,B[j].x,B[j].y,1.3+1.1*B[j].z,BONE); }
    }
  }

  var ctx=document.getElementById('c').getContext('2d');
  ctx.scale(2,2);
  window.renderHelix=function(t){
    ctx.clearRect(0,0,256,160);
    drawHelix(ctx,256,160,t);
  };
  window.renderHelix(0);
</script></body></html>"""


def capture_frames(frame_dir: Path) -> int:
    with sync_playwright() as pw:
        browser = pw.chromium.launch(headless=True)
        page = browser.new_page(
            viewport={"width": CSS_W, "height": CSS_H},
            device_scale_factor=SCALE,
        )
        page.set_content(PAGE, wait_until="networkidle")
        canvas = page.locator("#c")
        for k in range(FRAMES):
            t = (k / FRAMES) * PERIOD
            page.evaluate(f"window.renderHelix({t})")
            canvas.screenshot(path=str(frame_dir / f"frame_{k:04d}.png"))
        browser.close()
    return FRAMES


def build_gif(frame_dir: Path) -> None:
    palette = frame_dir / "palette.png"
    vf = f"fps={FPS},scale={OUT_W}:-1:flags=lanczos"
    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-framerate",
            str(FPS),
            "-i",
            str(frame_dir / "frame_%04d.png"),
            "-vf",
            f"{vf},palettegen=max_colors={MAX_COLORS}:stats_mode=diff",
            str(palette),
        ],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-framerate",
            str(FPS),
            "-i",
            str(frame_dir / "frame_%04d.png"),
            "-i",
            str(palette),
            "-lavfi",
            f"{vf}[x];[x][1:v]paletteuse=dither=sierra2_4a",
            "-loop",
            "0",
            str(OUT_GIF),
        ],
        check=True,
        capture_output=True,
    )


def main() -> None:
    if not shutil.which("ffmpeg"):
        sys.exit("ffmpeg not found on PATH.")
    frame_dir = Path(tempfile.mkdtemp(prefix="vr01_helix_gif_"))
    try:
        n = capture_frames(frame_dir)
        print(f"Captured {n} frames → assembling GIF…")
        build_gif(frame_dir)
        size_mb = OUT_GIF.stat().st_size / 1e6
        print(f"Wrote {OUT_GIF.name} ({size_mb:.2f} MB)")
    finally:
        shutil.rmtree(frame_dir, ignore_errors=True)


if __name__ == "__main__":
    main()

/**
 * Record the README hero GIF as a .webm via Playwright (convert with ffmpeg + gifski).
 * Serve the repo first, then run headless with ANGLE for reliable WebGL:
 *   python3 -m http.server 8766
 *   HEADLESS=1 CHROME_PATH=<chromium> NODE_PATH=<playwright> node scripts/record_demo.cjs
 * Env overrides: DEMO_URL, OUT_DIR, HEADLESS, CHROME_PATH.
 */
const { chromium } = require('playwright');

const URL = process.env.DEMO_URL || 'http://localhost:8766/cancer-explorer.html';
const OUT_DIR = process.env.OUT_DIR || '/tmp/demo_rec';
const W = 1280;
const H = 860;

(async () => {
  const headless = process.env.HEADLESS === '1';
  const launchOpts = headless
    ? { headless: true, args: ['--use-gl=angle', '--use-angle=default', '--enable-webgl'] }
    : { headless: false };
  if (process.env.CHROME_PATH) launchOpts.executablePath = process.env.CHROME_PATH;
  const browser = await chromium.launch(launchOpts);

  const context = await browser.newContext({
    viewport: { width: W, height: H },
    deviceScaleFactor: 2,
    recordVideo: { dir: OUT_DIR, size: { width: W, height: H } },
  });
  const page = await context.newPage();

  const t0 = Date.now();
  await page.goto(URL, { waitUntil: 'load', timeout: 60000 });
  await page.locator('#load').click();

  try {
    await page.waitForFunction(
      () => {
        const canvas = document.querySelector('#viewer canvas');
        const hot = document.querySelector('#hotspots');
        const about = document.querySelector('#about');
        const titled =
          (document.querySelector('#proteinTitle')?.textContent || '') !== 'Protein viewer';
        const hotReady = hot && !/Load a protein/i.test(hot.textContent || '');
        const aboutReady = about && !/Load a protein/i.test(about.textContent || '');
        return canvas && hotReady && aboutReady && titled;
      },
      { timeout: 60000 }
    );
  } catch (e) {
    console.error('  (content did not fully populate in time — continuing)');
  }

  await page.waitForTimeout(1500);

  console.log('SCROLL_START_SEC=' + ((Date.now() - t0) / 1000).toFixed(2));

  await page.evaluate(async () => {
    function pickScroller() {
      const cands = [document.scrollingElement, document.body, document.documentElement];
      let best = document.scrollingElement || document.body;
      let over = 0;
      for (const e of cands) {
        if (!e) continue;
        const o = e.scrollHeight - e.clientHeight;
        if (o > over) {
          over = o;
          best = e;
        }
      }
      return best;
    }
    // Override the page's smooth scroll so our own ease curve drives the motion.
    document.documentElement.style.scrollBehavior = 'auto';
    document.body.style.scrollBehavior = 'auto';

    const el = pickScroller();
    const max = el.scrollHeight - el.clientHeight;
    const duration = 5000;
    const ease = (t) => 0.5 - 0.5 * Math.cos(Math.PI * t);
    await new Promise((resolve) => {
      const start = performance.now();
      function frame(now) {
        const t = Math.min(1, (now - start) / duration);
        el.scrollTop = max * ease(t);
        if (t < 1) requestAnimationFrame(frame);
        else resolve();
      }
      requestAnimationFrame(frame);
    });
  });

  await page.waitForTimeout(800);

  const video = page.video();
  await context.close();
  const path = video ? await video.path() : null;
  await browser.close();
  console.log('VIDEO_PATH=' + path);
})().catch((e) => {
  console.error(e);
  process.exit(1);
});

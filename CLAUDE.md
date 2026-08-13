# DiDestimands — repo notes for Claude

This repo is served at https://didestimands.app/ via GitHub Pages
(`.github/workflows/static.yml`; any push to `main` redeploys in ~1 min). Served files:
`index.html` (the whole applet, a single self-contained HTML+CSS+JS file) and
`did-companion-replication.zip` (the replication package). The editable script sources are in
`replication/` — edit those, then re-embed them into `index.html`'s code boxes and rebuild the zip.

## Hard rules
- Keep this build **pre-SDID**: no SDID tab, no `synthdid` in the scripts or the zip. (SDID lives only
  on the separate Cornell WordPress build, not here.)
- `index.html` must stay **fully self-contained** — inline all CSS/JS; no external assets.
- The applet's DGP (`simulate`, `rng` mulberry32, `qnorm` Acklam, draw order) is **bit-identical** to
  the panel the replication scripts read. If you change the DGP, regenerate the reference panel
  (`did_lab_panel.csv`) and re-verify all three scripts before pushing.
- Standard errors are CRV1 clustered by matched pair, **fixest-exact**: the factors `G/(G-1)` and
  `(N-1)/(N-K)` (nested K = period FE + D) go **inside** the sqrt — in both the applet JS (`crSE`,
  levels delta-method, PPML) and the scripts' `se_theta`.
- Estimand names: "typical-unit %", "population-total %". Cornell red `#B31B1B` = true DGP effects only.
- Each script lives in **two places** — `replication/…` (and the zip) *and* the applet code box
  (`<div class="codebox" id="code-py"/"code-r"/"code-jl">`). Keep them in sync.

## Applet map (`index.html`, inside `<script>`)
- **DGP:** `simulate()`, `rng()`, `randn()`, `qnorm()`, `sigmaForShareN()`
- **Inference:** `estimateWithSE()`, `crSE()`, `ppml()`
- **Charts:** `drawPoints()` and the Tab-2 heatmap/trajectory
- **Embedded scripts:** the `code-py` / `code-r` / `code-jl` codeboxes
- **Tabs:** sections `t1`, `t3` (Levels-OLS bias / "Tab 2"), `tlit`, `trep`

## Verify before every push
- `node verify.js` — evaluates the applet headlessly and prints the default levels CI
  (expect `[-4.4694%, -3.2529%]`). (Harness in the handoff doc, Appendix 1.)
- Python replication: `uv run replication/did_lab_replication.py`.
- R: open `replication/did_lab_replication.R` (it pins `renv::use("fixest@0.11.2")`) and run.

## Deploy
`git pull --no-rebase` → `git add -A` → `git commit -m "…"` → `git push`. GitHub Pages redeploys in
~1 min; hard-refresh didestimands.app (Cmd/Ctrl-Shift-R).

## When editing a replication script
Edit `replication/<file>`, then ask Claude to **re-embed it into `index.html`'s matching code box and
rebuild `did-companion-replication.zip`** (HTML-escape the script text into the box; re-zip the
`replication/` contents to the repo root). Then verify and push.

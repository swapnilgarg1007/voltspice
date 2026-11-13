const $ = (s) => document.querySelector(s);
const COLORS = ["#2563eb", "#dc2626", "#16a34a", "#d97706", "#9333ea", "#0891b2", "#db2777", "#65a30d"];
let data = null;
let signals = [];          // [{name, y: column index, y2: phase column index (ac only)}]
const hidden = new Set();  // signal indices

const SI = ["f", "p", "n", "\u00b5", "m", "", "k", "M", "G"];
function eng(v) {
  if (v === 0 || !isFinite(v)) return String(v);
  const e = Math.max(-5, Math.min(3, Math.floor(Math.log10(Math.abs(v)) / 3)));
  return parseFloat((v / Math.pow(10, 3 * e)).toPrecision(3)) + SI[e + 5];
}

function niceTicks(lo, hi, n = 6) {
  const range = hi - lo || 1, raw = range / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw);
  const ticks = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + step * 1e-9; t += step) ticks.push(t);
  return ticks;
}

async function loadExamples() {
  const ex = await (await fetch("/api/examples")).json();
  const sel = $("#examples");
  for (const name of Object.keys(ex)) sel.add(new Option(name.replace(/_/g, " "), name));
  sel.onchange = () => { $("#netlist").value = ex[sel.value]; run(); };
  sel.onchange();
}

function status(msg, err) { const s = $("#status"); s.textContent = msg; s.className = err ? "err" : ""; }

async function run() {
  status("Running\u2026"); $("#stats").textContent = "";
  try {
    const t0 = performance.now();
    const res = await fetch("/api/simulate", { method: "POST", body: $("#netlist").value });
    const out = await res.json();
    if (out.error) { status("Error: " + out.error, true); return; }
    data = out; hidden.clear();
    status(`${out.analysis.toUpperCase()}: ${out.rows.length} point(s), ${out.columns.length} column(s), ${Math.round(performance.now() - t0)} ms round trip`);
    if (out.stats) $("#stats").textContent =
      `${out.stats.unknowns} unknowns \u00b7 nnz(A) ${out.stats.nnzA} \u00b7 nnz(LU) ${out.stats.nnzLU} \u00b7 ` +
      `factor ${(out.stats.factorSeconds * 1000).toFixed(1)} ms \u00b7 solve ${(out.stats.solveSeconds * 1000).toFixed(1)} ms`;
    buildSignals();
    render();
  } catch (e) { status("Request failed: " + e.message, true); }
}

function buildSignals() {
  signals = [];
  if (data.analysis === "ac") {  // columns come in pairs: vdb(x), vp(x)
    for (let i = 1; i + 1 < data.columns.length; i += 2)
      signals.push({ name: data.columns[i].replace("db(", "(") , y: i, y2: i + 1 });
  } else if (data.analysis !== "op") {
    for (let i = 1; i < data.columns.length; i++) signals.push({ name: data.columns[i], y: i });
  }
}

function render() {
  const wrap = $("#opWrap"), table = $("#opTable"), c1 = $("#plot"), c2 = $("#plot2");
  if (data.analysis === "op") {
    c1.style.display = c2.style.display = "none"; $("#legend").innerHTML = ""; $("#readout").textContent = "";
    wrap.style.display = "block";
    table.innerHTML = "<tr><th>Quantity</th><th>Value</th></tr>" +
      data.columns.map((c, i) => `<tr><td>${c}</td><td>${eng(data.rows[0][i])}${c.startsWith("i") ? "A" : "V"}</td></tr>`).join("");
    return;
  }
  wrap.style.display = "none"; c1.style.display = "block";
  c2.style.display = data.analysis === "ac" ? "block" : "none";
  $("#legend").innerHTML = signals.map((s, i) =>
    `<label><input type="checkbox" data-i="${i}" ${hidden.has(i) ? "" : "checked"}>` +
    `<i style="background:${COLORS[i % COLORS.length]}"></i>${s.name}</label>`).join("");
  document.querySelectorAll("#legend input").forEach((cb) => cb.onchange = () => {
    const i = +cb.dataset.i; cb.checked ? hidden.delete(i) : hidden.add(i); draw();
  });
  draw();
}

const M = { l: 62, r: 14, t: 12, b: 34 };

function plot(canvas, xcol, ycol, opts = {}) {
  const dpr = window.devicePixelRatio || 1;
  canvas.width = canvas.clientWidth * dpr; canvas.height = canvas.clientHeight * dpr;
  const g = canvas.getContext("2d"); g.scale(dpr, dpr);
  const W = canvas.clientWidth, H = canvas.clientHeight, css = getComputedStyle(document.body);
  const fg = css.color, grid = css.getPropertyValue("--line").trim(), muted = css.getPropertyValue("--muted").trim();
  const vis = signals.map((s, i) => [s, i]).filter(([, i]) => !hidden.has(i));
  const X = (r) => (opts.logx ? Math.log10(r[0]) : r[0]);
  const xmin = X(data.rows[0]), xmax = X(data.rows[data.rows.length - 1]);
  let ymin = Infinity, ymax = -Infinity;
  for (const r of data.rows) for (const [s] of vis) { const v = r[s[ycol]]; if (v < ymin) ymin = v; if (v > ymax) ymax = v; }
  if (!isFinite(ymin)) { ymin = 0; ymax = 1; }
  if (ymin === ymax) { ymin -= 1; ymax += 1; }
  const pad = (ymax - ymin) * 0.05; ymin -= pad; ymax += pad;
  const px = (x) => M.l + ((x - xmin) / (xmax - xmin || 1)) * (W - M.l - M.r);
  const py = (y) => H - M.b - ((y - ymin) / (ymax - ymin)) * (H - M.t - M.b);

  g.font = "11px system-ui"; g.lineWidth = 1; g.strokeStyle = grid; g.fillStyle = muted;
  g.textAlign = "right"; g.textBaseline = "middle";
  for (const t of niceTicks(ymin, ymax)) { g.beginPath(); g.moveTo(M.l, py(t)); g.lineTo(W - M.r, py(t)); g.stroke(); g.fillText(opts.ylabel ? t.toPrecision(3).replace(/\.?0+$/, "") : eng(t), M.l - 6, py(t)); }
  g.textAlign = "center"; g.textBaseline = "top";
  if (opts.logx) {  // decade gridlines
    for (let d = Math.ceil(xmin); d <= xmax; d++) { g.beginPath(); g.moveTo(px(d), M.t); g.lineTo(px(d), H - M.b); g.stroke(); g.fillText(eng(Math.pow(10, d)) + "Hz", px(d), H - M.b + 6); }
  } else {
    for (const t of niceTicks(xmin, xmax)) { g.beginPath(); g.moveTo(px(t), M.t); g.lineTo(px(t), H - M.b); g.stroke(); g.fillText(eng(t), px(t), H - M.b + 6); }
  }
  g.fillStyle = fg; g.textAlign = "left"; g.textBaseline = "top";
  g.fillText(opts.ylabel || "", M.l + 4, M.t + 2);
  if (!opts.logx) { g.textAlign = "right"; g.fillText(data.columns[0], W - M.r, H - 14); }

  for (const [s, i] of vis) {
    g.strokeStyle = COLORS[i % COLORS.length]; g.lineWidth = 1.8; g.beginPath();
    data.rows.forEach((r, k) => (k ? g.lineTo(px(X(r)), py(r[s[ycol]])) : g.moveTo(px(X(r)), py(r[s[ycol]]))));
    g.stroke();
  }
  if (opts.cursor !== undefined) {
    g.strokeStyle = muted; g.beginPath(); g.moveTo(px(opts.logx ? Math.log10(opts.cursor) : opts.cursor), M.t);
    g.lineTo(px(opts.logx ? Math.log10(opts.cursor) : opts.cursor), H - M.b); g.stroke();
  }
}

function draw(cursor) {
  if (!data || data.analysis === "op") return;
  if (data.analysis === "ac") {
    plot($("#plot"), 0, "y", { logx: true, ylabel: "magnitude (dB)", cursor });
    plot($("#plot2"), 0, "y2", { logx: true, ylabel: "phase (deg)", cursor });
  } else {
    plot($("#plot"), 0, "y", { cursor });
  }
}

function hover(e, canvas) {
  if (!data || data.analysis === "op") return;
  const rect = canvas.getBoundingClientRect();
  const f = (e.clientX - rect.left - M.l) / (rect.width - M.l - M.r);
  if (f < 0 || f > 1) return;
  const logx = data.analysis === "ac", n = data.rows.length;
  const a = logx ? Math.log10(data.rows[0][0]) : data.rows[0][0];
  const b = logx ? Math.log10(data.rows[n - 1][0]) : data.rows[n - 1][0];
  const target = a + f * (b - a);
  let lo = 0, hi = n - 1;  // nearest row by x (rows are sorted)
  while (lo < hi) { const mid = (lo + hi) >> 1; const mx = logx ? Math.log10(data.rows[mid][0]) : data.rows[mid][0]; mx < target ? (lo = mid + 1) : (hi = mid); }
  const row = data.rows[lo];
  const parts = signals.map((s, i) => hidden.has(i) ? "" :
    data.analysis === "ac" ? `${s.name}=${row[s.y].toFixed(2)}dB/${row[s.y2].toFixed(1)}\u00b0` : `${s.name}=${eng(row[s.y])}`);
  $("#readout").textContent = `${data.columns[0]}=${eng(row[0])}   ` + parts.filter(Boolean).slice(0, 8).join("   ");
  draw(row[0]);
}

for (const id of ["#plot", "#plot2"]) {
  $(id).addEventListener("mousemove", (e) => hover(e, $(id)));
  $(id).addEventListener("mouseleave", () => { $("#readout").textContent = ""; draw(); });
}
window.addEventListener("resize", () => draw());
$("#run").onclick = run;
$("#netlist").addEventListener("keydown", (e) => { if ((e.ctrlKey || e.metaKey) && e.key === "Enter") run(); });
loadExamples();

/* =============================================================================
 * app.js — VPD Dashboard interactive layer
 * Requires globals: Plotly, window.DASHBOARD_DATA, window.ReedFrost
 * Works from file:// — no fetch, no import, no bundler.
 * =========================================================================== */
(function () {
  "use strict";

  // ── Globals ────────────────────────────────────────────────────────────────
  const D  = window.DASHBOARD_DATA;
  const RF = window.ReedFrost;
  const M  = D.meta;          // N, p_default, generations, n_windows, herd, ages, counties

  // ── Palette ────────────────────────────────────────────────────────────────
  const NAVY      = "#1f4e79";
  const MIDBLUE   = "#2e75b6";
  const LIGHTBLUE = "#d6e4f0";
  const CRIMSON   = "#dc143c";
  const PLOTBG    = "#f9fafc";

  const COUNTY_COLORS = { Scruggs: NAVY, Simone: MIDBLUE, Watson: "#5ba4cf" };
  const COMP_COLORS   = { S: "#1565c0", E: "#e65100", I: "#c62828", R: "#2e7d32", cum: NAVY };

  // ── Helpers ─────────────────────────────────────────────────────────────────
  function hexToRgba(hex, a) {
    const r = parseInt(hex.slice(1, 3), 16);
    const g = parseInt(hex.slice(3, 5), 16);
    const b = parseInt(hex.slice(5, 7), 16);
    return "rgba(" + r + "," + g + "," + b + "," + a + ")";
  }

  function setEl(id, text) {
    const el = document.getElementById(id);
    if (el) el.textContent = text;
  }

  // Vertical shape at x = M.herd for histograms (frac on x-axis)
  function vHerdShape() {
    return {
      type: "line",
      x0: M.herd, x1: M.herd, xref: "x",
      y0: 0, y1: 1, yref: "paper",
      line: { color: CRIMSON, dash: "dash", width: 2 }
    };
  }

  // Dummy legend entry for the herd line (used with shape-based vertical lines)
  function herdLegendTrace() {
    return {
      x: [null], y: [null],
      type: "scatter", mode: "lines",
      name: "95% herd threshold",
      line: { color: CRIMSON, dash: "dash", width: 1.5 },
      showlegend: true, hoverinfo: "skip"
    };
  }

  // Horizontal scatter trace for continuous x-axis plots
  function herdLineTrace(x0, x1) {
    return {
      x: [x0, x1], y: [M.herd, M.herd],
      type: "scatter", mode: "lines",
      name: "95% herd threshold",
      line: { color: CRIMSON, dash: "dash", width: 1.5 },
      hoverinfo: "skip"
    };
  }

  // CI ribbon (filled area between q97_5 and q2_5)
  function ribbonTrace(ages, data, color, name) {
    const upper = data.map(d => d.q97_5);
    const lower = data.map(d => d.q2_5);
    return {
      x: ages.concat(ages.slice().reverse()),
      y: upper.concat(lower.slice().reverse()),
      fill: "toself",
      fillcolor: hexToRgba(color, 0.18),
      line: { color: "transparent" },
      name: name + " 95% CI",
      showlegend: false,
      type: "scatter",
      hoverinfo: "skip"
    };
  }

  // Base Plotly layout defaults
  const BASE = {
    paper_bgcolor: "white",
    plot_bgcolor:  PLOTBG,
    font: { family: "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif", size: 12 },
    margin: { t: 50, b: 65, l: 68, r: 25 },
    showlegend: true,
    legend: { orientation: "h", y: -0.22, x: 0, xanchor: "left", font: { size: 11 } }
  };

  function layout(extra) {
    return Object.assign({}, BASE, extra);
  }

  // ── Tab switching ────────────────────────────────────────────────────────────
  const tabInited = {};

  function switchTab(name) {
    document.querySelectorAll(".tab-btn").forEach(b => {
      b.classList.toggle("active", b.dataset.tab === name);
    });
    document.querySelectorAll(".tab-panel").forEach(p => {
      p.classList.toggle("active", p.id === "tab-" + name);
    });
    if (!tabInited[name]) {
      if (name === "eda")      initEDA();
      if (name === "coverage") initCoverage();
      if (name === "outbreak") initOutbreak();
      tabInited[name] = true;
    }
    window.dispatchEvent(new Event("resize"));
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // TAB 1 — EDA
  // ═══════════════════════════════════════════════════════════════════════════
  function initEDA() {
    renderEdaTrend();
    renderEdaHistogram();
    renderEdaBoxplot();
    renderEdaScatter();

    document.getElementById("cb-24m").addEventListener("change", renderEdaTrend);
    document.getElementById("cb-36m").addEventListener("change", renderEdaTrend);
    document.getElementById("eda-county").addEventListener("change", renderEdaHistogram);
    document.getElementById("eda-binwidth").addEventListener("input", function () {
      document.getElementById("eda-binwidth-val").textContent =
        parseFloat(this.value).toFixed(2);
      renderEdaHistogram();
    });
  }

  // ── EDA 1: State MCV1 trend ─────────────────────────────────────────────────
  function renderEdaTrend() {
    const show24 = document.getElementById("cb-24m").checked;
    const show36 = document.getElementById("cb-36m").checked;
    const child = D.eda.child;
    const ageGroups = ["24 months", "36 months"];
    const colors    = { "24 months": MIDBLUE, "36 months": NAVY };
    const shown     = { "24 months": show24, "36 months": show36 };

    const years = child.map(d => d.year);
    const minYear = Math.min.apply(null, years);
    const maxYear = Math.max.apply(null, years);

    const traces = ageGroups.map(age => {
      const rows = child.filter(d => d.age === age).slice().sort((a, b) => a.year - b.year);
      return {
        x: rows.map(d => d.year),
        y: rows.map(d => d.frac),
        type: "scatter", mode: "lines+markers",
        name: age,
        line: { color: colors[age], width: 2 },
        marker: { size: 5, color: colors[age] },
        visible: shown[age] ? true : "legendonly",
        hovertemplate: "Year: %{x}<br>Coverage: %{y:.1%}<extra>" + age + "</extra>"
      };
    });

    traces.push(herdLineTrace(minYear, maxYear));

    Plotly.react("plot-eda-trend", traces, layout({
      title: { text: "State MCV1 Coverage Trend", font: { color: NAVY, size: 14 } },
      xaxis: { title: "Year", dtick: 5 },
      yaxis: { title: "Coverage", tickformat: ".0%", range: [0.65, 1.01] },
      transition: { duration: 300 }
    }), { responsive: true });
  }

  // ── EDA 2: School coverage histogram ────────────────────────────────────────
  function renderEdaHistogram() {
    const county = document.getElementById("eda-county").value;
    const bw     = parseFloat(document.getElementById("eda-binwidth").value);
    let schools  = D.eda.schools;
    if (county !== "All") schools = schools.filter(s => s.county === county);
    const fracs = schools.map(s => s.frac);

    const traces = [
      {
        x: fracs,
        type: "histogram",
        xbins: { start: 0, end: 1.001, size: bw },
        marker: { color: hexToRgba(MIDBLUE, 0.55), line: { color: MIDBLUE, width: 0.8 } },
        name: "Schools",
        hovertemplate: "Coverage bin: %{x:.2f}<br>Count: %{y}<extra></extra>"
      },
      herdLegendTrace()
    ];

    const titleSuffix = county !== "All" ? " — " + county : "";
    Plotly.react("plot-eda-hist", traces, layout({
      title: { text: "School Two-Dose Coverage Distribution" + titleSuffix,
               font: { color: NAVY, size: 14 } },
      xaxis: { title: "Two-Dose Coverage (fraction)", tickformat: ".0%", range: [0, 1.05] },
      yaxis: { title: "Count" },
      barmode: "overlay",
      shapes: [vHerdShape()],
      transition: { duration: 300 }
    }), { responsive: true });
  }

  // ── EDA 3: County box plot ───────────────────────────────────────────────────
  function renderEdaBoxplot() {
    const schools  = D.eda.schools;
    const counties = M.counties;

    const traces = counties.map(county => ({
      y: schools.filter(s => s.county === county).map(s => s.frac),
      type: "box",
      name: county,
      boxmean: true,
      marker: { color: COUNTY_COLORS[county] || MIDBLUE, opacity: 0.75 },
      line:   { color: COUNTY_COLORS[county] || MIDBLUE },
      hovertemplate: county + "<br>Coverage: %{y:.1%}<extra></extra>"
    }));

    // Herd line across all county categories
    traces.push({
      x: counties,
      y: counties.map(() => M.herd),
      type: "scatter", mode: "lines",
      name: "95% herd threshold",
      line: { color: CRIMSON, dash: "dash", width: 1.5 },
      hoverinfo: "skip"
    });

    Plotly.react("plot-eda-box", traces, layout({
      title: { text: "County Two-Dose Coverage Distribution", font: { color: NAVY, size: 14 } },
      xaxis: { title: "County" },
      yaxis: { title: "Two-Dose Coverage", tickformat: ".0%", range: [0.35, 1.05] }
    }), { responsive: true });
  }

  // ── EDA 4: Enrollment vs coverage scatter ───────────────────────────────────
  function renderEdaScatter() {
    const schools   = D.eda.schools;
    const counties  = M.counties;
    const maxEnroll = Math.max.apply(null, schools.map(s => s.enrollment));

    const traces = counties.map(county => {
      const rows = schools.filter(s => s.county === county);
      return {
        x: rows.map(s => s.enrollment),
        y: rows.map(s => s.frac),
        text: rows.map(s => s.school),
        type: "scatter", mode: "markers",
        name: county,
        marker: { color: COUNTY_COLORS[county] || MIDBLUE, size: 7, opacity: 0.75,
                  line: { color: "white", width: 0.5 } },
        hovertemplate: "%{text}<br>Enrollment: %{x}<br>Coverage: %{y:.1%}<extra>"
                       + county + "</extra>"
      };
    });

    traces.push(herdLineTrace(0, maxEnroll + 5));

    Plotly.react("plot-eda-scatter", traces, layout({
      title: { text: "Enrollment vs. Two-Dose Coverage", font: { color: NAVY, size: 14 } },
      xaxis: { title: "Enrollment" },
      yaxis: { title: "Coverage", tickformat: ".0%", range: [0.3, 1.05] }
    }), { responsive: true });
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // TAB 2 — COVERAGE ESTIMATES
  // ═══════════════════════════════════════════════════════════════════════════
  let covScenario = "baseline";
  let covLevel    = "state";

  function initCoverage() {
    document.querySelectorAll("[data-scenario]").forEach(btn => {
      btn.addEventListener("click", function () {
        document.querySelectorAll("[data-scenario]").forEach(b => b.classList.remove("active"));
        this.classList.add("active");
        covScenario = this.dataset.scenario;
        renderCoverage();
      });
    });
    document.querySelectorAll("[data-level]").forEach(btn => {
      btn.addEventListener("click", function () {
        document.querySelectorAll("[data-level]").forEach(b => b.classList.remove("active"));
        this.classList.add("active");
        covLevel = this.dataset.level;
        renderCoverage();
      });
    });
    renderCoverage();
  }

  function renderCoverage() {
    const isReduced = covScenario === "reduced";
    const scen      = isReduced ? "reduced" : "baseline";
    const traces    = [];

    if      (covLevel === "state")  buildCovState(traces, scen, isReduced);
    else if (covLevel === "county") buildCovCounty(traces, scen, isReduced);
    else                            buildCovSchool(traces, scen, isReduced);

    const levelLabel = covLevel === "state" ? "State"
                      : covLevel === "county" ? "County"
                      : "School (Scruggs)";
    const scenLabel  = isReduced ? "Reduced" : "Baseline";
    const yRange     = covLevel === "school" ? [0.5, 1.02] : [0.72, 1.02];

    Plotly.react("plot-coverage", traces, layout({
      title: { text: levelLabel + " Two-Dose Coverage by Age — " + scenLabel,
               font: { color: NAVY, size: 14 } },
      xaxis: { title: "Age (years)", dtick: 1 },
      yaxis: { title: "Coverage", tickformat: ".0%", range: yRange },
      margin: { t: 50, b: 70, l: 68, r: 25 },
      transition: { duration: 400, easing: "cubic-in-out" }
    }), { responsive: true });

    updateCovKPIs(scen);
  }

  function buildCovState(traces, scen, isReduced) {
    const data = D.coverage.state[scen];
    const ages = data.map(d => d.age);

    traces.push(ribbonTrace(ages, data, MIDBLUE, "State"));
    traces.push({
      x: ages, y: data.map(d => d.q50),
      type: "scatter", mode: "lines+markers",
      name: isReduced ? "Reduced" : "Baseline",
      line: { color: MIDBLUE, width: 2.5 },
      marker: { size: 6, color: MIDBLUE },
      hovertemplate: "Age %{x}<br>Coverage: %{y:.1%}<extra></extra>"
    });

    if (isReduced) {
      const base = D.coverage.state.baseline;
      traces.push({
        x: base.map(d => d.age), y: base.map(d => d.q50),
        type: "scatter", mode: "lines",
        name: "Baseline (ghost)",
        line: { color: NAVY, dash: "dot", width: 1.5 },
        opacity: 0.4,
        hovertemplate: "Age %{x}<br>Baseline: %{y:.1%}<extra></extra>"
      });
    }

    traces.push(herdLineTrace(ages[0], ages[ages.length - 1]));
  }

  function buildCovCounty(traces, scen, isReduced) {
    M.counties.forEach(county => {
      const data = D.coverage.county[scen][county];
      const ages = data.map(d => d.age);
      const col  = COUNTY_COLORS[county] || MIDBLUE;

      traces.push(ribbonTrace(ages, data, col, county));
      traces.push({
        x: ages, y: data.map(d => d.q50),
        type: "scatter", mode: "lines+markers",
        name: county,
        line: { color: col, width: 2 },
        marker: { size: 5, color: col },
        hovertemplate: county + " Age %{x}<br>Coverage: %{y:.1%}<extra></extra>"
      });

      if (isReduced) {
        const base = D.coverage.county.baseline[county];
        traces.push({
          x: base.map(d => d.age), y: base.map(d => d.q50),
          type: "scatter", mode: "lines",
          name: county + " baseline",
          showlegend: false,
          line: { color: col, dash: "dot", width: 1 },
          opacity: 0.35,
          hoverinfo: "skip"
        });
      }
    });

    const ages = D.coverage.county.baseline[M.counties[0]].map(d => d.age);
    traces.push(herdLineTrace(ages[0], ages[ages.length - 1]));
  }

  function buildCovSchool(traces, scen, isReduced) {
    const data = D.coverage.school[scen];
    const ages = M.ages;

    ages.forEach(age => {
      const vals = data.filter(d => d.age === age).map(d => d.q50);
      traces.push({
        y: vals,
        name: String(age),
        type: "box", boxmean: true,
        marker: { color: hexToRgba(MIDBLUE, 0.65) },
        line:   { color: NAVY },
        showlegend: false,
        hovertemplate: "Age " + age + "<br>Coverage: %{y:.1%}<extra></extra>"
      });
    });

    if (isReduced) {
      const base  = D.coverage.school.baseline;
      const ghostY = ages.map(age => {
        const vals = base.filter(d => d.age === age).map(d => d.q50);
        return vals.reduce((a, b) => a + b, 0) / vals.length;
      });
      traces.push({
        x: ages.map(String), y: ghostY,
        type: "scatter", mode: "lines+markers",
        name: "Baseline median (ghost)",
        line: { color: NAVY, dash: "dot", width: 1.5 },
        opacity: 0.45,
        hovertemplate: "Age %{x}<br>Baseline: %{y:.1%}<extra></extra>"
      });
    }

    // Herd line across categorical age axis
    traces.push({
      x: ages.map(String), y: ages.map(() => M.herd),
      type: "scatter", mode: "lines",
      name: "95% herd threshold",
      line: { color: CRIMSON, dash: "dash", width: 1.5 },
      hoverinfo: "skip"
    });
  }

  function updateCovKPIs(scen) {
    const stateData = D.coverage.state[scen];
    const row18  = stateData.find(d => d.age === 18);
    const cov18  = row18 ? row18.q50 : null;
    const base18 = D.coverage.state.baseline.find(d => d.age === 18);
    const b18    = base18 ? base18.q50 : null;
    const gap    = (scen === "reduced" && b18 !== null && cov18 !== null)
                   ? ((b18 - cov18) * 100).toFixed(1) + " pp"
                   : "—";
    const nBelow = stateData.filter(d => d.q50 < M.herd).length;

    setEl("kpi-cov18",   cov18 !== null ? (cov18 * 100).toFixed(1) + "%" : "—");
    setEl("kpi-gap",     gap);
    setEl("kpi-below95", nBelow + " / " + stateData.length + " ages");
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // TAB 3 — OUTBREAK DYNAMICS
  // ═══════════════════════════════════════════════════════════════════════════
  let obSeed       = 42;
  let playInterval = null;
  let obDebounce   = null;

  function initOutbreak() {
    const winSlider = document.getElementById("ob-window");
    const pSlider   = document.getElementById("ob-p");
    winSlider.max   = M.n_windows;
    winSlider.min   = 1;
    winSlider.value = 1;
    pSlider.value   = M.p_default;
    document.getElementById("ob-p-val").textContent      = M.p_default.toFixed(2);
    document.getElementById("ob-window-val").textContent = "1";

    winSlider.addEventListener("input", function () {
      document.getElementById("ob-window-val").textContent = this.value;
      scheduleObUpdate();
    });
    pSlider.addEventListener("input", function () {
      document.getElementById("ob-p-val").textContent = parseFloat(this.value).toFixed(2);
      scheduleObUpdate();
    });
    document.getElementById("ob-sims").addEventListener("change", scheduleObUpdate);
    document.getElementById("ob-play").addEventListener("click", togglePlay);
    document.getElementById("ob-reroll").addEventListener("click", function () {
      obSeed = Math.floor(Math.random() * 0xFFFFFF);
      runObUpdate();
    });

    runObUpdate();
  }

  function scheduleObUpdate() {
    if (obDebounce) clearTimeout(obDebounce);
    obDebounce = setTimeout(runObUpdate, 150);
  }

  function togglePlay() {
    const btn = document.getElementById("ob-play");
    if (playInterval) {
      clearInterval(playInterval);
      playInterval = null;
      btn.innerHTML = "&#9654; Play";
    } else {
      btn.innerHTML = "&#9646;&#9646; Pause";
      playInterval = setInterval(function () {
        const sl  = document.getElementById("ob-window");
        const cur = parseInt(sl.value);
        const nxt = cur >= M.n_windows ? 1 : cur + 1;
        sl.value  = nxt;
        document.getElementById("ob-window-val").textContent = nxt;
        runObUpdate();
      }, 700);
    }
  }

  function runObUpdate() {
    const winIdx  = parseInt(document.getElementById("ob-window").value) - 1;
    const p       = parseFloat(document.getElementById("ob-p").value);
    const nSims   = parseInt(document.getElementById("ob-sims").value);
    const winData = D.windows[winIdx];

    const result = RF.runWindow({
      generations: M.generations,
      S0:          winData.S0,
      R0seed:      winData.R0_seed,
      p:           p,
      nSims:       nSims,
      seed:        obSeed,
      N:           M.N
    });

    renderTrajectories(result.bands, M.generations);
    renderOutbreakHist(result.finals, result.stats);
    updateObKPIs(winData, result.stats, winIdx + 1);
  }

  // ── Tab 3 Fig A: Trajectories (5-panel small multiples) ────────────────────
  function renderTrajectories(bands, gens) {
    const x = [];
    for (let g = 0; g <= gens; g++) x.push(g);

    const comps = [
      { key: "S",   label: "Susceptible",        color: COMP_COLORS.S,   row: 0, col: 0 },
      { key: "E",   label: "Exposed",             color: COMP_COLORS.E,   row: 0, col: 1 },
      { key: "I",   label: "Infectious",          color: COMP_COLORS.I,   row: 1, col: 0 },
      { key: "R",   label: "Recovered",           color: COMP_COLORS.R,   row: 1, col: 1 },
      { key: "cum", label: "Cumulative Infected", color: COMP_COLORS.cum, row: 2, col: 0 }
    ];

    // Domain layout: 3 rows × 2 cols
    const xDom = [[0, 0.46], [0.54, 1.0]];
    const yDom = [[0.69, 0.97], [0.36, 0.64], [0.03, 0.31]];

    const traces      = [];
    const trajLayout  = {
      paper_bgcolor: "white",
      plot_bgcolor:  PLOTBG,
      font: { family: "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif",
              size: 11 },
      margin:     { t: 15, b: 35, l: 58, r: 15 },
      showlegend: false,
      height:     490,
      annotations: []
    };

    comps.forEach((comp, idx) => {
      const r = comp.row;
      const c = comp.col;

      // Axis naming: "" for idx 0, "2" for idx 1, etc.
      const axN = idx === 0 ? "" : String(idx + 1);
      const xKey = "xaxis" + axN;
      const yKey = "yaxis" + axN;
      const xRef = "x" + axN;
      const yRef = "y" + axN;

      const bnd = bands[comp.key];
      const med = bnd.map(b => b.median);
      const q25 = bnd.map(b => b.q25);
      const q75 = bnd.map(b => b.q75);

      // IQR ribbon
      traces.push({
        x: x.concat(x.slice().reverse()),
        y: q75.concat(q25.slice().reverse()),
        fill: "toself",
        fillcolor: hexToRgba(comp.color, 0.18),
        line: { color: "transparent" },
        showlegend: false,
        xaxis: xRef, yaxis: yRef,
        hoverinfo: "skip",
        type: "scatter"
      });

      // Median line
      traces.push({
        x: x, y: med,
        type: "scatter", mode: "lines",
        name: comp.label,
        line: { color: comp.color, width: 2.2 },
        xaxis: xRef, yaxis: yRef,
        hovertemplate: "Gen %{x}: %{y:.1f}<extra>" + comp.label + "</extra>"
      });

      // Show x-axis title only on bottom panels of each column
      const isBottomOfCol = (r === 2 && c === 0) || (r === 1 && c === 1);
      trajLayout[xKey] = {
        domain:    xDom[c],
        tickfont:  { size: 10 },
        title:     isBottomOfCol ? { text: "Generation", font: { size: 10 } } : ""
      };
      trajLayout[yKey] = {
        domain:    yDom[r],
        tickfont:  { size: 10 },
        title:     { text: comp.label, font: { size: 10.5, color: comp.color } }
      };

      // Panel label annotation (top-left of panel)
      trajLayout.annotations.push({
        xref: "paper", yref: "paper",
        x: xDom[c][0] + 0.015,
        y: yDom[r][1] - 0.005,
        text: "<b>" + comp.label + "</b>",
        showarrow: false,
        xanchor: "left", yanchor: "top",
        font: { size: 10, color: comp.color }
      });
    });

    Plotly.react("plot-ob-traj", traces, trajLayout, { responsive: true });
  }

  // ── Tab 3 Fig B: Outbreak histogram ─────────────────────────────────────────
  function renderOutbreakHist(finals, stats) {
    const hist    = RF.histogram(finals, 1);
    const maxY    = Math.max.apply(null, hist.y.concat([1]));
    const median  = stats.median;
    const nExtinct = finals.filter(v => v < 5).length;
    const firstBinY = hist.y.length > 0 ? hist.y[0] : 0;

    const traces = [
      {
        x: hist.x, y: hist.y,
        type: "bar", width: 1,
        marker: { color: hexToRgba(MIDBLUE, 0.55), line: { color: MIDBLUE, width: 0.5 } },
        name: "Simulations",
        hovertemplate: "Infections: %{x:.0f}<br>Count: %{y}<extra></extra>"
      },
      {
        x: [median, median], y: [0, maxY * 1.08],
        type: "scatter", mode: "lines",
        name: "Median (" + median.toFixed(0) + ")",
        line: { color: NAVY, dash: "dash", width: 2 },
        hovertemplate: "Median: " + median.toFixed(1) + "<extra></extra>"
      }
    ];

    // Extinction spike annotation
    const annotations = [{
      x: hist.x.length > 0 ? hist.x[0] : 0.5,
      y: firstBinY,
      xanchor: "left", yanchor: "bottom",
      text: "Extinction<br>spike (" + nExtinct + ")",
      showarrow: true, arrowhead: 2, arrowcolor: CRIMSON,
      ax: 35, ay: -40,
      font: { color: CRIMSON, size: 10 }
    }];

    Plotly.react("plot-ob-hist", traces, layout({
      title: { text: "Outbreak Size Distribution", font: { color: NAVY, size: 14 } },
      xaxis: { title: "Total Infections (cumulative)", range: [0, 70] },
      yaxis: { title: "Count" },
      barmode: "overlay",
      annotations: annotations,
      transition: { duration: 200 }
    }), { responsive: true });
  }

  // ── Tab 3 KPIs ───────────────────────────────────────────────────────────────
  function updateObKPIs(winData, stats, windowNum) {
    setEl("kpi-ob-window",    String(windowNum));
    setEl("kpi-r0seed",       String(winData.R0_seed));
    setEl("kpi-median",       stats.median.toFixed(1));
    setEl("kpi-median-major", stats.medianMajor.toFixed(1));
    setEl("kpi-pmajor",       (stats.pMajor * 100).toFixed(1) + "%");
    setEl("kpi-pextinct",     (stats.pExtinct * 100).toFixed(1) + "%");
    setEl("kpi-pct-infected", stats.pctInfected.toFixed(1) + "%");
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // BOOT
  // ═══════════════════════════════════════════════════════════════════════════
  document.addEventListener("DOMContentLoaded", function () {
    document.querySelectorAll(".tab-btn").forEach(btn => {
      btn.addEventListener("click", function () {
        switchTab(this.dataset.tab);
      });
    });
    switchTab("eda");
  });

})();

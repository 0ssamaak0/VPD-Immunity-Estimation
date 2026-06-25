/* =============================================================================
 * reedfrost.js — client-side Reed-Frost SEIR engine.
 *
 * Exact port of simulate_reed_frost_seir() in
 * Phase_2_Dynamical_part/Reed_Frost_model.R, plus a seeded RNG so runs are
 * reproducible. Runs live in the browser (and under Node for tests).
 *
 *   S_{t+1} = S_t - newE
 *   E_{t+1} = newE ~ Binomial(S_t, 1 - (1-p)^{I_t})   [0 if I_t == 0]
 *   I_{t+1} = E_t
 *   R_{t+1} = R_t + I_t
 *   cumulative = (E + I + R) - R(0)                    [matches histogram.R]
 * ========================================================================== */
(function (root) {
  "use strict";

  // --- Seeded RNG (mulberry32): deterministic given a seed ------------------
  function makeRng(seed) {
    let a = seed >>> 0;
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // Binomial draw by summing Bernoullis (n = S is small, <~ 70). --------------
  function rbinom(n, p, rng) {
    if (p <= 0 || n <= 0) return 0;
    if (p >= 1) return n;
    let k = 0;
    for (let i = 0; i < n; i++) if (rng() < p) k++;
    return k;
  }

  // One stochastic trajectory. y0 = {S, E, I, R}. Returns array per generation.
  function simulate(generations, y0, p, rng) {
    let S = y0.S, E = y0.E, I = y0.I, R = y0.R;
    const R_init = R;
    const out = [{ gen: 0, S: S, E: E, I: I, R: R, cum: 0 }];
    for (let t = 1; t <= generations; t++) {
      const newE = I === 0 ? 0 : rbinom(S, 1 - Math.pow(1 - p, I), rng);
      const nS = S - newE, nE = newE, nI = E, nR = R + I;
      S = nS; E = nE; I = nI; R = nR;
      out.push({ gen: t, S: S, E: E, I: I, R: R, cum: (E + I + R) - R_init });
    }
    return out;
  }

  // Run `nSims` trajectories for one window. Returns:
  //   finals[]          : final cumulative infection per sim
  //   bands             : per-generation {median,q25,q75} for S,E,I,R,cum
  //   stats             : median, mean, pMajor (>=5), pExtinct (<5), pctInfected
  function runWindow(opts) {
    const { generations, S0, R0seed, p, nSims, seed, majorThreshold = 5 } = opts;
    const rng = makeRng(seed >>> 0);
    const G = generations;
    const comps = ["S", "E", "I", "R", "cum"];

    // collectors[comp][gen] = array over sims
    const collectors = {};
    comps.forEach(c => { collectors[c] = Array.from({ length: G + 1 }, () => []); });
    const finals = new Array(nSims);

    for (let s = 0; s < nSims; s++) {
      const traj = simulate(G, { S: S0, E: 1, I: 0, R: R0seed }, p, rng);
      for (let g = 0; g <= G; g++) {
        const row = traj[g];
        for (const c of comps) collectors[c][g].push(row[c]);
      }
      finals[s] = traj[G].cum;
    }

    const q = (arr, prob) => {
      const a = arr.slice().sort((x, y) => x - y);
      const idx = prob * (a.length - 1);
      const lo = Math.floor(idx), hi = Math.ceil(idx);
      return lo === hi ? a[lo] : a[lo] + (a[hi] - a[lo]) * (idx - lo);
    };

    const bands = {};
    comps.forEach(c => {
      bands[c] = collectors[c].map(arr => ({
        median: q(arr, 0.5), q25: q(arr, 0.25), q75: q(arr, 0.75)
      }));
    });

    const major = finals.filter(v => v >= majorThreshold);
    const stats = {
      R0seed: R0seed,
      median: q(finals, 0.5),
      mean: finals.reduce((a, b) => a + b, 0) / finals.length,
      medianMajor: major.length ? q(major, 0.5) : 0,
      pMajor: major.length / finals.length,
      pExtinct: 1 - major.length / finals.length,
      pctInfected: 100 * (q(finals, 0.5) / opts.N)
    };
    return { finals, bands, stats };
  }

  // Bin a finals[] array into {x[],y[]} counts for a histogram.
  function histogram(finals, binWidth) {
    const w = binWidth || 1;
    const max = Math.max(0, ...finals);
    const nBins = Math.floor(max / w) + 1;
    const y = new Array(nBins).fill(0);
    for (const v of finals) y[Math.floor(v / w)]++;
    const x = y.map((_, i) => i * w + w / 2);
    return { x: x, y: y };
  }

  const api = { makeRng, rbinom, simulate, runWindow, histogram };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.ReedFrost = api;
})(typeof window !== "undefined" ? window : globalThis);

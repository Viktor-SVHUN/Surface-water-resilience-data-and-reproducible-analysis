# -*- coding: utf-8 -*-
"""Surface-water resilience: annual physicochemical monitoring and network analysis.

Run from any working directory:
    python src/surface_water_resilience.py
    python src/surface_water_resilience.py --self-test-only
    python src/surface_water_resilience.py --output-dir path/to/new_results

Default inputs are in data/raw; results are written to a new outputs directory.
See README.md and docs/METHODS.md for data units, estimators and limitations.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import logging
import math
import platform
import random
import re
import sys
import time
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import linalg, stats
import networkx as nx
import igraph as ig
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = PROJECT_ROOT / "data" / "raw"
OUTPUT_DIR = PROJECT_ROOT / "outputs"

@dataclass
class Config:
    start_year: int = 2000
    end_year: int = 2022
    seed: int = 42
    min_years_local: int = 8
    min_years_hvg: int = 8
    min_parameters_ranking: int = 5
    critical_quantile: float = 0.90  # descriptive upper tail; not a legal cutoff
    redundancy_rho: float = 0.90
    correlation_min_n: int = 10
    correlation_min_basins: int = 8
    bootstrap_reps: int = 299
    permutations_small_n: int = 999
    sensitivity_reps: int = 20
    sensitivity_series_per_parameter: int = 50
    sensitivity_min_years: int = 12
    deletion_fractions: tuple = (0.10, 0.25, 0.50)
    # Direction-aware threshold perturbation: -0.05 = 5% stricter, +0.05 = 5% more lenient
    threshold_shifts: tuple = (-0.05, 0.05)
    permutations: int = 999
    top_associations: int = 20
    min_group_size_typology: int = 10
    min_upstream_exposed_years: int = 5
    tsne_perplexity: float = 30.0
    tsne_iterations: int = 1000
    tsne_learning_rate: float = 200.0
    attack_reps: int = 20
    run_tsne: bool = True
    run_sensitivity: bool = True
    selected_vor: str = "AIH049"

CFG = Config()
QUICK_RUN = False  # set by --quick; written into latex/analysis_values.tex
LOG = logging.getLogger("surface_water_resilience")
MASS_PARAMETERS = {
    "Ammonium (NH4-N)", "BOD (5)", "COD (Cr)", "Dissolved oxygen",
    "Nitrate (NO3-N)", "Nitrite (NO2-N)", "Orthophosphate (PO4-P)",
    "Phytoplankton (biomass  - chlorophyll-a)", "Total nitrogen", "Total phosphorus",
}
PARAM_RENAME = {"Zavarosság nefelometriás": "Turbidity (NTU)"}
LOCAL_METRICS = ["N_good", "N_bad", "Amax_bad", "Amax_good",
                 "Lmax_bad_years", "Lmax_good_years", "ARL_bad_years",
                 "ARL_good_years", "tR_first_observed_years"]

# "upstream pressure" descriptor and is declared before outcome-blind pruning.
HYDRO_PRIORITY = ["DCI", "UPS", "BC", "DCO", "DOWN", "HC_IN", "PRC", "AC_alpha09"]
HVG_PRIORITY = ["DZO", "LeM_U", "Tra_U", "GLE_D", "AVPL_D", "VAs_D"]
HVG_FULL = ["DZI", "DZO", "DZ_U", "BZ_D", "PRZ_D", "Dens_D", "Arcs",
            "Diam_D", "AVPL_D", "ALE_U", "GLE_D", "LeM_U", "IM_D",
            "Mot_D", "Tra_U", "RRes_U", "SRes_U", "VAs_D"]
HVG_FAST = ["DZO", "Dens_D", "AVPL_D", "GLE_D", "Tra_U", "VAs_D"]


def write_csv(df, path):
    """Excel-friendly UTF-8; unavailable values stay blank, never zero-filled."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(text)).strip("_")


def input_path(root, filename):
    p = Path(root) / filename
    if not p.is_file():
        raise FileNotFoundError(f"Required input missing: {p}")
    return p


def excel_read(path):
    try:
        import python_calamine  # noqa: F401
        return pd.read_excel(path, engine="calamine")
    except ImportError:
        return pd.read_excel(path, engine="openpyxl")


def discover_input(explicit=None):
    """Resolve a complete input directory independently of the working directory."""
    folder = Path(explicit or INPUT_DIR).expanduser().resolve()
    required = ("river_registry.csv", "lake_registry.csv", "water_quality_annual_site_summaries.xlsx", "operational_thresholds.xlsx")
    missing = [name for name in required if not (folder / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing inputs in {folder}: {', '.join(missing)}")
    return folder


def metadata_and_graph(rivers, lakes, raw):
    """Use named columns, sorted VOR identifiers, valid recipient links only."""
    rows, edges, rejected = [], [], []
    for source, typ in [(rivers, "River"), (lakes, "Lake")]:
        if source.VT_VOR_2015.duplicated().any():
            raise ValueError(f"Duplicate network VOR in {typ}; resolve before analysis.")
        for _, row in source.iterrows():
            v = str(row.VT_VOR_2015).strip()
            actual_type = "Lake" if typ == "River" and row.get("WB_CAT") == "LW" else typ
            rows.append({
                "VOR": v, "Type": actual_type, "name": row.get("VT_NEV_2015"),
                "metadata_catalogue": typ,
                "basin": ("Unknown" if pd.isna(row.get("ALEGYSEG")) else str(row.get("ALEGYSEG")).strip()),
                "catchment_km2": row.get("KOZV_VGY_MERET"),
                "total_catchment_km2": row.get("TELJES_VGY_MERET"),
                "river_length_km": row.get("HOSSZ_MERT", np.nan),
                "lake_area_km2": row.get("SUM_AREA_KM2", np.nan),
                "depth_m": row.get("ATLMELYS", np.nan),
                "artificial": {"Y": 1.0, "N": 0.0}.get(row.get("ARTIFICIAL"), np.nan),
                "modified": {"Y": 1.0, "N": 0.0}.get(row.get("MODIFIED"), np.nan),
                "altitude_class": row.get("ALT_CAT", "Unknown"),
                "slope_class": row.get("SLOPE_CAT", "Not applicable"),
                "substrate_class": row.get("GRAIN_CAT", "Not applicable"),
                "geochemistry": row.get("GEOL_CAT", row.get("GEOL_CHEM_CAT", "Unknown")),
                "typology_code": str(row.get("TIPUS_UJ", row.get("TIPUS", "Unknown"))),
            })
            target = row.get("BEF_VT_VOR")
            if pd.isna(target) or str(target).strip().lower() in {"", "nan", "null"}:
                rejected.append({"from": v, "to": target, "reason": "missing_recipient"})
            elif str(target).strip() == v:
                rejected.append({"from": v, "to": target, "reason": "self_loop"})
            else:
                edges.append((v, str(target).strip()))
    meta = pd.DataFrame(rows)
    if meta.VOR.duplicated().any():
        overlap = set(rivers.VT_VOR_2015) & set(lakes.VT_VOR_2015)
        for v in overlap:
            r = rivers.loc[rivers.VT_VOR_2015.eq(v)].iloc[0]
            l = lakes.loc[lakes.VT_VOR_2015.eq(v)].iloc[0]
            if r.get("WB_CAT") != "LW" or r.get("BEF_VT_VOR") != l.get("BEF_VT_VOR"):
                raise ValueError(f"Unresolved conflicting river/lake metadata for {v}.")
        # The Rivers catalogue contains LW rows also in Lakes. Their explicit
        # WB_CAT resolves the type; use the Lakes record with actual lake fields.
        meta = meta.sort_values("metadata_catalogue").drop_duplicates("VOR",keep="first")
    registered = set(meta.VOR)
    observed = set(raw.VOR.dropna().astype(str).str.strip())
    coords = raw.assign(VOR=raw.VOR.astype("string").str.strip()).groupby("VOR").agg(
        Latitude=("FOLDRAJZI_SZELESSEG", "median"),
        Longitude=("FOLDRAJZI_HOSSZUSAG", "median"),
        monitoring_sites=("KTJ", "nunique"))
    named = raw.dropna(subset=["VOR"]).groupby("VOR").VIZTEST_MEGNEVEZES.first()
    all_vor = sorted(registered | observed | {t for _, t in edges})
    meta = meta.set_index("VOR").reindex(all_vor)
    meta["Type"] = meta.Type.fillna("Other")
    meta["name"] = meta.name.fillna(named).fillna(pd.Series(all_vor, index=all_vor))
    meta["basin"] = meta.basin.fillna("Unknown")
    meta = meta.join(coords)
    meta["registered"] = meta.index.isin(registered)
    meta["observed"] = meta.index.isin(observed)
    meta["ID"] = np.arange(1, len(meta) + 1)  # stable sorted-code index

    meta["catalogue_priority_type"] = np.select([meta.index.isin(set(rivers.VT_VOR_2015)),
                                      meta.index.isin(set(lakes.VT_VOR_2015))],
                                     ["River","Lake"],default="Other")
    graph = nx.DiGraph()
    graph.add_nodes_from(all_vor)
    graph.add_edges_from(sorted(set(edges)))
    components = sorted(nx.weakly_connected_components(graph), key=lambda c: (-len(c), min(c)))
    component_map = {v: i for i, component in enumerate(components) for v in component}
    meta["weak_component"] = pd.Series(component_map)
    return meta.reset_index(names="VOR"), graph, pd.DataFrame(rejected)


def normalize_and_aggregate(raw, meta, limits):
    """One record/site/parameter/year; equal site weights within a water body.

    Copies that differ only by reporting organization are not independent.
    Conflicting versions are collapsed symmetrically within site/year; audited.
    MERESSZAM-weighted aggregation is a sensitivity, not the primary estimator.
    """
    x = raw.copy()
    x["VOR"] = x.VOR.astype("string").str.strip().replace({"": pd.NA, "null": pd.NA})
    x["parameter"] = x.ANYAG_NEV_ANGOL.replace(PARAM_RENAME)
    x["year"] = pd.to_numeric(x.DATUM_YEAR, errors="coerce")
    x["value"] = pd.to_numeric(x.ATLAG, errors="coerce")
    x["samples_reported"] = pd.to_numeric(x.MERESSZAM, errors="coerce")
    units = x.MERTEKEGYSEG_NEV.astype(str).str.strip().str.lower()
    x["unit"] = units.map({"gramm/liter": "g/L", "g/l": "g/L", "mg/l": "g/L",
                           "százalék": "%", "mikrosiemens/centiméter": "uS/cm", "ntu": "NTU"})
    mass = x.parameter.isin(MASS_PARAMETERS)
    x.loc[mass & units.eq("mg/l"), "value"] *= 0.001
    selected = mass | x.parameter.isin(["Conductivity", "Oxygen saturation", "Turbidity (NTU)"])
    unit_ok = ((mass & units.isin(["gramm/liter", "g/l", "mg/l"])) |
               (x.parameter.eq("Conductivity") & units.eq("mikrosiemens/centiméter")) |
               (x.parameter.eq("Oxygen saturation") & units.eq("százalék")) |
               (x.parameter.eq("Turbidity (NTU)") & units.eq("ntu")))
    x["exclusion_reason"] = np.select([
        ~selected, ~unit_ok, x.VOR.isna(), ~np.isfinite(x.value),
        ~x.year.between(CFG.start_year, CFG.end_year), x.KTJ.isna(), x.value.lt(0)],
        ["outside_physicochemical_scope", "incompatible_or_unknown_unit", "missing_VOR",
         "invalid_value", "outside_period_or_invalid_year", "missing_site_id", "negative_value"], default="")
    exclusions = x[x.exclusion_reason.ne("")].groupby(
        ["exclusion_reason", "parameter", "MERTEKEGYSEG_NEV"], dropna=False).size().reset_index(name="rows")
    x = x[x.exclusion_reason.eq("")].copy()
    x["year"] = x.year.astype(int)
    x["site"] = x.KTJ.astype(str)
    keep = ["VOR", "site", "parameter", "year", "value", "unit", "samples_reported"]
    versions = x[keep].drop_duplicates()  # source-org copies removed, site retained
    keys = ["VOR", "site", "parameter", "year", "unit"]
    site = versions.groupby(keys, as_index=False).agg(
        value=("value", "mean"), version_count=("value", "size"),
        conflicting_values=("value", "nunique"), samples_reported=("samples_reported", "max"))
    # Count is an approximate coverage annotation, never claimed as unique samples.
    site["reported_sample_count_approx"] = site.samples_reported
    site["weighted_value"] = site.value * site.samples_reported.clip(lower=1)
    site["weight"] = site.samples_reported.clip(lower=1)
    annual = site.groupby(["VOR", "parameter", "year", "unit"], as_index=False).agg(
        value=("value", "mean"), station_sd=("value", "std"),
        n_sites=("site", "nunique"), approximate_sample_count=("samples_reported", "sum"),
        weighted_sum=("weighted_value", "sum"), weight_sum=("weight", "sum"),
        conflicting_site_versions=("conflicting_values", lambda s: int((s > 1).sum())))
    annual["value_sample_weighted"] = annual.weighted_sum / annual.weight_sum
    annual = annual.drop(columns=["weighted_sum", "weight_sum"])
    annual = annual.merge(meta, on="VOR", how="left", validate="many_to_one")
    if limits.indicator.duplicated().any():
        raise ValueError("Duplicate indicator in operational_thresholds.xlsx.")
    annual = annual.merge(limits, left_on="parameter", right_on="indicator", how="left", validate="many_to_one")
    annual["threshold"] = np.select([annual.Type.eq("Lake"), annual.Type.eq("River")],
                                    [annual.Lake, annual.River], default=np.nan)
    annual["direction"] = annual.Boundary.map({"lower": 1.0, "higher": -1.0})
    annual.loc[annual.parameter.eq("Turbidity (NTU)"), "direction"] = 1.0
    # Threshold unit is assumed from the supplied workbook; not certified law.
    annual["threshold_source"] = np.where(annual.threshold.notna(), "supplied_type_level_operational", "unavailable")
    annual["adverse_deviation"] = annual.direction * (annual.value - annual.threshold)
    annual["eligible_local"] = annual.threshold.gt(0) & annual.direction.notna() & annual.Type.isin(["Lake", "River"])
    annual = annual.sort_values(["VOR", "parameter", "year"]).reset_index(drop=True)
    audit = {"raw_rows": len(raw), "selected_rows_before_dedup": len(x),
             "site_version_rows_after_copy_removal": len(versions), "site_parameter_year_rows": len(site),
             "conflicting_site_parameter_years": int((site.conflicting_values > 1).sum()),
             "annual_body_parameter_rows": len(annual),
             "annual_body_parameter_combinations": annual.groupby(["VOR", "parameter"]).ngroups,
             "aggregation": "equal weight to site/year summaries; weighted variant exported",
             "threshold_policy": "no mean substitution, no River threshold for Other",
             "turbidity": "NTU only; g/L turbidity records excluded as incompatible units"}
    return annual, site, exclusions, audit


def mask_runs(mask, years, gap=1, only_true=True):
    """Count runs, splitting at unobserved calendar years. No interpolation."""
    mask, years = np.asarray(mask, bool), np.asarray(years)
    if not len(mask):
        return []
    runs, size, previous = [], 0, None
    for i, state in enumerate(mask):
        continuation = i > 0 and state == previous and years[i]-years[i-1] == gap
        if continuation:
            size += 1
        else:
            if size and (not only_true or previous):
                runs.append(size)
            size = 1
        previous = state
    if size and (not only_true or previous):
        runs.append(size)
    return runs


def local_summary(years, values, threshold, direction):
    years, values = np.asarray(years, int), np.asarray(values, float)
    if len(np.unique(years)) != len(years):
        raise ValueError("Local metrics require one value per calendar year.")
    order = np.argsort(years)
    years, values = years[order], values[order]
    d = float(direction) * (values - float(threshold))
    tolerance = 1e-12 * max(abs(float(threshold)), 1e-12)
    bad = d > tolerance
    good = ~bad  # equality is compliant
    badruns, goodruns = mask_runs(bad, years), mask_runs(good, years)
    result = {"n_years": len(years), "N_good": int(good.sum()), "N_bad": int(bad.sum()),
              "noncompliance_fraction": float(bad.mean()),
              "Amax_bad": float(max(0.0, d.max())), "Amax_good": float(max(0.0, -d.min())),
              "Amax_bad_relative": float(max(0.0, d.max())/threshold),
              "Amax_good_relative": float(max(0.0, -d.min())/threshold),
              "Lmax_bad_years": max(badruns, default=0), "Lmax_good_years": max(goodruns, default=0),
              "ARL_bad_years": float(np.mean(badruns)) if badruns else np.nan,
              "ARL_good_years": float(np.mean(goodruns)) if goodruns else np.nan,
              "n_bad_runs": len(badruns), "n_good_runs": len(goodruns),
              "max_sampling_gap_years": int(np.diff(years).max()) if len(years)>1 else 0,
              "coverage_fraction": len(years)/(CFG.end_year-CFG.start_year+1),
              "tR_first_observed_years": np.nan, "tR_transition_lower_years": np.nan,
              "tR_transition_upper_years": np.nan, "tR_gap_uncertain": False,
              "recovery_status": "no_adverse_year", "worst_year": np.nan,
              "followup_years": np.nan}
    if bad.any():
        worst = int(np.argmax(d))
        result["worst_year"] = int(years[worst])
        result["followup_years"] = int(years[-1]-years[worst])
        after = np.flatnonzero((np.arange(len(years)) > worst) & good)
        if len(after):
            recovery = int(after[0])
            result["tR_first_observed_years"] = int(years[recovery]-years[worst])
            result["tR_transition_lower_years"] = int(years[recovery-1]-years[worst])
            result["tR_transition_upper_years"] = int(years[recovery]-years[worst])
            result["tR_gap_uncertain"] = bool(np.any(np.diff(years[worst:recovery+1]) > 1))
            result["recovery_status"] = "first_observed_compliance_after_worst_annual_state"
        else:
            result["recovery_status"] = "right_censored_no_observed_recovery"
    return result


def calculate_local(annual, value_col="value", threshold_shift=0.0):
    """Compute threshold-relative local indicators; negative threshold_shift is stricter."""
    records = []
    for (v, p), g in annual[annual.eligible_local].groupby(["VOR", "parameter"], sort=True):
        first = g.iloc[0]
        record = {"VOR": v, "parameter": p, "Type": first.Type, "unit": first.unit,
                  "threshold": first.threshold * (1.0 + threshold_shift * first.direction), "direction": first.direction,
                  "mean_sites_per_year": g.n_sites.mean(), "n_sites_max": g.n_sites.max(),
                  "basin": first.basin}
        record.update(local_summary(g.year, g[value_col], record["threshold"], record["direction"]))
        record["eligible_ranking"] = record["n_years"] >= CFG.min_years_local

        record["longest_bad_fraction_observed"] = record["Lmax_bad_years"] / record["n_years"]
        records.append(record)
    return pd.DataFrame(records)


def hvg_edges(values, reference_loop=False):
    """Strict directed HVG; reference_loop selects an independent direct-loop formulation. Gap lengths are not used; inputs must be finite."""
    x = np.asarray(values, float)
    if not np.isfinite(x).all():
        raise ValueError("HVG cannot receive nonfinite values.")
    n = len(x)
    if reference_loop:
        if n < 3:
            raise ValueError("Supplied R loops are not valid for n<3; excluded from parity test.")
        edges = []
        for i in range(n-2):
            alpha = 1
            edges.append((i, i+1))
            for j in range(2, n-i):
                crit, endpoint = x[i+alpha], x[i+j]
                if crit >= x[i]:
                    break
                if crit < endpoint:
                    edges.append((i, i+j))
                    alpha = j
        edges.append((n-2, n-1))
        return sorted(set(edges))
    edges = []
    for i in range(n-1):
        highest = -np.inf
        for j in range(i+1, n):
            if highest < min(x[i], x[j]):
                edges.append((i, j))
            highest = max(highest, x[j])
            if highest >= x[i]:
                break
    return edges


def hvg_transform(values, direction):
    """z_t = s_p (x_t - mean(x)); sign is relative to mean, not threshold."""
    x = np.asarray(values, float)
    return float(direction) * (x - x.mean())


def attack_auc(edges, n, seed, reps=20):
    """Random and adaptive targeted vertex-deletion AUC on the undirected HVG projection; component size is divided by the original node count."""
    graph = nx.Graph()
    graph.add_nodes_from(range(n)); graph.add_edges_from(edges)
    rng = np.random.default_rng(seed)
    def auc(order=None):
        h = graph.copy()
        sizes = [1.0]
        for step in range(n):
            target = order[step] if order is not None else min(h, key=lambda v: (-h.degree(v), v))
            h.remove_node(target)
            sizes.append(max((len(c) for c in nx.connected_components(h)), default=0)/n)
        trapezoid = getattr(np,"trapezoid",None) or np.trapz
        return float(trapezoid(sizes, dx=1/n))
    return float(np.mean([auc(rng.permutation(n)) for _ in range(reps)])), auc()


def hvg_properties(values, full=True, seed=42):
    x = np.asarray(values, float)
    n = len(x)
    if n < 3:
        raise ValueError("HVG properties require at least three values.")
    edges = hvg_edges(x)
    gd = ig.Graph(n=n, edges=edges, directed=True)
    gu = gd.as_undirected(mode="collapse")
    incoming, outgoing = np.array(gd.indegree()), np.array(gd.outdegree())
    dist = np.array(gd.distances(mode="out"), float)
    positive = np.isfinite(dist) & (dist > 0)
    inv = np.divide(1.0, dist, out=np.zeros_like(dist), where=positive)
    result = {
        "n_hvg_nodes": n,
        "DZO": float(np.sum(outgoing.max()-outgoing)/(n-1)**2),
        "Dens_D": len(edges)/(n*(n-1)),
        "AVPL_D": float(dist[positive].mean()),
        "GLE_D": float(inv.sum()/(n*(n-1))),
        "Tra_U": float(gu.transitivity_undirected(mode="zero")),
        "VAs_D": float(np.abs(outgoing-incoming).sum()/(2*len(edges))),
    }
    if full:
        total = incoming+outgoing
        bc = np.array(gd.betweenness(directed=True), float)
        pr = np.array(gd.pagerank(damping=0.85), float)
        local_eff = []
        for vertex in range(n):
            neighbors = gu.neighbors(vertex)
            if len(neighbors)<2:
                local_eff.append(0.0)
            else:
                sub = gu.induced_subgraph(neighbors)
                ds = np.array(sub.distances(), float)
                ok = np.isfinite(ds) & (ds > 0)
                local_eff.append(np.divide(1.,ds,out=np.zeros_like(ds),where=ok).sum()/(len(neighbors)*(len(neighbors)-1)))
        community = gu.community_leiden(objective_function="modularity", n_iterations=-1)
        infomap = gd.community_infomap(trials=10)
        rres, sres = attack_auc(edges, n, seed, CFG.attack_reps)
        result.update({
            "DZI": float(np.sum(incoming.max()-incoming)/(n-1)**2),
            "DZ_U": float(np.sum(total.max()-total)/((n-1)*(n-2))),
            "BZ_D": float(np.sum(bc.max()-bc)/((n-1)**2*(n-2))),
            "PRZ_D": float(np.sum(pr.max()-pr)),
            "Arcs": len(edges), "Diam_D": float(dist[positive].max()),
            "ALE_U": float(np.mean(local_eff)),
            "LeM_U": float(gu.modularity(community.membership)),
            "IM_D": float(gd.modularity(infomap.membership, directed=True)),
            "Mot_D": int(gd.motifs_randesu_no(size=3)),
            "RRes_U": rres, "SRes_U": sres,
        })
    return result


def alpha_centrality(graph, alpha=0.9):
    """(I-alpha A.T)c=1, with a unit exogenous vector.

    A sufficient absolute bound is checked numerically; nonconvergent cases
    are unavailable rather than silently replaced with different alpha.
    """
    nodes = list(graph)
    a = nx.to_numpy_array(graph, nodelist=nodes)
    spectral_radius = max(abs(linalg.eigvals(a))) if len(a) else 0
    if alpha * spectral_radius >= 1-1e-10:
        return {v: np.nan for v in nodes}, float(spectral_radius)
    vector = linalg.solve(np.eye(len(a))-alpha*a.T, np.ones(len(a)))
    return dict(zip(nodes, vector)), float(spectral_radius)


def steps_to_outlet(graph, v):
    """Number of recipient steps from v to the nearest reachable sink (outlet)."""
    lengths = nx.single_source_shortest_path_length(graph, v)
    sinks = [lengths[u] for u in lengths if graph.out_degree(u) == 0]
    return min(sinks) if sinks else np.nan


def hydro_properties(graph, meta):
    nodes = sorted(graph)
    g = graph.subgraph(nodes).copy()
    reverse = g.reverse(copy=False)
    ac, radius = alpha_centrality(g)

    idx = {v: i for i, v in enumerate(nodes)}
    gi = ig.Graph(n=len(nodes), edges=[(idx[u],idx[v]) for u,v in g.edges], directed=True)
    ig.set_random_number_generator(random.Random(CFG.seed))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        authorities = gi.authority_score(scale=True)
        hubs = gi.hub_score(scale=True)
    methods = {
        "DCI": dict(g.in_degree()), "DCO": dict(g.out_degree()), "DC": dict(g.degree()),
        "BC": nx.betweenness_centrality(g, normalized=True, weight=None),
        "CC_IN": nx.closeness_centrality(g, wf_improved=True),
        "CC_OUT": nx.closeness_centrality(reverse, wf_improved=True),
        "HC_IN": nx.harmonic_centrality(g), "HC_OUT": nx.harmonic_centrality(reverse),
        "PRC": nx.pagerank(g, alpha=0.85, max_iter=1000, tol=1e-12),
        "AC_alpha09": ac, "AUT": dict(zip(nodes,authorities)), "HBS": dict(zip(nodes,hubs)),

        "UPS": {v: len(nx.ancestors(g, v)) for v in nodes},
        "DOWN": {v: steps_to_outlet(g, v) for v in nodes},
    }
    table = pd.DataFrame({"VOR": nodes})
    for name, values in methods.items():
        table[name] = table.VOR.map(values)
    return table.merge(meta, on="VOR", validate="one_to_one"), radius


def calculate_hvgs(annual):
    records, rejected, edge_disagreements = [], [], []
    groups = list(annual.groupby(["VOR", "parameter"], sort=True))
    for count, ((v,p), group) in enumerate(groups):
        g = group.sort_values("year")
        row = g.iloc[0]
        if len(g)<CFG.min_years_hvg or pd.isna(row.direction) or row.Type not in {"Lake", "River"}:
            rejected.append({"VOR":v,"parameter":p,"n_years":len(g),
                             "reason":"short_series_or_unknown_direction_or_nonregistered_VOR"})
            continue
        z = hvg_transform(g.value, row.direction)
        if hvg_edges(z, reference_loop=True) != hvg_edges(z):
            edge_disagreements.append({"VOR":v,"parameter":p})
        seed = CFG.seed+int(hashlib.sha256(f"{v}|{p}".encode()).hexdigest()[:8],16)
        ig.set_random_number_generator(random.Random(seed))
        features = hvg_properties(z, full=True, seed=seed)
        records.append({"VOR":v,"parameter":p,"unit":row.unit,"Type":row.Type,
                        "n_years":len(g),"n_sites_max":g.n_sites.max(),"basin":row.basin,
                        "year_min":int(g.year.min()),"year_max":int(g.year.max()),
                        "n_missing_internal_years":int(g.year.max()-g.year.min()+1-len(g)),
                        "constant_series":bool(np.ptp(z)==0), **features})
        if count % 500 == 0:
            LOG.info("HVG progress: %s / %s body-parameter series", count, len(groups))
    return pd.DataFrame(records), pd.DataFrame(rejected), pd.DataFrame(edge_disagreements)


def coverage_tables(annual, site, raw, meta):
    groups = annual.groupby(["VOR", "parameter"], sort=True)
    rows = []
    for (v,p), g in groups:
        y = np.sort(g.year.to_numpy())
        r = g.iloc[0]
        rows.append({"VOR":v,"parameter":p,"Type":r.Type,"unit":r.unit,
                     "n_years":len(y),"year_min":y[0],"year_max":y[-1],
                     "n_internal_missing_years":y[-1]-y[0]+1-len(y),
                     "max_gap_years":int(np.diff(y).max()) if len(y)>1 else 0,
                     "sites_max_per_year":int(g.n_sites.max()),
                     "sites_mean_per_year":g.n_sites.mean(),
                     "threshold_available":bool(r.eligible_local),
                     "eligible_HVG":len(y)>=CFG.min_years_hvg and r.Type in {"River","Lake"},
                     "eligible_local_ranking":len(y)>=CFG.min_years_local and r.eligible_local})
    coverage = pd.DataFrame(rows)
    parameter = coverage.groupby(["parameter","unit"],as_index=False).agg(
        body_parameter_combinations=("VOR","size"), median_years=("n_years","median"),
        min_years=("n_years","min"), max_years=("n_years","max"),
        eligible_HVG=("eligible_HVG","sum"),
        eligible_local_ranking=("eligible_local_ranking","sum"))
    parameter["site_year_records"] = parameter.parameter.map(site.groupby("parameter").size())
    parameter["body_year_records"] = parameter.parameter.map(annual.groupby("parameter").size())
    rawcounts = raw.groupby(["ANYAG_NEV_ANGOL","MERTEKEGYSEG_NEV"],dropna=False).size().reset_index(name="raw_rows")
    return coverage, parameter, rawcounts


def select_nonredundant(df, priority, cutoff):
    """Outcome-blind pruning, ordered by predeclared hydrological rationale."""
    chosen, rows = [], []
    for metric in priority:
        if metric not in df or df[metric].dropna().nunique()<2:
            rows.append({"metric":metric,"retained":False,"reason":"constant_or_unavailable","with_metric":"","rho":np.nan})
            continue
        redundant = None
        for previous in chosen:
            pair = df[[metric,previous]].dropna()
            rho = stats.spearmanr(pair[metric],pair[previous]).statistic if len(pair)>2 else np.nan
            if np.isfinite(rho) and abs(rho)>=cutoff:
                redundant = (previous,rho); break
        if redundant:
            rows.append({"metric":metric,"retained":False,"reason":"redundant_at_prespecified_cutoff",
                         "with_metric":redundant[0],"rho":redundant[1]})
        else:
            chosen.append(metric)
            rows.append({"metric":metric,"retained":True,"reason":"theory_priority_and_nonredundancy","with_metric":"","rho":np.nan})
    return chosen, pd.DataFrame(rows)


def control_matrix(df, adjust=True):
    """Within-parameter control design. Missingness is explicit, not listwise hidden."""
    cols = [np.ones(len(df))]
    names = ["intercept"]
    if adjust:
        numeric = ["catchment_km2","river_length_km","lake_area_km2","depth_m",
                   "artificial","modified","n_years","n_sites_max"]
        for field in numeric:
            if field not in df:
                continue
            x = pd.to_numeric(df[field],errors="coerce").to_numpy(dtype=float,copy=True)
            missing = ~np.isfinite(x)
            if missing.all():
                continue
            if missing.any() and not missing.all():
                cols.append(missing.astype(float)); names.append(field+"_missing")
            x[missing] = np.nanmedian(x)
            if np.ptp(x)>0:
                cols.append(stats.rankdata(x)); names.append(field+"_rank")
        for field in ["Type","altitude_class","geochemistry","slope_class","substrate_class"]:
            if field not in df:
                continue
            dummy = pd.get_dummies(df[field].fillna("Unknown").astype(str), drop_first=True, dtype=float)
            for name in dummy:
                cols.append(dummy[name].to_numpy()); names.append(field+":"+name)
    design = np.column_stack(cols)
    # Rank-deficient dummy or type-specific covariates are removed explicitly.
    _, rr, pivot = linalg.qr(design,mode="economic",pivoting=True)
    tol = max(design.shape)*np.finfo(float).eps*abs(rr[0,0]) if rr.size else 0
    rank = int(np.sum(np.abs(np.diag(rr))>tol))
    selected = sorted(pivot[:rank])
    return design[:,selected], [names[i] for i in selected]


def rank_association(df, xname, yname, adjust, seed):
    """Partial rank correlation; cluster-robust slope test by RBMP sub-basin.

    CI: paired sub-basin bootstrap of FIXED residual ranks, conditional on
    fitted controls/ranks. It is not a fully refitted bootstrap Spearman CI.
    Shared river connectivity can cross basin boundaries: dependence is reduced
    by blocks but not assumed eliminated. Asymptotic raw p is supplementary.
    """
    valid = np.isfinite(pd.to_numeric(df[xname],errors="coerce")) & np.isfinite(pd.to_numeric(df[yname],errors="coerce"))
    d = df.loc[valid].copy()
    result = {"n":len(d),"rho":np.nan,"p_asymptotic_raw":np.nan,"p_cluster":np.nan,
              "ci_low":np.nan,"ci_high":np.nan,"n_basins":0,"n_controls":0,
              "status":"insufficient_n","ci_method":"fixed_rank_residual_subbasin_bootstrap",
              "controls":""}
    if len(d)<CFG.correlation_min_n:
        return result
    x, y = d[xname].to_numpy(float),d[yname].to_numpy(float)
    if np.ptp(x)==0 or np.ptp(y)==0:
        result["status"]="constant_metric"; return result
    raw = stats.spearmanr(x,y)
    result["p_asymptotic_raw"] = float(raw.pvalue)
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    design, names = control_matrix(d, adjust)
    result["n_controls"] = design.shape[1]-1
    result["controls"] = ";".join(names)
    if len(d)<=design.shape[1]+8:
        result["status"]="insufficient_residual_df"; return result
    rx -= design @ linalg.lstsq(design,rx)[0]
    ry -= design @ linalg.lstsq(design,ry)[0]
    if np.linalg.norm(rx)<1e-9 or np.linalg.norm(ry)<1e-9:
        result["status"]="constant_after_adjustment"; return result
    rho = float(np.dot(rx,ry)/np.sqrt(np.dot(rx,rx)*np.dot(ry,ry)))
    result["rho"] = rho
    basins = d.basin.fillna("Unknown").astype(str).to_numpy()
    labels, inverse = np.unique(basins,return_inverse=True)
    nbasins = len(labels)
    result["n_basins"] = nbasins

    # All unknown sub-basin observations share one cluster.
    if nbasins<CFG.correlation_min_basins:
        result["status"]="effect_only_insufficient_known_subbasins"
        return result
    # A binary outlet indicator may have only four zeros in three sub-basins.
    # Many total sub-basins do not supply independent support for that level.
    # Retain its descriptive effect, but do not manufacture precise p/q values.
    levels = np.unique(x)
    if len(levels) <= 3:
        support = [(int((x==level).sum()),len(np.unique(basins[x==level]))) for level in levels]
        if any(nlevel<10 or nbasinlevel<5 for nlevel,nbasinlevel in support):
            result["status"]="effect_only_rare_predictor_level"
            return result
    slope = np.dot(rx,ry)/np.dot(rx,rx)
    residual = ry-slope*rx
    cluster_scores = np.bincount(inverse,weights=rx*residual,minlength=nbasins)
    correction = nbasins/(nbasins-1)*(len(d)-1)/(len(d)-design.shape[1]-1)
    variance = correction*np.dot(cluster_scores,cluster_scores)/np.dot(rx,rx)**2
    if variance>0:
        result["p_cluster"] = float(2*stats.t.sf(abs(slope/np.sqrt(variance)),df=nbasins-1))
    elif slope != 0:
        result["p_cluster"] = 0.0
    moments = np.column_stack([np.bincount(inverse,weights=a,minlength=nbasins)
                              for a in [np.ones(len(d)),rx,ry,rx*rx,ry*ry,rx*ry]])
    rng = np.random.default_rng(seed)
    sample = rng.integers(0,nbasins,size=(CFG.bootstrap_reps,nbasins))
    boot = moments[sample].sum(axis=1)
    n,sx,sy,sxx,syy,sxy = boot.T
    den = np.sqrt(np.maximum(0,sxx-sx*sx/n)*np.maximum(0,syy-sy*sy/n))
    estimates = np.divide(sxy-sx*sy/n,den,out=np.full_like(den,np.nan),where=den>0)
    if np.isfinite(estimates).sum()>CFG.bootstrap_reps/2:
        result["ci_low"],result["ci_high"] = np.nanquantile(estimates,[.025,.975]).tolist()
    result["status"]="ok_cluster_robust"
    return result


def correlation_table(outcomes, hydro, xmetrics, ymetrics, family, adjust):
    additional = [x for x in hydro.columns if x not in outcomes.columns or x=="VOR"]
    df = outcomes.merge(hydro[additional],on="VOR",validate="many_to_one")
    records = []
    for p,group in df.groupby("parameter",sort=True):
        for x in xmetrics:
            for y in ymetrics:
                seed = CFG.seed+int(hashlib.sha256(f"{family}|{p}|{x}|{y}|{adjust}".encode()).hexdigest()[:8],16)
                result = rank_association(group,x,y,adjust,seed)
                result.update({"parameter":p,"hydro_metric":x,"outcome":y,
                               "family":family,"adjusted":adjust})
                records.append(result)
    table = pd.DataFrame(records)
    table["q_BH"] = np.nan
    if len(table):
        valid = np.isfinite(table.p_cluster)
        if valid.any():
            table.loc[valid,"q_BH"] = stats.false_discovery_control(table.loc[valid,"p_cluster"].to_numpy(),method="bh")
        table["FDR_significant"] = table.q_BH.lt(.05)
    return table


def ranking(local, meta):
    eligible = local[local.eligible_ranking].copy()
    if eligible.empty:
        return pd.DataFrame()
    rank = eligible.groupby("VOR",as_index=False).agg(
        score=("longest_bad_fraction_observed","mean"),
        mean_noncompliance=("noncompliance_fraction","mean"),
        n_parameters=("parameter","nunique"), minimum_years=("n_years","min"),
        mean_years=("n_years","mean"))
    rank["eligible"] = rank.n_parameters>=CFG.min_parameters_ranking
    rank["rank"] = np.nan
    rank.loc[rank.eligible,"rank"] = rank.loc[rank.eligible,"score"].rank(ascending=False,method="min")
    cutoff = rank.loc[rank.eligible,"score"].quantile(CFG.critical_quantile)
    rank["critical_cutoff"] = cutoff
    rank["critical"] = rank.eligible & rank.score.ge(cutoff)
    return rank.merge(meta[["VOR","name","Type","Latitude","Longitude","basin"]],on="VOR").sort_values(["eligible","rank","VOR"],ascending=[False,True,True])


def ranking_sensitivity(annual, local, base_rank, meta):
    """Threshold +/-5%, site weighting, year omissions, common parameter set.

    Year-omission intervals are sensitivity envelopes, not confidence intervals.
    Missing years split runs rather than fabricating continuity.
    """
    records = []
    for shift in CFG.threshold_shifts:
        label = f"threshold_{'stricter' if shift<0 else 'more_lenient'}_{abs(shift):.0%}"
        alt = ranking(calculate_local(annual,threshold_shift=shift),meta)
        for _,r in alt.iterrows():
            records.append({"scenario":label,"VOR":r.VOR,"score":r.score,"rank":r['rank'],"critical":r.critical})
    alt = ranking(calculate_local(annual,value_col="value_sample_weighted"),meta)
    for _,r in alt.iterrows():
        records.append({"scenario":"sample_count_weighted","VOR":r.VOR,"score":r.score,"rank":r['rank'],"critical":r.critical})
    eligible = local[local.eligible_ranking]
    frequencies = eligible.groupby("parameter").VOR.nunique()
    common = frequencies[frequencies>=.8*frequencies.max()].index.tolist()
    alt = ranking(local[local.parameter.isin(common)],meta)
    for _,r in alt.iterrows():
        records.append({"scenario":"common_parameters","VOR":r.VOR,"score":r.score,"rank":r['rank'],"critical":r.critical})
    series = [(v,p,g.year.to_numpy(int),g.adverse_deviation.to_numpy(float))
              for (v,p),g in annual[annual.eligible_local].groupby(["VOR","parameter"],sort=True)]
    for omit in range(CFG.start_year,CFG.end_year+1):
        scores = []
        for v,p,y,d in series:
            keep = y!=omit
            yy,dd = y[keep],d[keep]
            if len(yy)>=CFG.min_years_local:
                scores.append((v,p,max(mask_runs(dd>0,yy),default=0)/len(yy)))
        frame = pd.DataFrame(scores,columns=["VOR","parameter","score"])
        scenario = frame.groupby("VOR").agg(score=("score","mean"),n_parameters=("parameter","nunique"))
        scenario = scenario[scenario.n_parameters>=CFG.min_parameters_ranking]
        cutoff = scenario.score.quantile(CFG.critical_quantile)
        scenario["rank"] = scenario.score.rank(ascending=False,method="min")
        for v,r in scenario.iterrows():
            records.append({"scenario":f"omit_year_{omit}","VOR":v,"score":r.score,"rank":r['rank'],"critical":bool(r.score>=cutoff)})
    detail = pd.DataFrame(records)
    yr = detail[detail.scenario.str.startswith("omit_year_")]
    envelope = yr.groupby("VOR",as_index=False).agg(
        rank_min=("rank","min"),rank_max=("rank","max"),
        rank_median=("rank","median"),critical_year_omission_fraction=("critical","mean"),
        n_omission_scenarios=("scenario","size"))

    base = base_rank[base_rank.eligible].set_index("VOR")
    base_top5 = set(base.nsmallest(5,"rank").index)
    base_crit = set(base.index[base.critical])
    agreement = []
    for scenario,g in detail.groupby("scenario",sort=False):
        g = g.dropna(subset=["rank"]).set_index("VOR")
        common_ids = base.index.intersection(g.index)
        if len(common_ids)<3:
            continue
        rho = stats.spearmanr(base.loc[common_ids,"score"],g.loc[common_ids,"score"]).statistic
        crit = set(g.index[g.critical.astype(bool)])
        union = base_crit|crit
        top5 = set(g.nsmallest(5,"rank").index)
        agreement.append({"scenario":scenario,"n_bodies":len(common_ids),"spearman_score_vs_baseline":rho,
                          "jaccard_critical_sets":len(base_crit&crit)/len(union) if union else np.nan,
                          "top5_overlap":len(base_top5&top5)/5})
    return detail,envelope,common,pd.DataFrame(agreement)


def upstream_conditioned_assessment(annual, graph):
    """Empirical annual response conditional on observed upstream adversity.

    C(v,p)=sum[good(v,p,y)*adverse_observed_upstream(v,p,y)] / n_exposed.
    Direct predecessors only; no flow-volume weights, travel time or causal
    interpretation. Partial upstream coverage and complete coverage reported
    separately. Moving-block bootstrap CI preserves short temporal sequences.
    """
    eligible = annual[annual.eligible_local].copy()
    eligible["good"] = eligible.adverse_deviation <= 1e-12*eligible.threshold.abs().clip(lower=1e-12)
    lookup = {(v,p):g.set_index("year") for (v,p),g in eligible.groupby(["VOR","parameter"])}
    details, summaries = [], []
    for (v,p),g in sorted(lookup.items()):
        predecessors = sorted(graph.predecessors(v))
        if not predecessors:
            continue
        events = []
        for year,r in g.iterrows():
            states = []
            for u in predecessors:
                neighbor = lookup.get((u,p))
                if neighbor is not None and year in neighbor.index:
                    states.append(bool(neighbor.loc[year,"good"]))
            if not states:
                continue
            exposed = not all(states)
            complete = len(states)==len(predecessors)
            event = {"VOR":v,"parameter":p,"year":int(year),"downstream_good":bool(r.good),
                     "observed_upstream_adverse":exposed,"upstream_count":len(predecessors),
                     "observed_upstream_count":len(states),"complete_upstream_coverage":complete}
            events.append(event); details.append(event)
        if not events:
            continue
        frame = pd.DataFrame(events)
        a = frame.observed_upstream_adverse.to_numpy(bool)
        good = frame.downstream_good.to_numpy(bool)
        complete = frame.complete_upstream_coverage.to_numpy(bool)
        seed = CFG.seed+int(hashlib.sha256(f"upstream|{v}|{p}".encode()).hexdigest()[:8],16)
        rng = np.random.default_rng(seed)
        # Reindex to calendar years to retain actual missingness and gap pattern.
        calendar = frame.set_index("year").reindex(range(CFG.start_year,CFG.end_year+1))
        exposures = calendar.observed_upstream_adverse.eq(True).to_numpy()
        outcomes = calendar.downstream_good.eq(True).to_numpy()
        estimates = []
        length = len(calendar); blocklen = 3
        for _ in range(CFG.bootstrap_reps):
            starts = rng.integers(0,length-blocklen+1,size=math.ceil(length/blocklen))
            indices = np.concatenate([np.arange(s,s+blocklen) for s in starts])[:length]
            ae = exposures[indices]
            if ae.sum()>0:
                estimates.append(float(outcomes[indices][ae].mean()))
        n_exposed = int(a.sum())
        lo,hi = np.quantile(estimates,[.025,.975]) if n_exposed>=CFG.min_upstream_exposed_years and len(estimates)>CFG.bootstrap_reps/2 else (np.nan,np.nan)
        clean = (~a) & complete  # unobserved upstream bodies cannot be called clean
        exposedcomplete = a & complete
        summaries.append({"VOR":v,"parameter":p,"n_matched_years":len(frame),
                          "n_exposed_years":n_exposed,
                          "conditional_compliance":float(good[a].mean()) if n_exposed else np.nan,
                          "ci_block_low":lo,"ci_block_high":hi,
                          "n_complete_exposed_years":int(exposedcomplete.sum()),
                          "complete_coverage_conditional_compliance":float(good[exposedcomplete].mean()) if exposedcomplete.any() else np.nan,
                          "n_completely_observed_clean_years":int(clean.sum()),
                          "compliance_when_all_upstream_observed_clean":float(good[clean].mean()) if clean.any() else np.nan,
                          "mean_upstream_observation_fraction":float((frame.observed_upstream_count/frame.upstream_count).mean()),
                          "adequate_exposed_years":n_exposed>=CFG.min_upstream_exposed_years,
                          "interpretation":"annual empirical compliance under observed upstream adversity; noncausal"})
    return pd.DataFrame(summaries),pd.DataFrame(details)


def hvg_sensitivity(annual, hvg):
    """Stratified random subset; random and contiguous-block omissions.

    Export replicate metrics and between-series Spearman concordance plus
    absolute changes, NOT a single example's visual clustering agreement.
    """
    rng = np.random.default_rng(CFG.seed)
    good = hvg[hvg.n_years>=CFG.sensitivity_min_years]
    selected = []
    for p,g in good.groupby("parameter",sort=True):
        if len(g)>CFG.sensitivity_series_per_parameter:
            g = g.iloc[np.sort(rng.choice(len(g),CFG.sensitivity_series_per_parameter,replace=False))]
        selected.extend(zip(g.VOR,g.parameter))
    lookup = {(v,p):g for (v,p),g in annual.groupby(["VOR","parameter"])}
    reps = []
    for count,(v,p) in enumerate(selected):
        g = lookup[(v,p)].sort_values("year")
        z = hvg_transform(g.value,g.direction.iloc[0])
        base = hvg_properties(z,full=False)
        for frac in CFG.deletion_fractions:
            remove = max(1,int(round(len(z)*frac)))
            if len(z)-remove<3:
                continue
            for method in ["random","contiguous_block"]:
                for replicate in range(CFG.sensitivity_reps):
                    if method=="random":
                        deleted = rng.choice(len(z),remove,replace=False)
                    else:
                        start = rng.integers(0,len(z)-remove+1)
                        deleted = np.arange(start,start+remove)
                    mask = np.ones(len(z),bool);mask[deleted]=False
                    alt = hvg_properties(z[mask],full=False)
                    for metric in HVG_FAST:
                        reps.append({"VOR":v,"parameter":p,"method":method,
                                     "deletion_fraction":frac,"replicate":replicate,"metric":metric,
                                     "baseline_n":len(z),"retained_n":int(mask.sum()),
                                     "baseline":base[metric],"perturbed":alt[metric],
                                     "absolute_change":abs(alt[metric]-base[metric])})
        if count % 100==0:
            LOG.info("Sensitivity progress: %s / %s selected series",count,len(selected))
    detail = pd.DataFrame(reps)
    result = []
    if not detail.empty:
        averaged = detail.groupby(["parameter","method","deletion_fraction","metric","VOR"],as_index=False).agg(
            baseline=("baseline","first"),perturbed=("perturbed","median"),
            absolute_change=("absolute_change","median"))
        for keys,g in averaged.groupby(["parameter","method","deletion_fraction","metric"],sort=True):
            rho = stats.spearmanr(g.baseline,g.perturbed).statistic if len(g)>2 and g.baseline.nunique()>1 and g.perturbed.nunique()>1 else np.nan
            result.append(dict(zip(["parameter","method","deletion_fraction","metric"],keys),
                               n_series=len(g),rho_baseline_vs_median_perturbed=rho,
                               median_absolute_change=g.absolute_change.median(),
                               p95_absolute_change=g.absolute_change.quantile(.95)))
    # Annual station-weight sensitivity and two-year bin aggregation.
    aggregation = []
    for _,row in hvg.iterrows():
        g = lookup[(row.VOR,row.parameter)].sort_values("year")
        for mode in ["sample_count_weighted","two_year_bins"]:
            if mode=="sample_count_weighted":
                vals = g.value_sample_weighted.to_numpy()
            else:
                tmp = g.assign(bin=(g.year-CFG.start_year)//2).groupby("bin").value.mean()
                vals = tmp.to_numpy()
            if len(vals)<3:
                continue
            alt = hvg_properties(hvg_transform(vals,row.get("direction",g.direction.iloc[0])),full=False)
            for metric in HVG_FAST:
                aggregation.append({"VOR":row.VOR,"parameter":row.parameter,"scheme":mode,"metric":metric,
                                    "baseline":row[metric],"perturbed":alt[metric],"n_baseline":len(g),"n_perturbed":len(vals)})
    return detail,pd.DataFrame(result),pd.DataFrame(aggregation),pd.DataFrame(selected,columns=["VOR","parameter"])


def optional_rbmp(root, rank):
    """Join optional official status by VOR. Required columns: VOR, status, assessment_year. Absence of the optional input does not create a result table."""
    p = Path(root)/"official_water_body_status.csv"
    if not p.is_file():
        return pd.DataFrame(),"not included: no official-status input is distributed"
    rbmp = pd.read_csv(p)
    required = {"VOR","status","assessment_year"}
    if not required.issubset(rbmp):
        raise ValueError(f"official_water_body_status.csv requires columns {sorted(required)}")
    rbmp["assessment_year"] = pd.to_numeric(rbmp.assessment_year,errors="raise")
    if rbmp.VOR.duplicated().any():
        raise ValueError("Supply one RBMP assessment row per VOR, for an explicitly selected assessment period.")
    return rank.merge(rbmp,on="VOR",how="left",validate="one_to_one"),"official status joined by VOR; assessment_year identifies the supplied assessment period"


def series_length_sensitivity(hvg, hydro, xmetrics, annual):
    tables = []

    # already-truncated >=8-year table evaluates a six-year minimum.
    extra=[]
    for (v,p),g in annual.groupby(["VOR","parameter"],sort=True):
        if 6<=len(g)<CFG.min_years_hvg and g.Type.iloc[0] in {"River","Lake"} and pd.notna(g.direction.iloc[0]):
            extra.append({"VOR":v,"parameter":p,"n_years":len(g),"n_sites_max":g.n_sites.max(),
                          "Type":g.Type.iloc[0],"basin":g.basin.iloc[0],
                          **hvg_properties(hvg_transform(g.value,g.direction.iloc[0]),full=False)})
    all_lengths=pd.concat([hvg,pd.DataFrame(extra)],ignore_index=True)
    for cutoff in [6,8,10,15]:
        table = correlation_table(all_lengths[all_lengths.n_years>=cutoff],hydro,xmetrics,
                                  HVG_FAST,"HVG_length_sensitivity",True)
        table["minimum_years"] = cutoff
        tables.append(table)
    return pd.concat(tables,ignore_index=True) if tables else pd.DataFrame()


def quality_profile_tsne(annual, meta, output):
    """Exploratory embedding; typology is used only as a plotting label."""
    output = Path(output)
    pivot, features, imputed_fraction = quality_profile_matrix(annual)
    settings = {"n_components": 2, "random_state": CFG.seed,
                "perplexity": CFG.tsne_perplexity,
                "learning_rate": CFG.tsne_learning_rate,
                "max_iter": CFG.tsne_iterations, "init": "pca",
                "metric": "euclidean", "method": "barnes_hut", "angle": .5}
    import inspect
    if "max_iter" not in inspect.signature(TSNE).parameters:
        settings["n_iter"] = settings.pop("max_iter")
    if len(pivot) <= CFG.tsne_perplexity or len(pivot.columns) < 2:
        return {"executed": False, "reason": "insufficient_profile_size"}
    coordinates = TSNE(**settings).fit_transform(features)
    embedding = pd.DataFrame(coordinates, index=pivot.index,
                             columns=["TSNE1", "TSNE2"]).reset_index()
    embedding = embedding.merge(meta, on="VOR", how="left", validate="one_to_one")
    embedding["n_profile_parameters"] = pivot.notna().sum(axis=1).reindex(embedding.VOR).to_numpy()
    write_csv(embedding, output.parent / "tables" / "quality_profile_tsne.csv")
    write_csv(pivot.reset_index(), output.parent / "tables" / "quality_profile_input_matrix.csv")
    plot_tsne(embedding, output / "quality_profile_tsne",
              "Physicochemical profiles; typology used only for colouring")
    return {"corrected": {"n_bodies": len(pivot),
            "input_parameters": pivot.columns.tolist(),
            "imputed_fraction": imputed_fraction,
            "preparation": "mean oriented relative deviation per body/parameter; at least 3 parameters per body; parameters in at least 10 bodies; parameter-median imputation; z-scaling",
            "settings": settings,
            "interpretation": "exploratory; not a significance or causal test"}}


def save_figure(fig, stem):
    fig.savefig(str(stem)+".pdf",bbox_inches="tight",pad_inches=.12)
    fig.savefig(str(stem)+".png",bbox_inches="tight",dpi=180,pad_inches=.12)
    plt.close(fig)


def plot_tsne(df, stem, title):
    fig,ax = plt.subplots(figsize=(7.5,6.2))
    palette = {"River":"#0072B2","Lake":"#D55E00","Other":"#666666"}
    for typ,g in df.groupby("Type",sort=True):
        ax.scatter(g.TSNE1,g.TSNE2,s=13 if len(df)>2000 else 24,c=palette.get(typ,"gray"),
                   alpha=.55,label=f"{typ} (n={len(g)})",rasterized=len(df)>2000)
    ax.set(xlabel="t-SNE coordinate 1 (arbitrary units)",ylabel="t-SNE coordinate 2 (arbitrary units)",title=title)
    ax.legend(frameon=False);ax.grid(alpha=.15)
    save_figure(fig,stem)


def example_vor(annual):
    """CFG.selected_vor (Lake Balaton, AIH049) if present, else the body with most BOD years."""
    if annual.VOR.eq(CFG.selected_vor).any():
        return CFG.selected_vor
    b = annual[annual.parameter.eq("BOD (5)") & annual.eligible_local]
    return b.groupby("VOR").size().idxmax() if len(b) else annual.VOR.iloc[0]


def figures(annual, local, hydro, ranking_df, rawcorr, partialcorr, root):
    root = Path(root);root.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":12,"axes.titlesize":14,"axes.labelsize":12,
                         "xtick.labelsize":11,"ytick.labelsize":11,"pdf.fonttype":42,
                         "ps.fonttype":42})

    fig,ax = plt.subplots(figsize=(13,8));ax.set_axis_off()
    boxes = {
        "inputs":(.5,.94,"Water-quality annual summaries + sites\nHydrological recipients + typology + operational thresholds"),
        "pre":(.5,.76,"Unit validation; source-copy removal; site/year aggregation\nStable VOR keys; coverage and exclusions"),
        "local":(.18,.51,"Annual compliance and recovery\nNine local indicators"),
        "hydro":(.50,.51,"Hydrological graph\nVOR nodes; upstream → recipient edges\nCentrality + redundancy selection"),
        "hvg":(.83,.51,"Annual water-quality series → HVG\nObservation nodes; earlier → later edges\nTemporal graph metrics"),
        "stats":(.50,.27,"Associations; typology/coverage controls; sub-basin uncertainty; FDR\nDeletion and aggregation sensitivity; observed upstream adversity"),
        "outputs":(.50,.07,"Reproducible rankings and critical upper tail\nMonitoring priorities; RBMP comparison when official input is supplied")}
    for key,(x,y,label) in boxes.items():
        ax.text(x,y,label,ha="center",va="center",fontsize=11,
                bbox={"boxstyle":"round,pad=.7","facecolor":"#eef5fa","edgecolor":"#386b86"},transform=ax.transAxes)
    for source,target in [("inputs","pre"),("pre","local"),("pre","hydro"),("pre","hvg"),
                          ("local","stats"),("hydro","stats"),("hvg","stats"),("stats","outputs")]:
        sx,sy,_=boxes[source];tx,ty,_=boxes[target]
        ax.annotate("",xy=(tx,ty+.065),xytext=(sx,sy-.065),xycoords="axes fraction",
                    arrowprops={"arrowstyle":"->","color":"#386b86","lw":1.4})
    save_figure(fig,root/"workflow")
    # Annual example: threshold deviation and HVG on the SAME data representation.
    candidates = annual[annual.VOR.eq(example_vor(annual)) & annual.eligible_local]
    if len(candidates):
        counts = candidates.groupby("parameter").size().sort_values(ascending=False)
        p = "BOD (5)" if "BOD (5)" in counts else counts.index[0]
        g = candidates[candidates.parameter.eq(p)].sort_values("year")
        fig,ax = plt.subplots(figsize=(11,5))
        ax.plot(g.year,g.adverse_deviation,"o-",color="#0072B2",lw=1.5)
        ax.axhline(0,color="#D55E00",ls="--",label="Operational threshold; equality compliant")
        sc,ul=(1000.,"mg/L") if g.unit.iloc[0]=="g/L" else (1.,g.unit.iloc[0])
        ax.lines[0].set_ydata(g.adverse_deviation*sc)
        ax.relim();ax.autoscale_view()
        ax.set(xlabel="Calendar year",ylabel=f"Oriented deviation d = s(x − T) [{ul}]",
               title=f"{example_vor(annual)}: {p}; annual water-body values (d > 0: non-compliant)")
        ax.legend(frameon=False);ax.grid(alpha=.2)
        save_figure(fig,root/"annual_threshold_deviation")
        z = hvg_transform(g.value,g.direction.iloc[0]);edges=hvg_edges(z)

        sc,ul=(1000.,"mg/L") if g.unit.iloc[0]=="g/L" else (1.,g.unit.iloc[0])
        zz=z*sc;yrs=g.year.to_numpy()
        fig,(top,bottom)=plt.subplots(2,1,figsize=(10,7.2),sharex=True,gridspec_kw={"height_ratios":[1.3,1]})
        base=zz.min()-.08*np.ptp(zz) if np.ptp(zz)>0 else zz.min()-1
        for i,j in edges:
            h=min(zz[i],zz[j])
            top.plot([yrs[i],yrs[j]],[h,h],color="#D55E00",lw=.9,alpha=.75,zorder=1)
        top.vlines(yrs,base,zz,color="#0b3d5c",lw=5,zorder=2)
        top.scatter(yrs,zz,s=30,color="#0b3d5c",zorder=3)
        top.axhline(0,color="gray",ls="--",lw=1)
        top.set_ylabel(f"z = s(x − mean x) [{ul}]",fontsize=13)
        top.set_title("(a) Annual series z; orange: horizontal lines of sight (HVG arcs)",fontsize=13)
        for i,j in edges:
            mid,r=(yrs[i]+yrs[j])/2,(yrs[j]-yrs[i])/2
            t=np.linspace(0,np.pi,60)
            bottom.plot(mid+r*np.cos(t[::-1]),r*np.sin(t),color="#386b86",lw=1)
        bottom.scatter(yrs,np.zeros(len(yrs)),s=110,color="#c5deeb",edgecolor="#0b3d5c",zorder=3)
        indeg=np.bincount([j for _,j in edges],minlength=len(z));outdeg=np.bincount([i for i,_ in edges],minlength=len(z))
        for k,yv in enumerate(yrs):
            bottom.text(yv,-.9,f"{outdeg[k]}/{indeg[k]}",ha="center",va="top",fontsize=10)
        bottom.set_ylim(-1.8,max(1,(yrs.max()-yrs.min())/2+.5));bottom.set_yticks([])
        bottom.set_xlabel("Calendar year (node = one annual value; arcs earlier → later)",fontsize=13)
        bottom.set_title("(b) Directed HVG; labels: out-degree/in-degree",fontsize=13)
        for spine in ["left","right","top"]: bottom.spines[spine].set_visible(False)
        fig.suptitle(f"{example_vor(annual)}: {p}",fontsize=15);fig.tight_layout()
        save_figure(fig,root/"annual_visibility_graph")
        # Individual measurements unavailable. Figure 4 is an explicit
        # site-annual comparison, never labelled raw observations.
    # Monitoring-centroid network schematic; no false river geometry.
    observed = hydro[hydro.observed & hydro.Latitude.notna() & hydro.Longitude.notna()]
    fig,ax=plt.subplots(figsize=(9,7))
    positions=observed.set_index("VOR")[["Longitude","Latitude"]].to_dict("index")
    edgefile=Path(root).parent/"tables"/"recipient_links.csv"
    if edgefile.is_file():
        for _,e in pd.read_csv(edgefile).iterrows():
            if e['from'] in positions and e['to'] in positions:
                u,w=positions[e['from']],positions[e['to']]
                ax.plot([u['Longitude'],w['Longitude']],[u['Latitude'],w['Latitude']],color="#aaaaaa",lw=.5,alpha=.45,zorder=1)
    for typ,g in observed.groupby("Type"):
        colour={"River":"#0072B2","Lake":"#D55E00","Other":"#666666"}.get(typ,"gray")
        ax.scatter(g.Longitude,g.Latitude,s=20,label=f"{typ} (n={len(g)})",c=colour,zorder=2,alpha=.8)
    mid=observed.Latitude.median();ax.set_aspect(1/max(.1,np.cos(np.deg2rad(mid))))
    ax.autoscale_view();xmin,xmax=ax.get_xlim();ymin,ymax=ax.get_ylim()
    ax.annotate("N",xy=(xmax-.15,ymax-.20),xytext=(xmax-.15,ymax-.65),ha="center",fontsize=14,
                arrowprops={"arrowstyle":"-|>","color":"black"})
    km50=50/(111.32*np.cos(np.deg2rad(mid)))
    ax.plot([xmin+.2,xmin+.2+km50],[ymin+.15,ymin+.15],color="black",lw=3)
    ax.text(xmin+.2+km50/2,ymin+.20,"50 km (approx.)",ha="center")
    ax.set(xlabel="Longitude [degrees E]",ylabel="Latitude [degrees N]",
           title="Hungarian monitoring locations and recipient connectivity\nLinks between monitoring centroids are schematic")
    ax.legend(frameon=False);ax.grid(alpha=.2)
    save_figure(fig,root/"monitoring_network")
    # Histograms keep all annual body summaries; outlier display is transparent.
    params=sorted(annual.parameter.unique());ncols=3;nrows=math.ceil(len(params)/ncols)
    fig,axes=plt.subplots(nrows,ncols,figsize=(12,3.3*nrows))
    for ax,p in zip(axes.flat,params):
        g=annual[annual.parameter.eq(p)]
        # Display mass concentrations in mg/L (analysis itself stays in g/L).
        scale,ulabel=(1000.,"mg/L") if g.unit.iloc[0]=="g/L" else (1.,g.unit.iloc[0])
        vals=g.value.to_numpy()*scale
        logged=np.max(vals)>100*max(np.median(vals[vals>0]) if np.any(vals>0) else 1,1e-12)
        tr=(lambda v: np.log10(1+v)) if logged else (lambda v: v)
        if logged:
            # Transform both data and labels explicitly rather than clipped bins.
            ax.hist(tr(vals),bins=35,color="#436b84",alpha=.85)
            ax.set_xlabel(f"log10(1 + x / 1 {ulabel}) [–]")
        else:
            ax.hist(vals,bins=35,color="#436b84",alpha=.85);ax.set_xlabel(f"Value [{ulabel}]")

        e=g[g.eligible_local]
        share=np.mean(e.adverse_deviation>1e-12*e.threshold.abs()) if len(e) else np.nan
        for typ,ls in [("River","--"),("Lake",":")]:
            th=e[e.Type.eq(typ)].threshold
            if len(th):
                ax.axvline(tr(th.iloc[0]*scale),color="#D55E00",ls=ls,lw=1.6,label=f"{typ} threshold")
        if len(e): ax.legend(frameon=False,fontsize=9)
        ax.set_title(f"{PLOT_PARAM.get(p,p)}: n={len(g)}, {g.VOR.nunique()} bodies\nnon-compliant {100*share:.0f}%",fontsize=12)
        ax.set_ylabel("Body-years");ax.grid(alpha=.15)
    for ax in list(axes.flat)[len(params):]: ax.axis("off")
    fig.tight_layout();save_figure(fig,root/"parameter_distributions")
    if not ranking_df.empty:
        fig,ax=plt.subplots(figsize=(9,7))
        m=ranking_df[ranking_df.eligible & ranking_df.Latitude.notna()]
        scatter=ax.scatter(m.Longitude,m.Latitude,c=m.score,cmap="cividis",s=35,edgecolor="white",linewidth=.4)
        crit=m[m.critical]
        ax.scatter(crit.Longitude,crit.Latitude,s=120,facecolor="none",edgecolor="#c0392b",linewidth=1.6,
                   label=f"Critical: score ≥ {CFG.critical_quantile:.0%} quantile (n={len(crit)})")
        ax.legend(frameon=False,fontsize=11,loc="lower left")
        for _,r in m.nsmallest(5,"rank").iterrows():
            ax.annotate(r.VOR,(r.Longitude,r.Latitude),xytext=(4,5),textcoords="offset points",fontsize=13,fontweight="bold")
        ax.set(xlabel="Longitude [degrees E]",ylabel="Latitude [degrees N]",title="Persistence score S; five highest-ranked labelled")
        fig.colorbar(scatter,ax=ax,label="Mean longest adverse annual run / observed years")
        ax.set_aspect(1/np.cos(np.deg2rad(mid)));ax.grid(alpha=.2)
        save_figure(fig,root/"critical_water_bodies")

    for family,table in [("local_associations",rawcorr),("hvg_associations",partialcorr)]:
        for p,g in table.groupby("parameter",sort=True):
            matrix=g.pivot(index="hydro_metric",columns="outcome",values="rho")
            if matrix.empty or not matrix.notna().any().any():continue
            qs=g.pivot(index="hydro_metric",columns="outcome",values="q_BH").reindex(index=matrix.index,columns=matrix.columns)
            ns=g.pivot(index="hydro_metric",columns="outcome",values="n").reindex(index=matrix.index,columns=matrix.columns)
            annotations=np.empty(matrix.shape,dtype=object)
            for i in range(matrix.shape[0]):
                for j in range(matrix.shape[1]):
                    rho=matrix.iloc[i,j]
                    annotations[i,j]="" if pd.isna(rho) else f"{rho:.2f}{'*' if qs.iloc[i,j]<.05 else ''}\nn={int(ns.iloc[i,j])}"
            fig,ax=plt.subplots(figsize=(max(10,1.6*matrix.shape[1]),max(5,1.1*matrix.shape[0])))
            sns.heatmap(matrix,ax=ax,cmap="RdBu_r",vmin=-1,vmax=1,annot=annotations,fmt="",
                        annot_kws={"fontsize":11},linewidths=.4,cbar_kws={"label":"Rank association"})
            kind="Controlled" if g.adjusted.iloc[0] else "Unadjusted"
            ax.set_title(f"{kind} associations: {p}\n* BH q < 0.05 using sub-basin cluster tests; unavailable effects blank")
            ax.set_xlabel("Annual local / temporal HVG metric");ax.set_ylabel("Hydrological node metric")
            ax.tick_params(axis="x",rotation=45);fig.tight_layout()
            save_figure(fig,root/f"{family}_{safe_name(p)}")


# ======================================================================

# ======================================================================

def mantel_haenszel(strata):
    """Mantel-Haenszel common OR, Robins-Breslow-Greenland 95% CI, CMH test.

    strata: iterable of (a,b,c,d) = (exposed&bad, exposed&good,
    unexposed&bad, unexposed&good). Strata with n<2 are ignored.
    """
    R = S = 0.0
    var_terms = [0.0, 0.0, 0.0]
    sum_a = sum_ea = sum_va = 0.0
    informative = 0
    for a, b, c, d in strata:
        n = a + b + c + d
        if n < 2:
            continue
        r_i, s_i = a * d / n, b * c / n
        p_i, q_i = (a + d) / n, (b + c) / n
        R += r_i; S += s_i
        var_terms[0] += p_i * r_i
        var_terms[1] += p_i * s_i + q_i * r_i
        var_terms[2] += q_i * s_i
        row1, row0, col1, col0 = a + b, c + d, a + c, b + d
        if row1 and row0 and col1 and col0:
            informative += 1
            sum_a += a
            sum_ea += row1 * col1 / n
            sum_va += row1 * row0 * col1 * col0 / (n * n * (n - 1))
    out = {"OR_MH": np.nan, "OR_ci_low": np.nan, "OR_ci_high": np.nan,
           "CMH_chi2": np.nan, "p_CMH": np.nan, "n_informative_strata": informative}
    if R > 0 and S > 0:
        orv = R / S
        var = var_terms[0] / (2 * R * R) + var_terms[1] / (2 * R * S) + var_terms[2] / (2 * S * S)
        se = math.sqrt(var)
        out["log_OR_se"] = se
        out.update(OR_MH=orv, OR_ci_low=math.exp(math.log(orv) - 1.96 * se),
                   OR_ci_high=math.exp(math.log(orv) + 1.96 * se))
    if sum_va > 0:
        chi2 = (abs(sum_a - sum_ea) - 0.5) ** 2 / sum_va
        out.update(CMH_chi2=chi2, p_CMH=float(stats.chi2.sf(chi2, 1)))
    return out


def placebo_upstream_events(annual, graph):
    """Negative control: same design with hydrologically UNCONNECTED pseudo-upstream bodies.

    For every body v with k direct predecessors, k bodies are drawn (seeded)
    among those with no directed path to or from v, within the same RBMP
    sub-basin where possible. Their annual adversity cannot propagate to v
    through the recipient network; an association therefore measures shared
    regional / common-year variation (weather, monitoring campaigns).
    """
    eligible = annual[annual.eligible_local].copy()
    eligible["good"] = eligible.adverse_deviation <= 1e-12*eligible.threshold.abs().clip(lower=1e-12)
    lookup = {(v,p):g.set_index("year") for (v,p),g in eligible.groupby(["VOR","parameter"])}
    basin = eligible.groupby("VOR").basin.first().to_dict()
    bodies_by_param = eligible.groupby("parameter").VOR.unique().to_dict()
    rows = []
    for (v,p),g in sorted(lookup.items()):
        k = graph.in_degree(v) if v in graph else 0
        if not k:
            continue
        connected = nx.ancestors(graph,v) | nx.descendants(graph,v) | {v}
        pool = [u for u in bodies_by_param[p] if u not in connected]
        same = [u for u in pool if basin.get(u) == basin.get(v)]
        pool = same if len(same) >= k else pool
        if len(pool) < k:
            continue
        rng = np.random.default_rng(CFG.seed+int(hashlib.sha256(f"placebo|{v}|{p}".encode()).hexdigest()[:8],16))
        pseudo = sorted(rng.choice(pool, size=k, replace=False))
        for year,r in g.iterrows():
            states = [bool(lookup[(u,p)].loc[year,"good"]) for u in pseudo if year in lookup[(u,p)].index]
            if states:
                rows.append({"VOR":v,"parameter":p,"year":int(year),"downstream_good":bool(r.good),
                             "observed_upstream_adverse":not all(states)})
    return pd.DataFrame(rows)


def topological_resilience_tests(upstream, events, placebo=None):
    """Operational topological resilience TR(v,p) and its network-level test.

    TR(v,p) = P(downstream compliant | >=1 observed direct upstream body adverse).
    R_clean(v,p) = P(compliant | all observed direct upstream bodies compliant).
    Buffering contrast = TR - R_clean (<=0 expected if upstream adversity matters).
    Per parameter, a Mantel-Haenszel test stratified by water body asks whether
    downstream non-compliance is more frequent in years with upstream adversity.
    Stratification by body removes all time-invariant body-level confounding
    (typology, catchment); common-year (weather) shocks are NOT removed.
    """
    if upstream.empty or events.empty:
        return pd.DataFrame(), pd.DataFrame()
    up = upstream.copy()
    up["buffering_contrast"] = up.conditional_compliance - up.compliance_when_all_upstream_observed_clean
    rows = []
    for p, g in events.groupby("parameter", sort=True):
        strata = []
        for v, h in g.groupby("VOR"):
            e = h.observed_upstream_adverse.to_numpy(bool)
            bad = ~h.downstream_good.to_numpy(bool)
            strata.append((int((e & bad).sum()), int((e & ~bad).sum()),
                           int((~e & bad).sum()), int((~e & ~bad).sum())))
        mh = mantel_haenszel(strata)
        if placebo is not None and len(placebo):
            pstrata = []
            for v, h in placebo[placebo.parameter.eq(p)].groupby("VOR"):
                e = h.observed_upstream_adverse.to_numpy(bool); bad = ~h.downstream_good.to_numpy(bool)
                pstrata.append((int((e & bad).sum()), int((e & ~bad).sum()), int((~e & bad).sum()), int((~e & ~bad).sum())))
            pm = mantel_haenszel(pstrata)
            mh.update(OR_MH_placebo=pm["OR_MH"], OR_placebo_ci_low=pm["OR_ci_low"],
                      OR_placebo_ci_high=pm["OR_ci_high"], p_CMH_placebo=pm["p_CMH"])
            # Ratio of odds ratios: network-specific association beyond shared variation.
            if np.isfinite(mh["OR_MH"]) and np.isfinite(pm["OR_MH"]) and pm["OR_MH"] > 0:
                mh["ROR_network_vs_placebo"] = mh["OR_MH"]/pm["OR_MH"]
                # Independence approximation: covariance between the two log odds
                # ratios is omitted; conservativeness is not established by this formula.
                se = math.sqrt(mh.get("log_OR_se", np.nan)**2 + pm.get("log_OR_se", np.nan)**2)
                if np.isfinite(se) and se > 0:
                    z = math.log(mh["ROR_network_vs_placebo"])/se
                    mh.update(ROR_ci_low=math.exp(math.log(mh["ROR_network_vs_placebo"])-1.96*se),
                              ROR_ci_high=math.exp(math.log(mh["ROR_network_vs_placebo"])+1.96*se),
                              p_ROR=float(2*stats.norm.sf(abs(z))))
        u = up[up.parameter.eq(p)]
        ua = u[u.adequate_exposed_years]
        rows.append({"parameter": p, "n_bodies_with_upstream_data": int(u.VOR.nunique()),
                     "n_bodies_adequate_exposure": int(ua.VOR.nunique()),
                     "matched_body_years": int(len(g)),
                     "exposed_body_years": int(g.observed_upstream_adverse.sum()),
                     "median_TR": float(ua.conditional_compliance.median()) if len(ua) else np.nan,
                     "TR_q25": float(ua.conditional_compliance.quantile(.25)) if len(ua) else np.nan,
                     "TR_q75": float(ua.conditional_compliance.quantile(.75)) if len(ua) else np.nan,
                     "median_compliance_clean_upstream": float(ua.compliance_when_all_upstream_observed_clean.median()) if ua.compliance_when_all_upstream_observed_clean.notna().any() else np.nan,
                     "median_buffering_contrast": float(ua.buffering_contrast.median()) if ua.buffering_contrast.notna().any() else np.nan,
                     **mh})
    table = pd.DataFrame(rows)
    table["q_BH"] = np.nan
    ok = np.isfinite(table.p_CMH)
    if ok.any():
        table.loc[ok, "q_BH"] = stats.false_discovery_control(table.loc[ok, "p_CMH"].to_numpy(), method="bh")
    if "p_ROR" in table:
        table["q_BH_ROR"] = np.nan
        ok = np.isfinite(table.p_ROR)
        if ok.any():
            table.loc[ok, "q_BH_ROR"] = stats.false_discovery_control(table.loc[ok, "p_ROR"].to_numpy(), method="bh")
    return table, up


def quality_profile_matrix(annual):
    """Body x parameter matrix of mean oriented relative deviation s(x/T-1).

    Same preparation as the quality-profile t-SNE: >=3 parameters per body,
    parameters observed in >=10 bodies, parameter-median imputation, z-scaling.
    """
    profile = annual[annual.eligible_local].copy()
    profile["oriented_relative"] = profile.direction*(profile.value/profile.threshold-1)
    pivot = profile.groupby(["VOR","parameter"]).oriented_relative.mean().unstack()
    pivot = pivot[pivot.notna().sum(axis=1)>=3]
    pivot = pivot[pivot.columns[pivot.notna().sum()>=10]]
    imputed_fraction = float(pivot.isna().to_numpy().mean()) if pivot.size else np.nan
    features = StandardScaler().fit_transform(pivot.fillna(pivot.median()).to_numpy()) if pivot.size else np.empty((0,0))
    return pivot, features, imputed_fraction


def permanova(features, labels, permutations, seed):
    """One-way PERMANOVA on Euclidean distances (Anderson 2001) + PERMDISP.

    For Euclidean distances the pseudo-F equals the centroid sums-of-squares
    ratio. PERMDISP (Anderson 2006): ANOVA F of distances to group centroids.
    """
    x = np.asarray(features, float)
    lab = np.asarray(labels)
    groups = np.unique(lab)
    n, a = len(lab), len(groups)
    def stat(lv):
        grand = x.mean(axis=0)
        sst = float(((x-grand)**2).sum())
        ssw, disp = 0.0, []
        for gname in groups:
            m = lv == gname
            cent = x[m].mean(axis=0)
            ssw += float(((x[m]-cent)**2).sum())
            disp.append(np.sqrt(((x[m]-cent)**2).sum(axis=1)))
        ssa = sst-ssw
        f = (ssa/(a-1))/(ssw/(n-a))
        return f, ssa/sst
    def dispersion_f(lv):
        dist, codes = [], []
        for k, gname in enumerate(groups):
            m = lv == gname
            cent = x[m].mean(axis=0)
            dist.append(np.sqrt(((x[m]-cent)**2).sum(axis=1))); codes.append(np.full(m.sum(), k))
        return stats.f_oneway(*dist).statistic
    f_obs, r2 = stat(lab)
    d_obs = dispersion_f(lab)
    rng = np.random.default_rng(seed)
    f_perm, d_perm = [], []
    for _ in range(permutations):
        lp = rng.permutation(lab)
        f_perm.append(stat(lp)[0]); d_perm.append(dispersion_f(lp))
    return {"n": n, "groups": a, "pseudo_F": f_obs, "R2": r2,
            "p_perm": (1+np.sum(np.array(f_perm) >= f_obs))/(permutations+1),
            "PERMDISP_F": d_obs,
            "p_PERMDISP": (1+np.sum(np.array(d_perm) >= d_obs))/(permutations+1)}


def typology_tests(annual, meta):
    """Formal test of typological separation of physicochemical profiles."""
    pivot, features, imputed = quality_profile_matrix(annual)
    rows, per_param = [], []
    if pivot.empty:
        return pd.DataFrame(), pd.DataFrame(), imputed
    m = meta.set_index("VOR").reindex(pivot.index)
    for grouping in ["Type", "typology_code", "geochemistry", "altitude_class"]:
        lab = m[grouping].fillna("Unknown").astype(str)
        counts = lab.value_counts()
        keep_groups = counts[counts >= CFG.min_group_size_typology].index
        keep = lab.isin(keep_groups) & lab.ne("Unknown")
        if keep_groups.size < 2 or keep.sum() < 2*CFG.min_group_size_typology:
            rows.append({"grouping": grouping, "status": "insufficient_groups"}); continue
        res = permanova(features[keep.to_numpy()], lab[keep].to_numpy(), CFG.permutations, CFG.seed)
        rows.append({"grouping": grouping, "status": "ok", **res})
    table = pd.DataFrame(rows)
    if "p_perm" in table:
        ok = table.p_perm.notna()
        table.loc[ok, "q_BH"] = stats.false_discovery_control(table.loc[ok, "p_perm"].to_numpy(), method="bh")
    # Univariate lake vs river contrasts on the non-imputed body means.
    types = m.Type
    for p in pivot.columns:
        lake = pivot.loc[types.eq("Lake"), p].dropna(); river = pivot.loc[types.eq("River"), p].dropna()
        if len(lake) >= 5 and len(river) >= 5:
            test = stats.mannwhitneyu(lake, river, alternative="two-sided")
            per_param.append({"parameter": p, "n_lake": len(lake), "n_river": len(river),
                              "median_lake": lake.median(), "median_river": river.median(),
                              "rank_biserial": 1-2*test.statistic/(len(lake)*len(river)),
                              "p_MWU": test.pvalue})
    per_param = pd.DataFrame(per_param)
    if len(per_param):
        per_param["q_BH"] = stats.false_discovery_control(per_param.p_MWU.to_numpy(), method="bh")
    return table, per_param, imputed


def aggregation_summary(agg):
    """Concordance of HVG metrics under alternative annual aggregation schemes."""
    rows = []
    if agg.empty:
        return pd.DataFrame()
    for (p, scheme, metric), g in agg.groupby(["parameter", "scheme", "metric"], sort=True):
        g = g.dropna(subset=["baseline", "perturbed"])
        rho = (stats.spearmanr(g.baseline, g.perturbed).statistic
               if len(g) > 2 and g.baseline.nunique() > 1 and g.perturbed.nunique() > 1 else np.nan)
        rows.append({"parameter": p, "scheme": scheme, "metric": metric, "n_series": len(g),
                     "rho_baseline_vs_alternative": rho,
                     "median_absolute_change": float((g.perturbed-g.baseline).abs().median())})
    return pd.DataFrame(rows)


def hvg_pca(hvg):
    """Dimensionality of the HVG metric set: PCA on z-scored metrics."""
    from sklearn.decomposition import PCA
    cols = [c for c in HVG_FULL if c in hvg and hvg[c].nunique() > 1]
    x = hvg[cols].replace([np.inf, -np.inf], np.nan).dropna()
    if len(x) < 10:
        return pd.DataFrame(), pd.DataFrame()
    z = StandardScaler().fit_transform(x)
    pca = PCA(random_state=CFG.seed).fit(z)
    var = pd.DataFrame({"component": [f"PC{i+1}" for i in range(len(cols))],
                        "explained_variance_ratio": pca.explained_variance_ratio_,
                        "cumulative": np.cumsum(pca.explained_variance_ratio_)})
    load = pd.DataFrame(pca.components_[:4].T, index=cols, columns=[f"PC{i+1}" for i in range(min(4, len(cols)))])
    return var, load.reset_index(names="metric")


# ---------------------------- figures ---------------------------------

PLOT_PARAM = {"Phytoplankton (biomass  - chlorophyll-a)": "Chlorophyll-a", "Turbidity (NTU)": "Turbidity",
              "Orthophosphate (PO4-P)": "PO4-P", "Ammonium (NH4-N)": "NH4-N", "Nitrate (NO3-N)": "NO3-N",
              "Nitrite (NO2-N)": "NO2-N", "BOD (5)": "BOD5", "COD (Cr)": "CODcr"}
PLOT_METRIC = {"N_good": "N+", "N_bad": "N−", "Amax_bad_relative": "A−max/T", "Amax_good_relative": "A+max/T",
               "Lmax_bad_years": "L−max", "Lmax_good_years": "L+max", "ARL_bad_years": "ARL−",
               "ARL_good_years": "ARL+", "tR_first_observed_years": "tR", "noncompliance_fraction": "N−/n"}


def association_forest(table, stem, title, top=None):
    """Strongest controlled associations with conditional 95% intervals."""
    top = top or CFG.top_associations
    t = table[table.rho.notna()].copy()
    if t.empty:
        return
    t["abs"] = t.rho.abs()
    t = t.sort_values("abs", ascending=False).head(top).iloc[::-1]
    labels = [f"{PLOT_PARAM.get(r.parameter, r.parameter)}:  {r.hydro_metric} ~ {PLOT_METRIC.get(r.outcome, r.outcome)}"
              for r in t.itertuples()]
    fig, ax = plt.subplots(figsize=(9.5, max(4, .40*len(t)+1.4)))
    y = np.arange(len(t))
    sig = t.q_BH.lt(.05).to_numpy()
    lo = (t.rho-t.ci_low).clip(lower=0).fillna(0); hi = (t.ci_high-t.rho).clip(lower=0).fillna(0)
    ax.errorbar(t.rho, y, xerr=[lo, hi], fmt="none", ecolor="#7a7a7a", lw=1.2, capsize=3)
    ax.scatter(t.rho[~sig], y[~sig], s=55, facecolor="white", edgecolor="#386b86", lw=1.5, zorder=3, label="BH q ≥ 0.05 or not tested")
    ax.scatter(t.rho[sig], y[sig], s=60, color="#c0392b", zorder=3, label="BH q < 0.05 (sub-basin cluster test)")
    ax.axvline(0, color="black", lw=.8)
    for thr in (-.3, .3):
        ax.axvline(thr, color="#bbbbbb", ls=":", lw=1)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlim(-1, 1)
    ax.set_xlabel("Controlled rank association ρ [–]; bars: 95% interval", fontsize=13)
    ax.set_title(title, fontsize=13)
    ax.legend(frameon=False, fontsize=11, loc="upper center", bbox_to_anchor=(.5, -.09), ncol=2)
    ax.grid(axis="x", alpha=.2)
    fig.tight_layout(); save_figure(fig, stem)


def topological_figure(up, mh, stem):
    if up.empty or mh.empty:
        return
    ua = up[up.adequate_exposed_years]
    params = sorted(mh.parameter)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, max(5, .5*len(params)+2.2)),
                                 gridspec_kw={"width_ratios": [1.3, 1]}, sharey=True)
    y = np.arange(len(params))
    for i, p in enumerate(params):
        g = ua[ua.parameter.eq(p)]
        if len(g):
            kw = dict(positions=[i], widths=.55, patch_artist=True, boxprops={"facecolor": "#c5deeb"},
                      medianprops={"color": "#0b3d5c"})
            try:
                a1.boxplot(g.conditional_compliance.dropna(), orientation="horizontal", **kw)
            except TypeError:  # matplotlib < 3.10
                a1.boxplot(g.conditional_compliance.dropna(), vert=False, **kw)
            a1.scatter(g.compliance_when_all_upstream_observed_clean, np.full(len(g), i)+.3, s=10,
                       color="#D55E00", alpha=.5)
    a1.scatter([], [], s=14, color="#D55E00", label="R_clean: compliance when upstream compliant")
    a1.set_yticks(y); a1.set_yticklabels([PLOT_PARAM.get(q, q) for q in params], fontsize=12)
    a1.set_xlabel("TR = P(compliant | upstream non-compliant) [–]", fontsize=12)
    a1.set_xlim(-.02, 1.02); a1.grid(axis="x", alpha=.2)
    a1.legend(frameon=False, fontsize=11, loc="upper center", bbox_to_anchor=(.5, -.13))
    m = mh.set_index("parameter").reindex(params)
    ok = m.OR_MH.notna().to_numpy()
    lo = (m.OR_MH-m.OR_ci_low).clip(lower=0).fillna(0); hi = (m.OR_ci_high-m.OR_MH).clip(lower=0).fillna(0)
    a2.errorbar(m.OR_MH[ok], y[ok], xerr=[lo[ok], hi[ok]], fmt="o", color="#0b3d5c", capsize=3)
    sig = (m.q_BH < .05).to_numpy()
    a2.scatter(m.OR_MH[sig], y[sig], s=70, color="#c0392b", zorder=3, label="Upstream water bodies, q < 0.05")
    if "OR_MH_placebo" in m:
        okp = m.OR_MH_placebo.notna().to_numpy()
        lp = (m.OR_MH_placebo-m.OR_placebo_ci_low).clip(lower=0).fillna(0)
        hp = (m.OR_placebo_ci_high-m.OR_MH_placebo).clip(lower=0).fillna(0)
        a2.errorbar(m.OR_MH_placebo[okp], y[okp]-.28, xerr=[lp[okp], hp[okp]], fmt="s", color="#999999",
                    ms=5, capsize=2, label="Placebo: unconnected water bodies")
    a2.axvline(1, color="black", lw=.8); a2.set_xscale("log")
    a2.set_xlabel("Mantel–Haenszel odds ratio [–] (log scale)", fontsize=12)
    a2.legend(frameon=False, fontsize=11, loc="upper center", bbox_to_anchor=(.5, -.13), ncol=1)
    a2.grid(axis="x", alpha=.2)
    fig.tight_layout(); save_figure(fig, stem)


def hvg_sensitivity_figure(summary, stem):
    if summary.empty:
        return
    med = summary.groupby(["method", "deletion_fraction", "metric"], as_index=False).rho_baseline_vs_median_perturbed.median()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), sharey=True)
    for ax, method in zip(axes, ["random", "contiguous_block"]):
        g = med[med.method.eq(method)]
        for metric, h in g.groupby("metric"):
            ax.plot(100*h.deletion_fraction, h.rho_baseline_vs_median_perturbed, "o-", lw=1.8, label=metric)
        ax.set_title(f"{'Random' if method=='random' else 'Contiguous-block'} deletion", fontsize=13)
        ax.set_xlabel("Deleted share of annual values [%]", fontsize=12); ax.grid(alpha=.2)
        ax.set_ylim(min(0, med.rho_baseline_vs_median_perturbed.min()-.05), 1.02)
    axes[0].set_ylabel("Spearman ρ, baseline vs perturbed [–]", fontsize=12)
    axes[1].legend(frameon=False, fontsize=11, title="HVG metric")
    fig.tight_layout(); save_figure(fig, stem)


# ---------------------------- LaTeX export ----------------------------

def tex_escape(text):
    text = "" if text is None or (isinstance(text, float) and np.isnan(text)) else str(text)
    for a, b in [("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("$", r"\$"), ("#", r"\#"),
                 ("_", r"\_"), ("{", r"\{"), ("}", r"\}"), ("~", r"\textasciitilde{}"), ("^", r"\^{}")]:
        text = text.replace(a, b)
    return text


def fmt_int(x):
    try:
        if x is None or (isinstance(x, float) and not np.isfinite(x)):
            return "--"
        return f"{int(round(float(x))):,}".replace(",", "{,}")
    except (TypeError, ValueError):
        return "--"


def fmt_num(x, digits=2):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "--"
    if not np.isfinite(x):
        return "--"
    return f"{x:.{digits}f}"


def fmt_p(x):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "--"
    if not np.isfinite(x):
        return "--"
    if x < 0.001:
        return "$<$0.001"
    return f"{x:.3f}"


def tabular(path, colspec, header, rows, note=None):
    lines = [r"\begin{tabular}{" + colspec + "}", r"\toprule", " & ".join(header) + r" \\", r"\midrule"]
    lines += [" & ".join(r) + r" \\" for r in rows] or [r"\multicolumn{" + str(len(header)) + r"}{c}{no eligible rows} \\"]
    lines += [r"\bottomrule", r"\end{tabular}"]
    if note:
        lines.append(r"\par\vspace{2pt}{\footnotesize " + note + "}")
    Path(path).write_text("% Generated by surface_surface_water_resilience.py - do not edit by hand\n" + "\n".join(lines) + "\n", encoding="utf-8")


PARAM_SHORT = {"Phytoplankton (biomass  - chlorophyll-a)": "Chlorophyll-a", "Orthophosphate (PO4-P)": "PO$_4$-P",
               "Ammonium (NH4-N)": "NH$_4$-N", "Nitrate (NO3-N)": "NO$_3$-N", "Nitrite (NO2-N)": "NO$_2$-N",
               "BOD (5)": "BOD$_5$", "COD (Cr)": "COD$_\\mathrm{Cr}$", "Turbidity (NTU)": "Turbidity"}


def pname(p):
    return PARAM_SHORT.get(p, tex_escape(p))


METRIC_TEX = {"N_good": "$N^+$", "N_bad": "$N^-$", "Amax_bad_relative": "$A^-_{\\max}/T$",
              "Amax_good_relative": "$A^+_{\\max}/T$", "Lmax_bad_years": "$L^-_{\\max}$",
              "Lmax_good_years": "$L^+_{\\max}$", "ARL_bad_years": "$ARL^-$", "ARL_good_years": "$ARL^+$",
              "tR_first_observed_years": "$t_R$", "noncompliance_fraction": "$N^-/n$"}


def mname(m):
    return METRIC_TEX.get(m, tex_escape(m))


def export_generated(gen, figs, ctx):
    """Export numerical macros, tables and figures for LaTeX documents."""
    import shutil
    gen = Path(gen); gen.mkdir(parents=True, exist_ok=True)
    macros = {}
    m = ctx["manifest"]; audit = m["data_audit"]
    meta = ctx["meta"]; local = ctx["local"]; hvg = ctx["hvg"]; ranks = ctx["ranks"]
    reg = meta[meta.registered].Type.value_counts(); obs = meta[meta.observed].Type.value_counts()
    macros.update(
        resRawRows=fmt_int(audit["raw_rows"]), resRawVOR=fmt_int(m["raw_VOR_count"]),
        resRawSites=fmt_int(m["raw_monitoring_KTJ_count"]),
        resRegRivers=fmt_int(reg.get("River", 0)), resRegLakes=fmt_int(reg.get("Lake", 0)),
        resObsRivers=fmt_int(obs.get("River", 0)), resObsLakes=fmt_int(obs.get("Lake", 0)),
        resObsOther=fmt_int(obs.get("Other", 0)),
        resNetNodes=fmt_int(m["network_nodes"]), resNetEdges=fmt_int(m["network_edges"]),
        resNetDAG="acyclic" if m["network_DAG"] else "not acyclic",
        resNetComponents=fmt_int(meta.weak_component.nunique()),
        resSelectedRows=fmt_int(audit["selected_rows_before_dedup"]),
        resSiteVersionRows=fmt_int(audit["site_version_rows_after_copy_removal"]),
        resSiteYearRows=fmt_int(audit["site_parameter_year_rows"]),
        resConflicts=fmt_int(audit["conflicting_site_parameter_years"]),
        resBodyYearRows=fmt_int(audit["annual_body_parameter_rows"]),
        resBodyParamComb=fmt_int(audit["annual_body_parameter_combinations"]),
        resLocalRows=fmt_int(len(local)),
        resLocalEligible=fmt_int((local.n_years >= CFG.min_years_local).sum()),
        resLocalBodies=fmt_int(local.VOR.nunique()),
        resHVGRows=fmt_int(len(hvg)), resHVGExcluded=fmt_int(len(ctx["hvg_rejected"])),
        resParams=fmt_int(ctx["annual"].parameter.nunique()),
        resMedianYears=fmt_num(ctx["coverage"].n_years.median(), 0),
        resYearMin=str(CFG.start_year), resYearMax=str(CFG.end_year),
        resMinYears=str(CFG.min_years_local), resMinParams=str(CFG.min_parameters_ranking),
        resCriticalQ=f"{CFG.critical_quantile*100:.0f}", resRedundancyRho=fmt_num(CFG.redundancy_rho, 2),
        resBootReps=fmt_int(CFG.bootstrap_reps), resPerms=fmt_int(CFG.permutations),
        resSeed=str(CFG.seed), resTsnePerplexity=fmt_num(CFG.tsne_perplexity, 0),
        resTsneIter=fmt_int(CFG.tsne_iterations), resTsneLR=fmt_num(CFG.tsne_learning_rate, 0),
        resMinCorrN=str(CFG.correlation_min_n), resMinBasins=str(CFG.correlation_min_basins),
        resMinExposed=str(CFG.min_upstream_exposed_years), resSensReps=fmt_int(CFG.sensitivity_reps),
        resSensMinYears=str(CFG.sensitivity_min_years),
        resHydroSelected=", ".join(map(tex_escape,ctx["selected_hydro"])), resHVGSelected=", ".join(map(tex_escape,ctx["selected_hvg"])),
        resBasins=fmt_int(meta[meta.registered].basin.nunique()),
    )
    for k, v in m["versions"].items():
        macros["resV" + re.sub(r"[^A-Za-z]", "", k.title())] = tex_escape(v)
    macros["resQuickRun"] = "QUICK TEST RUN -- not for publication" if ctx.get("quick") else "full configuration"
    # Ranking
    el = ranks[ranks.eligible] if len(ranks) else ranks
    macros["resRankEligible"] = fmt_int(len(el))
    macros["resCriticalN"] = fmt_int(el.critical.sum()) if len(el) else "--"
    macros["resCriticalCutoff"] = fmt_num(el.critical_cutoff.iloc[0], 3) if len(el) else "--"
    top5 = el.nsmallest(5, "rank") if len(el) else el
    macros["resTopFive"] = ", ".join(f"{r.VOR} ({tex_escape(r['name'])})" for _, r in top5.iterrows()) or "--"
    agr = ctx.get("rank_agreement", pd.DataFrame())
    if len(agr):
        yo = agr[agr.scenario.str.startswith("omit_year_")]
        other = agr[~agr.scenario.str.startswith("omit_year_")]
        macros.update(resRankRhoMin=fmt_num(agr.spearman_score_vs_baseline.min()),
                      resRankRhoYearMedian=fmt_num(yo.spearman_score_vs_baseline.median()),
                      resRankJaccardMin=fmt_num(agr.jaccard_critical_sets.min()),
                      resRankJaccardMedian=fmt_num(agr.jaccard_critical_sets.median()),
                      resRankJaccardYearMedian=fmt_num(yo.jaccard_critical_sets.median()),
                      resRankTopFiveMin=fmt_num(agr.top5_overlap.min()*100, 0))
    # Associations
    for key, table in [("Local", ctx["corr_local_adj"]), ("HVG", ctx["corr_hvg_adj"]),
                       ("LocalRaw", ctx["corr_local_raw"]), ("HVGRaw", ctx["corr_hvg_raw"])]:
        tested = table.p_cluster.notna() if len(table) else pd.Series(dtype=bool)
        macros[f"res{key}Pairs"] = fmt_int(table.rho.notna().sum()) if len(table) else "0"
        macros[f"res{key}Tests"] = fmt_int(tested.sum())
        macros[f"res{key}Sig"] = fmt_int(table.q_BH.lt(.05).sum()) if len(table) else "0"
        if len(table) and table.rho.notna().any():
            sig = table[table.q_BH.lt(.05)]
            macros[f"res{key}MaxAbsRho"] = fmt_num(table.rho.abs().max())
            macros[f"res{key}MedAbsRhoSig"] = fmt_num(sig.rho.abs().median()) if len(sig) else "--"
            best = (sig if len(sig) else table[table.rho.notna()]).assign(a=lambda d: d.rho.abs()).sort_values("a").iloc[-1]
            macros[f"res{key}TopPair"] = f"{pname(best.parameter)}: {tex_escape(best.hydro_metric)}--{mname(best.outcome)}"
            macros[f"res{key}TopRho"] = fmt_num(best.rho)
            macros[f"res{key}TopCI"] = f"$[{fmt_num(best.ci_low)},\\ {fmt_num(best.ci_high)}]$"
            macros[f"res{key}TopQ"] = fmt_p(best.q_BH)
            macros[f"res{key}TopN"] = fmt_int(best.n)
    # Typology
    typ, typp, imputed = ctx["typology"], ctx["typology_params"], ctx["typology_imputed"]
    macros["resTsneImputed"] = fmt_num(100*imputed, 1) if imputed == imputed else "--"
    if len(typ) and "p_perm" in typ:
        t = typ.set_index("grouping")
        if "Type" in t.index and t.loc["Type", "status"] == "ok":
            r = t.loc["Type"]
            macros.update(resPermN=fmt_int(r.n), resPermF=fmt_num(r.pseudo_F), resPermRsq=fmt_num(r.R2, 3),
                          resPermP=fmt_p(r.p_perm), resPermDispP=fmt_p(r.p_PERMDISP))
    if len(typp):
        macros.update(resMWUTests=fmt_int(len(typp)), resMWUSig=fmt_int(typp.q_BH.lt(.05).sum()))
    tsne_meta = m.get("tsne", {}).get("corrected", {})
    macros["resTsneN"] = fmt_int(tsne_meta.get("n_bodies", np.nan))
    macros["resTsneParams"] = fmt_int(len(tsne_meta.get("input_parameters", [])))
    # Topological resilience
    mh, up = ctx["mh"], ctx["upstream_ext"]
    if len(mh):
        macros.update(resTRParams=fmt_int(mh.OR_MH.notna().sum()), resTRSig=fmt_int(mh.q_BH.lt(.05).sum()),
                      resTRPairs=fmt_int(up.adequate_exposed_years.sum()),
                      resTRBodies=fmt_int(up[up.adequate_exposed_years].VOR.nunique()),
                      resTRMedian=fmt_num(up[up.adequate_exposed_years].conditional_compliance.median()),
                      resTRContrastMedian=fmt_num(up[up.adequate_exposed_years].buffering_contrast.median()),
                      resMHORMin=fmt_num(mh.OR_MH.min()), resMHORMax=fmt_num(mh.OR_MH.max()))
        if "OR_MH_placebo" in mh:
            macros.update(resPlaceboORMin=fmt_num(mh.OR_MH_placebo.min()), resPlaceboORMax=fmt_num(mh.OR_MH_placebo.max()),
                          resRORMedian=fmt_num(mh.ROR_network_vs_placebo.median()),
                          resRORAboveOne=fmt_int((mh.ROR_network_vs_placebo > 1).sum()),
                          resRORSig=fmt_int(((mh.ROR_network_vs_placebo > 1) & mh.get("q_BH_ROR", pd.Series(np.nan, index=mh.index)).lt(.05)).sum()))
    # HVG sensitivity / aggregation / PCA
    ss = ctx["sens_summary"]
    if len(ss):
        macros["resSensSeries"] = fmt_int(len(ctx["sens_sample"]))
        for frac, name in [(.10, "Ten"), (.25, "TwentyFive"), (.50, "Fifty")]:
            for method, mn in [("random", "Rand"), ("contiguous_block", "Block")]:
                v = ss[(ss.deletion_fraction.round(2) == frac) & ss.method.eq(method)].rho_baseline_vs_median_perturbed.median()
                macros[f"resSens{mn}{name}"] = fmt_num(v)
    ag = ctx["agg_summary"]
    if len(ag):
        for scheme, nm in [("sample_count_weighted", "Weighted"), ("two_year_bins", "Bins")]:
            macros[f"resAgg{nm}"] = fmt_num(ag[ag.scheme.eq(scheme)].rho_baseline_vs_alternative.median())
    pv = ctx["pca_var"]
    if len(pv):
        macros["resPCAone"] = fmt_num(100*pv.explained_variance_ratio.iloc[0], 0)
        macros["resPCAk"] = fmt_int(int((pv.cumulative < .8).sum()+1))
        macros["resPCAmetrics"] = fmt_int(len(pv))
    lines = ["% Generated by surface_surface_water_resilience.py - do not edit by hand",
             f"% run: {time.strftime('%Y-%m-%d %H:%M')}; configuration: {macros['resQuickRun']}"]
    for k, v in sorted(macros.items()):
        lines.append(rf"\def\{k}{{{v}}}")
    (gen/"analysis_values.tex").write_text("\n".join(lines)+"\n", encoding="utf-8")

    # ---- Tables
    cov = ctx["parameter_counts"]; annual = ctx["annual"]
    nc = annual[annual.eligible_local].assign(bad=lambda d: d.adverse_deviation > 1e-12*d.threshold.abs()).groupby("parameter").bad.mean()
    rows = []
    for _, r in cov.sort_values("parameter").iterrows():
        rows.append([pname(r.parameter), tex_escape(r.unit), fmt_int(r.site_year_records), fmt_int(r.body_year_records),
                     fmt_int(r.body_parameter_combinations), f"{fmt_int(r.median_years)} ({fmt_int(r.min_years)}--{fmt_int(r.max_years)})",
                     fmt_int(r.eligible_local_ranking), fmt_int(r.eligible_HVG), fmt_num(100*nc.get(r.parameter, np.nan), 1)])
    tabular(gen/"table_parameter_coverage.tex", "lllrrrrrr",
            ["Parameter", "Unit", "Site-years", "Body-years", "Bodies", "Years, median (range)",
             r"Local ($\geq$" + str(CFG.min_years_local) + " y)", r"HVG ($\geq$" + str(CFG.min_years_hvg) + " y)", r"Non-compl. [\%]"], rows)
    bal = local[local.VOR.eq(CFG.selected_vor)].sort_values("parameter")
    rows = [[pname(r.parameter), fmt_int(r.n_years), fmt_int(r.N_good), fmt_int(r.N_bad),
             fmt_num(r.Amax_good_relative), fmt_num(r.Amax_bad_relative), fmt_int(r.Lmax_good_years),
             fmt_int(r.Lmax_bad_years), fmt_num(r.ARL_good_years, 1), fmt_num(r.ARL_bad_years, 1),
             fmt_int(r.tR_first_observed_years) if pd.notna(r.tR_first_observed_years) else
             ("n.a." if r.recovery_status == "no_adverse_year" else "cens.")] for _, r in bal.iterrows()]
    tabular(gen/"table_lake_balaton.tex", "lrrrrrrrrrr",
            ["Parameter", "$n$ [y]", "$N^+$ [y]", "$N^-$ [y]", r"$A^+_{\max}/T$ [--]", r"$A^-_{\max}/T$ [--]",
             r"$L^+_{\max}$ [y]", r"$L^-_{\max}$ [y]", "$ARL^+$ [y]", "$ARL^-$ [y]", "$t_R$ [y]"], rows)
    env = ctx.get("rank_envelope", pd.DataFrame())
    envi = env.set_index("VOR") if len(env) else pd.DataFrame()
    rows = []
    for _, r in el.nsmallest(10, "rank").iterrows() if len(el) else []:
        e = envi.loc[r.VOR] if len(envi) and r.VOR in envi.index else None
        rows.append([fmt_int(r["rank"]), tex_escape(r.VOR), tex_escape(r["name"]), tex_escape(r.Type), fmt_num(r.score, 3),
                     fmt_num(r.mean_noncompliance, 2), fmt_int(r.n_parameters),
                     f"{fmt_int(e.rank_min)}--{fmt_int(e.rank_max)}" if e is not None else "--",
                     fmt_num(100*e.critical_year_omission_fraction, 0) if e is not None else "--"])
    tabular(gen/"table_top_ranked_water_bodies.tex", "rllllrrrr",
            ["Rank", "VOR", "Name", "Type", "Score [--]", r"$\overline{N^-/n}$", "Param.", "Rank range$^a$", r"Critical [\%]$^a$"], rows)
    if len(agr):
        yo = agr[agr.scenario.str.startswith("omit_year_")]
        show = agr[~agr.scenario.str.startswith("omit_year_")].copy()
        rows = [[tex_escape(r.scenario.replace("_", " ")), fmt_int(r.n_bodies), fmt_num(r.spearman_score_vs_baseline),
                 fmt_num(r.jaccard_critical_sets), fmt_num(100*r.top5_overlap, 0)] for _, r in show.iterrows()]
        if len(yo):
            rows.append([f"single-year omission ({len(yo)} scenarios), median [min]", fmt_int(yo.n_bodies.median()),
                         f"{fmt_num(yo.spearman_score_vs_baseline.median())} [{fmt_num(yo.spearman_score_vs_baseline.min())}]",
                         f"{fmt_num(yo.jaccard_critical_sets.median())} [{fmt_num(yo.jaccard_critical_sets.min())}]",
                         f"{fmt_num(100*yo.top5_overlap.median(), 0)} [{fmt_num(100*yo.top5_overlap.min(), 0)}]"])
        tabular(gen/"table_ranking_sensitivity.tex", "lrrrr",
                ["Scenario", "Bodies", r"Spearman $\rho$ vs baseline", "Jaccard (critical sets)", r"Top-5 overlap [\%]"], rows)
    rows = []
    for graph_name, sel in [("Hydrological", ctx["hydro_selection"]), ("Temporal (HVG)", ctx["hvg_selection"])]:
        for _, r in sel.iterrows():
            why = ("retained" if r.retained else
                   (f"redundant with {tex_escape(r.with_metric)} ($|\\rho|$={fmt_num(abs(r.rho))})" if r.reason.startswith("redundant") else "constant/unavailable"))
            rows.append([graph_name, tex_escape(r.metric), why])
    tabular(gen/"table_metric_selection.tex", "lll", ["Graph", "Metric (priority order)", "Outcome-blind decision"], rows)
    rows = []
    for fam, table in [("Local", ctx["corr_local_adj"]), ("HVG", ctx["corr_hvg_adj"])]:
        t = table[table.rho.notna()].assign(a=lambda d: d.rho.abs()).sort_values(["a"], ascending=False).head(8)
        for _, r in t.iterrows():
            rows.append([fam, pname(r.parameter), tex_escape(r.hydro_metric), mname(r.outcome), fmt_int(r.n), fmt_int(r.n_basins),
                         fmt_num(r.rho), f"{fmt_num(r.ci_low)}, {fmt_num(r.ci_high)}", fmt_p(r.p_cluster), fmt_p(r.q_BH)])
    tabular(gen/"table_principal_associations.tex", "lllllrrrrr",
            ["Family", "Parameter", "Hydro.", "Outcome", "$n$", "Sub-basins", r"$\rho_\mathrm{ctrl}$", "95\\% interval", "$p_\\mathrm{cluster}$", "$q_\\mathrm{BH}$"], rows)
    rows = []
    if len(typ) and "p_perm" in typ:
        for _, r in typ.iterrows():
            if r.status != "ok":
                rows.append([tex_escape(r.grouping), "--", "--", "--", "--", "--", "insufficient groups"]); continue
            rows.append([tex_escape(r.grouping.replace("_", " ")), fmt_int(r.n), fmt_int(r.groups), fmt_num(r.pseudo_F),
                         fmt_num(r.R2, 3), fmt_p(r.p_perm), fmt_p(r.p_PERMDISP)])
    tabular(gen/"table_typology_tests.tex", "lrrrrrr",
            ["Grouping", "Bodies", "Groups", "Pseudo-$F$", "$R^2$", "$p_\\mathrm{perm}$", "$p_\\mathrm{PERMDISP}$"], rows)
    if len(mh):
        rows = [[pname(r.parameter), fmt_int(r.n_bodies_adequate_exposure), fmt_int(r.exposed_body_years),
                 fmt_num(r.median_TR), fmt_num(r.median_compliance_clean_upstream),
                 f"{fmt_num(r.OR_MH)} ({fmt_num(r.OR_ci_low)}--{fmt_num(r.OR_ci_high)})", fmt_p(r.q_BH),
                 f"{fmt_num(r.get('OR_MH_placebo', np.nan))} ({fmt_num(r.get('OR_placebo_ci_low', np.nan))}--{fmt_num(r.get('OR_placebo_ci_high', np.nan))})",
                 f"{fmt_num(r.get('ROR_network_vs_placebo', np.nan))} ({fmt_num(r.get('ROR_ci_low', np.nan))}--{fmt_num(r.get('ROR_ci_high', np.nan))})",
                 fmt_p(r.get('q_BH_ROR', np.nan))]
                for _, r in mh.iterrows()]
        tabular(gen/"table_topological_resilience.tex", "lrrrrlrllr",
                ["Parameter", "Bodies$^a$", "Exposed body-years", r"Median $TR$", r"Median $R^\mathrm{clean}$",
                 r"$OR_\mathrm{MH}$ (95\% CI)", "$q_\\mathrm{BH}$", r"Placebo $OR_\mathrm{MH}$ (95\% CI)",
                 r"$ROR$ (95\% CI)$^b$", "$q_\\mathrm{BH}$"], rows)
    if len(ss):
        med = ss.groupby(["metric", "method", "deletion_fraction"]).rho_baseline_vs_median_perturbed.median().unstack()
        agm = ag.groupby(["metric", "scheme"]).rho_baseline_vs_alternative.median().unstack() if len(ag) else pd.DataFrame()
        rows = []
        for (metric, method), r in med.iterrows():
            extra = []
            if method == "random":
                extra = [fmt_num(agm.loc[metric, "sample_count_weighted"]) if metric in agm.index and "sample_count_weighted" in agm else "--",
                         fmt_num(agm.loc[metric, "two_year_bins"]) if metric in agm.index and "two_year_bins" in agm else "--"]
            else:
                extra = ["", ""]
            rows.append([tex_escape(metric) if method == "random" else "", "random" if method == "random" else "block",
                         *[fmt_num(r.get(f, np.nan)) for f in CFG.deletion_fractions], *extra])
        tabular(gen/"table_hvg_sensitivity.tex", "llrrrrr",
                ["Metric", "Deletion", *[f"{int(100*f)}\\%" for f in CFG.deletion_fractions], "Weighted", "2-year bins"], rows)

    # ---- Figures with fixed names
    mapping = {"annual_threshold_deviation": "annual_threshold_deviation", "annual_visibility_graph": "annual_visibility_graph",
               "site_annual_comparison": "site_annual_comparison", "monitoring_network": "monitoring_network",
               "parameter_distributions": "parameter_distributions", "quality_profile_tsne": "quality_profile_tsne",
               "critical_water_bodies": "critical_water_bodies", "local_association_forest": "local_association_forest",
               "hvg_association_forest": "hvg_association_forest", "topological_resilience": "topological_resilience",
               "hvg_deletion_sensitivity": "hvg_deletion_sensitivity"}
    copied = []
    for src, dst in mapping.items():
        f = Path(figs)/f"{src}.pdf"
        if f.is_file():
            shutil.copyfile(f, gen/f"{dst}.pdf"); copied.append(dst)
    # Multi-page supplementary heatmap books
    from matplotlib.backends.backend_pdf import PdfPages
    for book, prefix in [("supplementary_local_association_heatmaps", "local_associations_"), ("supplementary_hvg_association_heatmaps", "hvg_associations_")]:
        files = sorted(p for p in Path(figs).glob(prefix+"*.png") if "forest" not in p.name)
        if files:
            with PdfPages(gen/f"{book}.pdf") as pdf:
                for f in files:
                    img = plt.imread(f)
                    fig = plt.figure(figsize=(img.shape[1]/150, img.shape[0]/150))
                    fig.figimage(img, resize=True); pdf.savefig(fig); plt.close(fig)
    return {"macros": len(macros), "figures": copied}


def metric_dictionary():
    rows=[]
    def add(name,graph,formula,unit,note):
        rows.append({"metric":name,"graph_or_level":graph,"definition":formula,"unit":unit,"interpretation":note})
    for name,formula,note in [
        ("DCI","sum_u A_uv","number of direct upstream recipient connections"),
        ("DCO","sum_u A_vu","number of outgoing recipient connections"),
        ("DC","DCI+DCO","all incident directed connections"),
        ("BC","sum_{s!=v!=t} sigma_st(v)/sigma_st divided by (N-1)(N-2)","directed normalized shortest-path betweenness"),
        ("CC_IN","NetworkX WF improved inward closeness","disconnected-network correction; incoming paths"),
        ("CC_OUT","WF improved closeness on reversed graph","outgoing paths"),
        ("HC_IN","sum_{u!=v} 1/d(u,v); unreachable contributes 0","inward unnormalized harmonic centrality"),
        ("HC_OUT","sum_{u!=v} 1/d(v,u); unreachable contributes 0","outward unnormalized harmonic centrality"),
        ("PRC","PR_v=(1-.85)/N+.85 sum_u PR_u/outdegree_u; dangling mass uniform","damping .85, unit-sum scores; hydrological connectivity, not pollutant mass"),
        ("AC_alpha09","(I-.9*A.T)c=1","unit exogenous vector; require .9*spectral_radius<1"),
        ("AUT","leading right singular vector of A, max-scaled","HITS authority: connections from hubs; degeneracy possible"),
        ("HBS","leading left singular vector of A, max-scaled","HITS hub: connections towards authorities; degeneracy possible"),
        ("UPS","|ancestors(v)|","number of water bodies draining (directly or indirectly) into v; upstream contributing network size"),
        ("DOWN","min_{sink t reachable from v} d(v,t)","recipient steps from v to the network outlet")]:
        add(name,"hydrological directed VOR graph",formula,"dimensionless",note)
    for name,formula,unit in [
        ("N_good","count observed annual values with s*(x-threshold)<=0","annual observations"),
        ("N_bad","count observed annual values with s*(x-threshold)>0","annual observations"),
        ("Amax_bad","max(0,max s*(x-threshold))","parameter unit"),
        ("Amax_good","max(0,max -s*(x-threshold))","parameter unit"),
        ("Lmax_bad_years","longest run of observed adverse consecutive calendar years","observed consecutive annual summaries"),
        ("Lmax_good_years","longest run of observed compliant consecutive calendar years","observed consecutive annual summaries"),
        ("ARL_bad_years","mean length of adverse runs only; absent state NA","observed consecutive annual summaries"),
        ("ARL_good_years","mean length of compliant runs only; absent state NA","observed consecutive annual summaries"),
        ("tR_first_observed_years","first observed compliant year after worst adverse annual year minus worst year","calendar-year difference")]:
        add(name,"water-body/parameter annual series",formula,unit,"Not duration of uninterrupted conditions between actual samples; gaps and censoring explicit")
    for name,formula,note in [
        ("DZI","sum_i(max indegree-indegree_i)/(n-1)^2","directed Freeman degree centralization"),
        ("DZO","sum_i(max outdegree-outdegree_i)/(n-1)^2","directed Freeman degree centralization"),
        ("DZ_U","sum_i(max degree-degree_i)/((n-1)(n-2))","undirected projection"),
        ("BZ_D","sum_i(max BC_raw-BC_raw_i)/((n-1)^2*(n-2))","directed betweenness centralization"),
        ("PRZ_D","sum_i(max PR-PR_i)","unnormalised PageRank concentration, not average PageRank"),
        ("Dens_D","m/(n(n-1))","directed graph; not undirected 2m formula"),
        ("Arcs","m","series-length dependent; secondary, not selected primary"),
        ("Diam_D","max finite positive directed geodesic distance","unreachable paths excluded"),
        ("AVPL_D","mean finite positive directed geodesic distances","reachable pairs only"),
        ("GLE_D","sum_{i!=j,reachable} 1/d(i,j)/(n(n-1))","unreachable pairs contribute 0"),
        ("ALE_U","mean undirected efficiency of induced neighbor subgraphs","explicit undirected projection; low-degree vertices contribute 0"),
        ("LeM_U","undirected modularity of Leiden partition, objective modularity, n_iterations=-1","partition and modularity both use the undirected projection"),
        ("IM_D","directed modularity of Infomap membership; trials=10","partition modularity, not Infomap codelength; fixed seed"),
        ("Mot_D","igraph motifs_randesu_no(size=3)","connected induced three-node motif count; series-length dependent"),
        ("Tra_U","3*triangles/connected_triplets","igraph transitivity ignores direction, made explicit; zero if no triplets"),
        ("VAs_D","sum_i|outdegree_i-indegree_i|/(2*m)","adjacency row/column imbalance"),
        ("RRes_U","mean trapezoid AUC(LCC/n against removed fraction), random permutation","area under the largest-component retention curve under random vertex deletion"),
        ("SRes_U","trapezoid AUC(LCC/n against removed fraction), adaptive max-degree attack","adaptive highest-degree deletion with deterministic ties; not ecological recovery")]:
        add(name,"annual temporal HVG (D=directed,U=undirected projection)",formula,"dimensionless or graph edge steps",note)
    add("conditional_compliance","water body with monitored direct upstream neighbors",
        "sum_y I(good_vy) I(any_observed_upstream_adverse_y) / sum_y I(any_observed_upstream_adverse_y)",
        "fraction","Annual empirical upstream-conditioned compliance, not a causal intrinsic capacity; complete-coverage sensitivity and block-bootstrap intervals")
    return pd.DataFrame(rows)


def self_tests():
    results=[]
    def check(name,condition):
        if not condition: raise AssertionError(name)
        results.append({"test":name,"passed":True})
    s=local_summary([2000,2001,2002,2003,2004,2005],[0,0,2,2,2,0],1,1)
    check("Separate ARL good/bad runs",s['ARL_good_years']==1.5 and s['ARL_bad_years']==3)
    check("Recovery occurs on below-threshold value, not just equality",s['tR_first_observed_years']==3)
    s=local_summary([2000,2001,2003],[0,0,0],1,1)
    check("Missing calendar year splits longest run",s['Lmax_good_years']==2)
    check("Absent adverse run and no recovery event stay NA",pd.isna(s['ARL_bad_years']) and pd.isna(s['tR_first_observed_years']))
    s=local_summary([2000,2001],[8,6],7,-1)
    check("Higher-is-better oxygen classification",s['N_good']==1 and s['N_bad']==1 and s['Amax_bad']==1)
    check("Unrecovered episode is censored",s['recovery_status']=='right_censored_no_observed_recovery')
    s=local_summary([2000,2001],[1,1],1,1)
    check("Threshold equality is compliant",s['N_good']==2 and s['N_bad']==0)
    check("Strict HVG blocks equal-height endpoints",hvg_edges([1,1,1])==[(0,1),(1,2)])
    check("HVG uses both endpoint heights",(0,2) not in hvg_edges([5,4,3]))
    check("HVG bridges lower intermediate node",(0,2) in hvg_edges([5,2,3]))
    rng=np.random.default_rng(42)
    series=[rng.normal(size=int(rng.integers(3,30))) for _ in range(150)]
    series += [rng.integers(-3,4,size=int(rng.integers(3,30))) for _ in range(150)]
    check("Reference-loop and running-maximum HVG agree on 300 series",all(hvg_edges(x,True)==hvg_edges(x) for x in series))
    check("HVG invariant to positive scale and translation",all(hvg_edges(x)==hvg_edges(3*x+10) for x in series))
    try:
        from ts2vg import HorizontalVG
        check("Independent ts2vg edge parity on continuous and tied series",all(sorted(HorizontalVG(directed='left_to_right').build(x).edges)==hvg_edges(x) for x in series))
    except ImportError:
        results.append({"test":"Independent ts2vg edge parity","passed":False,"reason":"optional ts2vg package unavailable"})
    p=hvg_properties([1,1,1],full=False)
    check("Directed density uses m/(n(n-1))",abs(p['Dens_D']-1/3)<1e-12)
    check("Forward efficiency counts unreachable pairs as zero",abs(p['GLE_D']-(1+1+.5)/6)<1e-12)
    check("BH multiplicity adjustment",np.allclose(stats.false_discovery_control([.01,.04,.2]),[.03,.06,.2]))
    graph=nx.DiGraph([('u','v')]);graph.add_node('w')
    c,_=alpha_centrality(graph)
    check("Alpha centrality orientation and unit exogenous vector",abs(c['u']-1)<1e-12 and abs(c['v']-1.9)<1e-12)
    mh=mantel_haenszel([(10,20,5,40)])
    check("Mantel-Haenszel OR reproduces single-stratum OR",abs(mh["OR_MH"]-4.0)<1e-12)
    mh=mantel_haenszel([(10,20,5,40),(10,20,5,40)])
    check("Mantel-Haenszel OR stable for identical strata",abs(mh["OR_MH"]-4.0)<1e-12 and mh["p_CMH"]<.01)
    rr=np.random.default_rng(1)
    xs=np.vstack([rr.normal(0,1,(30,3)),rr.normal(2,1,(30,3))]);lab=np.repeat(["a","b"],30)
    pm=permanova(xs,lab,199,1)
    check("PERMANOVA detects shifted groups",pm["p_perm"]<.01 and pm["R2"]>.3)
    pm=permanova(rr.normal(0,1,(60,3)),lab,199,1)
    check("PERMANOVA null not rejected",pm["p_perm"]>.01)
    g=nx.DiGraph([("a","b"),("b","c"),("d","c")])
    check("Upstream count and outlet distance",len(nx.ancestors(g,"c"))==3 and steps_to_outlet(g,"a")==2)
    tiny=pd.DataFrame({"VOR":["x"]*3,"parameter":["DO"]*3,"year":[2000,2001,2002],"value":[6.,8.,8.],
                       "threshold":[7.]*3,"direction":[-1.]*3,"eligible_local":[True]*3,"Type":["River"]*3,
                       "unit":["g/L"]*3,"n_sites":[1]*3,"basin":["b"]*3})
    strict=calculate_local(tiny,threshold_shift=-.2).iloc[0];lenient=calculate_local(tiny,threshold_shift=.2).iloc[0]
    check("Direction-aware threshold shift (higher-is-better)",abs(strict.threshold-8.4)<1e-12 and abs(lenient.threshold-5.6)<1e-12 and strict.N_bad==3 and lenient.N_bad==0)
    try:
        hvg_edges([1,np.nan,2]);bad=False
    except ValueError:bad=True
    check("Nonfinite HVG input rejected",bad)
    return pd.DataFrame(results)


def main(input_dir=None,output_dir=None):
    started=time.time()
    root=discover_input(input_dir or INPUT_DIR)
    out=Path(output_dir or OUTPUT_DIR or root/"outputs").expanduser().resolve()
    if out == root or root in out.parents or out in root.parents:
        raise ValueError("Input and output directories must not overlap.")
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {out}. Choose a new --output-dir.")
    out.mkdir(parents=True,exist_ok=True)
    tables=out/"tables";tables.mkdir(exist_ok=True)
    figs=out/"figures";figs.mkdir(exist_ok=True)
    logging.basicConfig(level=logging.INFO,format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(),logging.FileHandler(out/"analysis.log",mode="w",encoding="utf-8")],force=True)
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    LOG.info("Surface-water resilience analysis; relative output paths use the selected output directory.")
    random.seed(CFG.seed);np.random.seed(CFG.seed);ig.set_random_number_generator(random.Random(CFG.seed))
    tests=self_tests();write_csv(tests,out/"self_tests.csv")
    LOG.info("Core scientific self-checks passed: %s",int(tests.passed.sum()))
    paths=[input_path(root,f) for f in ["river_registry.csv","lake_registry.csv","water_quality_annual_site_summaries.xlsx","operational_thresholds.xlsx"]]
    rivers=pd.read_csv(paths[0]);lakes=pd.read_csv(paths[1]);raw=excel_read(paths[2]);limits=excel_read(paths[3])
    required={"VOR","KTJ","ANYAG_NEV_ANGOL","ATLAG","DATUM_YEAR","MERTEKEGYSEG_NEV","MERESSZAM",
              "FOLDRAJZI_SZELESSEG","FOLDRAJZI_HOSSZUSAG","VIZTEST_MEGNEVEZES"}
    if not required.issubset(raw):raise ValueError(f"Workbook lacks {sorted(required-set(raw))}")
    if not {"indicator","Lake","River","Boundary"}.issubset(limits):raise ValueError("Threshold workbook schema invalid.")
    if not limits.Boundary.isin(["lower","higher"]).all():raise ValueError("Unknown threshold direction in Boundary.")
    meta,graph,rejected_edges=metadata_and_graph(rivers,lakes,raw)
    write_csv(meta,tables/"hydrological_nodes.csv")
    write_csv(meta[meta.Type.ne(meta.catalogue_priority_type)],tables/"catalogue_type_reconciliation.csv")
    write_csv(pd.DataFrame(sorted(graph.edges),columns=["from","to"]),tables/"recipient_links.csv")
    write_csv(rejected_edges,tables/"excluded_recipient_links.csv")
    nx.write_graphml(graph,out/"hydrological_network.graphml")
    write_csv(limits,tables/"operational_thresholds.csv")
    LOG.info("Network: %s VOR nodes, %s edges; DAG=%s",len(graph),graph.number_of_edges(),nx.is_directed_acyclic_graph(graph))
    annual,site,exclusions,audit=normalize_and_aggregate(raw,meta,limits)
    write_csv(exclusions,tables/"excluded_monitoring_records.csv")
    write_csv(site,tables/"monitoring_site_annual_summaries.csv")
    write_csv(site[site.conflicting_values>1],tables/"conflicting_site_annual_records.csv")
    write_csv(annual,tables/"water_body_annual_summaries.csv")
    coverage,parameter_counts,raw_counts=coverage_tables(annual,site,raw,meta)
    write_csv(coverage,tables/"water_body_parameter_coverage.csv")
    write_csv(parameter_counts,tables/"parameter_effective_sample_sizes.csv")
    write_csv(raw_counts,tables/"source_parameter_unit_counts.csv")
    write_csv(raw.groupby("DATUM_YEAR").size().reset_index(name="raw_rows"),tables/"source_records_by_year.csv")
    local=calculate_local(annual)
    write_csv(local,tables/"local_resilience_indicators.csv")
    write_csv(local[local.VOR.eq(CFG.selected_vor)],tables/"lake_balaton_local_indicators.csv")
    LOG.info("Local resilience indicators: %s body-parameter combinations.",len(local))
    hydro,radius=hydro_properties(graph,meta)
    write_csv(hydro,tables/"hydrological_descriptors.csv")
    # Selection uses observed registered hydrological nodes, before outcomes.
    selected_hydro,hydro_selection=select_nonredundant(hydro[hydro.observed & hydro.registered],HYDRO_PRIORITY,CFG.redundancy_rho)
    write_csv(hydro_selection,tables/"hydrological_metric_selection.csv")
    write_csv(hydro[HYDRO_PRIORITY].corr(method="spearman").reset_index(names="metric"),tables/"hydrological_metric_redundancy.csv")
    LOG.info("HVG extraction and reference-edge parity checks.")
    hvg,hvg_rejected,hvg_differences=calculate_hvgs(annual)
    write_csv(hvg,tables/"hvg_annual_metrics.csv")
    write_csv(hvg_rejected,tables/"hvg_excluded_series.csv")
    write_csv(hvg_differences.reindex(columns=["VOR","parameter"]),tables/"hvg_reference_edge_disagreements.csv")
    if len(hvg_differences):raise AssertionError("Reference-loop HVG parity failed on actual annual series.")
    selected_hvg,hvg_selection=select_nonredundant(hvg,HVG_PRIORITY,CFG.redundancy_rho)
    write_csv(hvg_selection,tables/"hvg_metric_selection.csv")
    write_csv(hvg[HVG_FULL].corr(method="spearman").reset_index(names="metric"),tables/"hvg_metric_redundancy.csv")
    LOG.info("Association analysis; selected hydro=%s, HVG=%s",selected_hydro,selected_hvg)
    local_outcomes=["N_good","N_bad","Amax_bad_relative","Amax_good_relative","Lmax_bad_years","Lmax_good_years",
                    "ARL_bad_years","ARL_good_years","tR_first_observed_years","noncompliance_fraction"]
    eligible_local=local[local.n_years>=CFG.min_years_local]
    corr_local_raw=correlation_table(eligible_local,hydro,selected_hydro,local_outcomes,"local_unadjusted",False)
    corr_local_adj=correlation_table(eligible_local,hydro,selected_hydro,local_outcomes,"local_controlled",True)
    corr_hvg_raw=correlation_table(hvg,hydro,selected_hydro,selected_hvg,"HVG_unadjusted",False)
    corr_hvg_adj=correlation_table(hvg,hydro,selected_hydro,selected_hvg,"HVG_controlled",True)
    for name,table in [("local_associations_unadjusted",corr_local_raw),("local_associations_controlled",corr_local_adj),
                       ("hvg_associations_unadjusted",corr_hvg_raw),("hvg_associations_controlled",corr_hvg_adj)]:
        write_csv(table,tables/(name+".csv"))
    ranks=ranking(local,meta);write_csv(ranks,tables/"water_body_ranking.csv")
    upstream,upstream_events=upstream_conditioned_assessment(annual,graph)
    write_csv(upstream,tables/"upstream_conditioned_compliance.csv")
    write_csv(upstream_events,tables/"upstream_matched_annual_states.csv")

    placebo_events=placebo_upstream_events(annual,graph)
    write_csv(placebo_events,tables/"placebo_annual_states.csv")
    mh_table,upstream_ext=topological_resilience_tests(upstream,upstream_events,placebo_events)
    write_csv(mh_table,tables/"topological_resilience_tests.csv")
    write_csv(upstream_ext,tables/"topological_resilience_by_parameter.csv")
    if len(upstream_ext):
        body_tr=upstream_ext[upstream_ext.adequate_exposed_years].groupby("VOR",as_index=False).agg(
            TR_mean=("conditional_compliance","mean"),buffering_contrast_mean=("buffering_contrast","mean"),
            n_parameters=("parameter","nunique"))
        write_csv(body_tr,tables/"topological_resilience_by_water_body.csv")

    typology_table,typology_params,typology_imputed=typology_tests(annual,meta)
    write_csv(typology_table,tables/"typology_permanova.csv")
    write_csv(typology_params,tables/"typology_parameter_tests.csv")
    pca_var,pca_load=hvg_pca(hvg)
    write_csv(pca_var,tables/"hvg_pca_variance.csv");write_csv(pca_load,tables/"hvg_pca_loadings.csv")
    sensitivity_summary=pd.DataFrame();ranking_envelope=pd.DataFrame();common=[]
    rank_agreement=pd.DataFrame();agg_summary=pd.DataFrame();selected_series=pd.DataFrame()
    if CFG.run_sensitivity:
        LOG.info("Ranking and HVG sensitivity analysis.")
        rank_details,ranking_envelope,common,rank_agreement=ranking_sensitivity(annual,local,ranks,meta)
        write_csv(rank_agreement,tables/"ranking_scenario_agreement.csv")
        write_csv(rank_details,tables/"ranking_sensitivity_scenarios.csv")
        write_csv(ranking_envelope,tables/"ranking_year_omission_envelopes.csv")
        sensitivity_detail,sensitivity_summary,agg_sensitivity,selected_series=hvg_sensitivity(annual,hvg)
        write_csv(sensitivity_detail,tables/"hvg_deletion_replicates.csv")
        write_csv(sensitivity_summary,tables/"hvg_deletion_concordance.csv")
        write_csv(agg_sensitivity,tables/"hvg_aggregation_sensitivity.csv")
        agg_summary=aggregation_summary(agg_sensitivity)
        write_csv(agg_summary,tables/"hvg_aggregation_concordance.csv")
        write_csv(selected_series,tables/"hvg_sensitivity_sample.csv")
        write_csv(series_length_sensitivity(hvg,hydro,selected_hydro,annual),tables/"hvg_series_length_sensitivity.csv")
    rbmp,rbmp_status=optional_rbmp(root,ranks)
    if not rbmp.empty:
        write_csv(rbmp,tables/"official_status_comparison.csv")
    write_csv(metric_dictionary(),tables/"metric_dictionary.csv")
    LOG.info("Exporting analysis figures.")
    # Per-parameter heatmaps (Supplementary) use the CONTROLLED tables for both families.
    figures(annual,local,hydro,ranks,corr_local_adj,corr_hvg_adj,figs)
    association_forest(corr_local_adj,figs/"local_association_forest",
                       f"Strongest controlled associations: hydrological position vs local indicators (top {CFG.top_associations})")
    association_forest(corr_hvg_adj,figs/"hvg_association_forest",
                       f"Strongest controlled associations: hydrological position vs temporal HVG metrics (top {CFG.top_associations})")
    topological_figure(upstream_ext,mh_table,figs/"topological_resilience")
    hvg_sensitivity_figure(sensitivity_summary,figs/"hvg_deletion_sensitivity")
    # Figure 4: actual annual monitoring-site summaries, not invented raw samples.
    ex=site[site.VOR.eq(example_vor(annual)) & site.parameter.eq("BOD (5)")]
    if len(ex):
        fig,ax=plt.subplots(figsize=(9,5.2))
        for station,g in ex.groupby("site"):
            ax.plot(g.year,g.value,"o-",lw=.8,alpha=.45,label=station)
        aggregated=annual[annual.VOR.eq(example_vor(annual)) & annual.parameter.eq("BOD (5)")]
        ax.plot(aggregated.year,aggregated.value,"o-",color="black",lw=2.5,label="Equal-site annual body mean")
        ax.set(xlabel="Calendar year",ylabel="BOD$_5$ [g/L]",
               title=f"{example_vor(annual)}: monitoring-site annual summaries\nIndividual sample dates unavailable")
        ax.legend(title="KTJ monitoring site",fontsize=9,ncol=3,frameon=False);ax.grid(alpha=.15)
        save_figure(fig,figs/"site_annual_comparison")
    tsne_meta={"executed":False}
    if CFG.run_tsne:
        LOG.info("Exploratory physicochemical-profile embedding.")
        tsne_meta=quality_profile_tsne(annual,meta,figs)
    versions={"Python":platform.python_version()}
    for package in ["numpy","pandas","scipy","networkx","igraph","scikit-learn","matplotlib","seaborn","openpyxl","python-calamine","ts2vg"]:
        try:versions[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:versions[package]="unavailable"
    source_files=paths
    manifest={"configuration":asdict(CFG),"versions":versions,"elapsed_seconds":round(time.time()-started,2),
              "input_files":[{"filename":p.name,"size":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in source_files],
              "data_audit":audit,
              "raw_VOR_count":int(raw.VOR.nunique()),"raw_monitoring_KTJ_count":int(raw.KTJ.nunique()),
              "raw_missing_VOR_rows":int(raw.VOR.isna().sum()),
              "raw_parameter_count":int(raw.ANYAG_NEV_ANGOL.nunique()),
              "river_catalogue_rows":len(rivers),"lake_catalogue_rows":len(lakes),
              "registered_rivers":int((meta.registered & meta.Type.eq("River")).sum()),
              "registered_lakes":int((meta.registered & meta.Type.eq("Lake")).sum()),
              "observed_types":meta[meta.observed].Type.value_counts().to_dict(),
              "overlapping_river_lake_catalogue_VORs":len(set(rivers.VT_VOR_2015)&set(lakes.VT_VOR_2015)),
              "observed_reclassified_LW_VORs":int((meta.observed & meta.Type.ne(meta.catalogue_priority_type)).sum()),
              "network_nodes":len(graph),"network_edges":graph.number_of_edges(),
              "network_spectral_radius":radius,"network_DAG":nx.is_directed_acyclic_graph(graph),
              "local_indicator_rows":len(local),"HVG_rows":len(hvg),
              "selected_hydrological_metrics":selected_hydro,"selected_HVG_metrics":selected_hvg,
              "common_ranking_parameters":common,"tsne":tsne_meta,"RBMP_status":rbmp_status,
              "reference_edge_disagreements":len(hvg_differences),
              "analysis_units":"one registered VOR/parameter with annual summaries; not individual samples",
              "inference_scope":"annual physicochemical associations and operational threshold compliance; no causal or full ecological claim",
              "program_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out/"manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2,default=str),encoding="utf-8")

    export_info=export_generated(out/"latex",figs,dict(
        manifest=manifest,meta=meta,local=local,hvg=hvg,hvg_rejected=hvg_rejected,ranks=ranks,annual=annual,
        coverage=coverage,parameter_counts=parameter_counts,selected_hydro=selected_hydro,selected_hvg=selected_hvg,
        hydro_selection=hydro_selection,hvg_selection=hvg_selection,
        corr_local_adj=corr_local_adj,corr_hvg_adj=corr_hvg_adj,corr_local_raw=corr_local_raw,corr_hvg_raw=corr_hvg_raw,
        rank_agreement=rank_agreement,rank_envelope=ranking_envelope,typology=typology_table,
        typology_params=typology_params,typology_imputed=typology_imputed,mh=mh_table,upstream_ext=upstream_ext,
        sens_summary=sensitivity_summary,sens_sample=selected_series,agg_summary=agg_summary,pca_var=pca_var,
        rbmp=rbmp,quick=QUICK_RUN))
    LOG.info("LaTeX export: %s macros; figures %s -> %s",export_info["macros"],export_info["figures"],out/"latex")
    LOG.info("Analysis completed in %.1f minutes; outputs: tables, figures, latex, manifest.",(time.time()-started)/60)
    LOG.info("Official RBMP comparison: %s",rbmp_status)
    return manifest


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir");parser.add_argument("--output-dir")
    parser.add_argument("--self-test-only",action="store_true")
    parser.add_argument("--quick",action="store_true",help="20 bootstrap and 3 deletion/attack reps; max 10 series per parameter")
    parser.add_argument("--skip-tsne",action="store_true");parser.add_argument("--skip-sensitivity",action="store_true")
    parser.add_argument("--min-years",type=int,default=None)
    args=parser.parse_args()
    if args.quick:
        QUICK_RUN=True;CFG.permutations=99
        CFG.bootstrap_reps=20;CFG.sensitivity_reps=3;CFG.attack_reps=3;CFG.sensitivity_series_per_parameter=10
    if args.skip_tsne:CFG.run_tsne=False
    if args.skip_sensitivity:CFG.run_sensitivity=False
    if args.min_years is not None:
        if args.min_years<3:parser.error("--min-years must be >=3")
        CFG.min_years_local=CFG.min_years_hvg=args.min_years
    if args.self_test_only:
        print(self_tests().to_string(index=False))
    else:
        main(args.input_dir,args.output_dir)

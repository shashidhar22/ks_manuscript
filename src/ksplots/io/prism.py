"""GraphPad Prism file parser.

Prism 9+ files are ZIP archives containing JSON metadata and CSV data tables.
This module extracts data tables into pandas DataFrames.
"""
from __future__ import annotations

import csv
import io
import json
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def list_tables(prism_path: str | Path) -> list[dict[str, Any]]:
    """List all data tables in a Prism file.

    Returns a list of dicts with keys: title, format, dataFormat,
    replicatesCount, columns.
    """
    prism_path = Path(prism_path)
    results = []
    with zipfile.ZipFile(prism_path) as z:
        sheets = [
            f for f in z.namelist()
            if f.startswith("data/sheets/") and f.endswith("sheet.json")
        ]
        for s in sheets:
            sheet = json.loads(z.read(s))
            title = sheet.get("title")
            if not title:
                continue
            table = sheet.get("table", {})
            table_uid = table.get("uid", "")
            csv_path = f"data/tables/{table_uid}/data.csv"
            if csv_path not in z.namelist():
                continue

            ds_titles = _get_dataset_titles(z, table.get("dataSets", []))
            results.append({
                "title": title,
                "format": table.get("format", ""),
                "dataFormat": table.get("dataFormat", ""),
                "replicatesCount": table.get("replicatesCount", 1),
                "columns": ds_titles,
                "_uid": table_uid,
            })
    return results


def extract_table(
    prism_path: str | Path,
    title: str,
) -> pd.DataFrame:
    """Extract a single data table by title.

    For **XY** tables (dose-response curves):
        Returns DataFrame with columns: x, {col_name}_r1, {col_name}_r2, ...
        where col_name comes from the dataset title.

    For **grouped** tables (bar charts):
        Returns DataFrame with columns: label, {col_name}_r1, {col_name}_r2, ...

    For **column** tables:
        Returns DataFrame with columns: label, {col_name}_r1, ...
    """
    prism_path = Path(prism_path)
    with zipfile.ZipFile(prism_path) as z:
        sheet = _find_sheet(z, title)
        if sheet is None:
            available = [t["title"] for t in list_tables(prism_path)]
            raise KeyError(
                f"Table {title!r} not found. "
                f"Available ({len(available)}): {available[:10]}..."
            )

        table = sheet["table"]
        table_uid = table["uid"]
        fmt = table.get("format", "")
        n_reps = table.get("replicatesCount", 1)
        ds_titles = _get_dataset_titles(z, table.get("dataSets", []))

        csv_data = z.read(f"data/tables/{table_uid}/data.csv").decode()
        rows = list(csv.reader(io.StringIO(csv_data)))

        if fmt == "xy":
            return _parse_xy(rows, ds_titles, n_reps)
        elif fmt in ("grouped", "column"):
            return _parse_grouped(rows, ds_titles, n_reps)
        else:
            return _parse_generic(rows, ds_titles, n_reps)


def extract_xy_mean_sem(
    prism_path: str | Path,
    title: str,
) -> pd.DataFrame:
    """Extract XY table and compute mean +/- SEM per dataset.

    Returns DataFrame with columns: x, {name}_mean, {name}_sem for each
    dataset column. Convenience wrapper for dose-response plotting.
    """
    df = extract_table(prism_path, title)
    if "x" not in df.columns:
        raise ValueError(f"Table {title!r} is not XY format (no 'x' column)")

    result = pd.DataFrame({"x": df["x"]})
    # Group replicate columns by base name
    base_names = []
    for col in df.columns:
        if col == "x":
            continue
        base = col.rsplit("_r", 1)[0]
        if base not in base_names:
            base_names.append(base)

    for base in base_names:
        rep_cols = [c for c in df.columns if c.startswith(base + "_r")]
        vals = df[rep_cols].astype(float)
        result[f"{base}_mean"] = vals.mean(axis=1)
        result[f"{base}_sem"] = vals.sem(axis=1)

    return result


def extract_grouped_mean_sem(
    prism_path: str | Path,
    title: str,
) -> pd.DataFrame:
    """Extract grouped table and compute mean +/- SEM per condition.

    Returns DataFrame with columns: label, {condition}_mean, {condition}_sem.
    """
    df = extract_table(prism_path, title)
    if "label" not in df.columns:
        raise ValueError(f"Table {title!r} has no 'label' column")

    # Drop empty rows
    data_cols = [c for c in df.columns if c != "label"]
    mask = df[data_cols].notna().any(axis=1)
    df = df[mask].copy()

    result = pd.DataFrame({"label": df["label"]})
    base_names = []
    for col in data_cols:
        base = col.rsplit("_r", 1)[0]
        if base not in base_names:
            base_names.append(base)

    for base in base_names:
        rep_cols = [c for c in df.columns if c.startswith(base + "_r")]
        vals = df[rep_cols].astype(float)
        result[f"{base}_mean"] = vals.mean(axis=1).values
        result[f"{base}_sem"] = vals.sem(axis=1).values

    return result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_dataset_titles(z: zipfile.ZipFile, ds_ids: list[str]) -> list[str]:
    titles = []
    seen: dict[str, int] = {}
    for ds_id in ds_ids:
        ds_path = f"data/sets/{ds_id}.json"
        if ds_path in z.namelist():
            ds = json.loads(z.read(ds_path))
            t = ds.get("title", {})
            if isinstance(t, dict):
                # Take first line of multi-line title
                raw = t.get("string", "")
                name = raw.split("\r")[0].strip()
            elif isinstance(t, str):
                name = t.strip()
            else:
                name = f"Col_{len(titles) + 1}"
            name = name if name else f"Col_{len(titles) + 1}"
        else:
            name = f"Col_{len(titles) + 1}"

        # Disambiguate duplicate names
        if name in seen:
            seen[name] += 1
            name = f"{name} ({seen[name]})"
        else:
            seen[name] = 1
        titles.append(name)
    return titles


def _find_sheet(z: zipfile.ZipFile, title: str) -> dict | None:
    sheets = [
        f for f in z.namelist()
        if f.startswith("data/sheets/") and f.endswith("sheet.json")
    ]
    for s in sheets:
        sheet = json.loads(z.read(s))
        if sheet.get("title") == title:
            return sheet
    return None


def _safe_float(val: str) -> float | None:
    if not val or val.strip() == "":
        return None
    try:
        return float(val)
    except ValueError:
        return None


def _parse_xy(
    rows: list[list[str]],
    ds_titles: list[str],
    n_reps: int,
) -> pd.DataFrame:
    """Parse XY format: col0=row_title, col1=x, then n_reps per dataset."""
    records = []
    for row in rows:
        if not row or all(v.strip() == "" for v in row):
            continue
        # First two columns are row title and x value
        x_val = _safe_float(row[1]) if len(row) > 1 else None
        if x_val is None:
            continue

        rec: dict[str, Any] = {"x": x_val}
        col_idx = 2
        for ds_name in ds_titles:
            for r in range(n_reps):
                col_name = f"{ds_name}_r{r + 1}"
                if col_idx < len(row):
                    rec[col_name] = _safe_float(row[col_idx])
                else:
                    rec[col_name] = None
                col_idx += 1
        records.append(rec)

    return pd.DataFrame(records)


def _parse_grouped(
    rows: list[list[str]],
    ds_titles: list[str],
    n_reps: int,
) -> pd.DataFrame:
    """Parse grouped/column format: col0=label, then n_reps per dataset."""
    records = []
    for row in rows:
        if not row or all(v.strip() == "" for v in row):
            continue
        label = row[0].strip().strip('"').replace("\r", " ") if row else ""
        if not label:
            continue

        rec: dict[str, Any] = {"label": label}
        col_idx = 1
        for ds_name in ds_titles:
            for r in range(n_reps):
                col_name = f"{ds_name}_r{r + 1}"
                if col_idx < len(row):
                    rec[col_name] = _safe_float(row[col_idx])
                else:
                    rec[col_name] = None
                col_idx += 1
        records.append(rec)

    return pd.DataFrame(records)


def _parse_generic(
    rows: list[list[str]],
    ds_titles: list[str],
    n_reps: int,
) -> pd.DataFrame:
    """Fallback parser — treat as grouped."""
    return _parse_grouped(rows, ds_titles, n_reps)

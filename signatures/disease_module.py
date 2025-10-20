"""Disease module assembly and meta-analysis utilities."""
from __future__ import annotations

import json
import logging
import math
import pathlib
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


@dataclass
class ExpressionStudy:
    """Container describing an expression study used for disease signatures."""

    study_id: str
    disease_id: str
    expression: pd.DataFrame
    case_samples: Sequence[str]
    control_samples: Sequence[str]
    gene_column: str = "gene_symbol"

    def run(self) -> pd.DataFrame:
        """Run differential expression for the study.

        Returns a dataframe with columns:
        `disease_id`, `study_id`, `gene_symbol`, `logfc`, `se`, `pvalue`.
        """

        matrix = self.expression.copy()
        if self.gene_column not in matrix.columns:
            raise KeyError(f"Expression table missing gene column '{self.gene_column}'")

        missing = [
            sample
            for sample in list(self.case_samples) + list(self.control_samples)
            if sample not in matrix.columns
        ]
        if missing:
            raise KeyError(f"Expression table missing sample columns: {missing}")

        matrix = matrix.set_index(self.gene_column).astype(float)
        case = matrix[list(self.case_samples)]
        control = matrix[list(self.control_samples)]

        case_mean = case.mean(axis=1)
        control_mean = control.mean(axis=1)
        logfc = np.log2((case_mean + 1e-6) / (control_mean + 1e-6))

        case_var = case.var(axis=1, ddof=1)
        control_var = control.var(axis=1, ddof=1)
        case_n = len(self.case_samples)
        control_n = len(self.control_samples)

        se = np.sqrt(case_var / case_n + control_var / control_n)
        se.replace(0, np.nan, inplace=True)

        with np.errstate(divide="ignore", invalid="ignore"):
            t_stat = logfc / se

        dof_num = (case_var / case_n + control_var / control_n) ** 2
        dof_den = (
            (case_var**2) / ((case_n**2) * max(case_n - 1, 1))
            + (control_var**2) / ((control_n**2) * max(control_n - 1, 1))
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            dof = dof_num / dof_den

        try:
            from scipy.stats import t as student_t  # type: ignore

            pvalues = 2 * student_t.sf(np.abs(t_stat), dof)
        except Exception:  # pragma: no cover - SciPy optional
            pvalues = 2 * (1 - 0.5 * (1 + np.vectorize(_erf)(np.abs(t_stat) / math.sqrt(2))))

        result = pd.DataFrame(
            {
                "disease_id": self.disease_id,
                "study_id": self.study_id,
                "gene_symbol": matrix.index,
                "logfc": logfc.values,
                "se": se.values,
                "pvalue": np.asarray(pvalues),
            }
        )
        return result.dropna(subset=["logfc", "se", "pvalue"])


def _erf(x: float) -> float:
    """Numerically stable error function fallback."""

    # Polynomial approximation (Abramowitz & Stegun 7.1.26)
    sign = 1 if x >= 0 else -1
    x = abs(x)
    t = 1.0 / (1.0 + 0.5 * x)
    tau = t * math.exp(
        -x * x
        - 1.26551223
        + 1.00002368 * t
        + 0.37409196 * t**2
        + 0.09678418 * t**3
        - 0.18628806 * t**4
        + 0.27886807 * t**5
        - 1.13520398 * t**6
        + 1.48851587 * t**7
        - 0.82215223 * t**8
        + 0.17087277 * t**9
    )
    return sign * (1 - tau)


def run_study_analysis(studies: Iterable[ExpressionStudy]) -> List[pd.DataFrame]:
    """Execute differential expression across all studies."""

    outputs: List[pd.DataFrame] = []
    for study in studies:
        try:
            outputs.append(study.run())
        except Exception as exc:  # pragma: no cover - defensive guard
            logger.error("Study %s failed: %s", study.study_id, exc)
    return outputs


def meta_analyse_studies(study_tables: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """Perform DerSimonian–Laird random-effects meta-analysis across studies."""

    if not study_tables:
        return pd.DataFrame(
            columns=[
                "disease_id",
                "gene_symbol",
                "meta_logfc",
                "meta_se",
                "meta_pvalue",
                "tau2",
                "i2",
                "q",
                "n_studies",
            ]
        )

    combined = pd.concat(study_tables, ignore_index=True, sort=False)
    combined = combined.dropna(subset=["disease_id", "gene_symbol", "logfc", "se"])
    if combined.empty:
        return pd.DataFrame(
            columns=[
                "disease_id",
                "gene_symbol",
                "meta_logfc",
                "meta_se",
                "meta_pvalue",
                "tau2",
                "i2",
                "q",
                "n_studies",
            ]
        )

    def _per_group(df: pd.DataFrame) -> pd.Series:
        df = df.copy()
        df["variance"] = df["se"] ** 2
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["variance"])
        df = df[df["variance"] > 0]
        if df.empty:
            return pd.Series(
                {
                    "meta_logfc": np.nan,
                    "meta_se": np.nan,
                    "meta_pvalue": np.nan,
                    "tau2": np.nan,
                    "i2": np.nan,
                    "q": np.nan,
                    "n_studies": 0,
                }
            )

        weights = 1.0 / df["variance"]
        effects = df["logfc"]
        weighted_mean = (weights * effects).sum() / weights.sum()
        q = (weights * (effects - weighted_mean) ** 2).sum()
        df_deg = max(len(df) - 1, 0)

        c = weights.sum() - (weights**2).sum() / weights.sum() if weights.sum() else 0.0
        tau2 = max((q - df_deg) / c, 0.0) if c > 0 else 0.0

        random_weights = 1.0 / (df["variance"] + tau2)
        weight_sum = random_weights.sum()
        meta = (random_weights * effects).sum() / weight_sum
        se = math.sqrt(1.0 / weight_sum)
        z = meta / se if se > 0 else 0.0
        pvalue = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
        i2 = max((q - df_deg) / q, 0.0) if q > 0 else 0.0

        return pd.Series(
            {
                "meta_logfc": meta,
                "meta_se": se,
                "meta_pvalue": pvalue,
                "tau2": tau2,
                "i2": i2,
                "q": q,
                "n_studies": len(df),
            }
        )

    aggregated = (
        combined.groupby(["disease_id", "gene_symbol"], as_index=False)
        .apply(_per_group)
        .reset_index()
    )
    aggregated = aggregated.drop(columns=["level_2"], errors="ignore")
    return aggregated


def integrate_disease_evidence(
    meta_expression: pd.DataFrame,
    seeds: Optional[pd.DataFrame] = None,
    expression_weight: float = 0.6,
    seed_weight: float = 0.4,
) -> pd.DataFrame:
    """Combine meta-analysed expression with genetics/curated seeds into disease_gene scores."""

    if meta_expression.empty:
        return pd.DataFrame(
            columns=[
                "disease_id",
                "gene_symbol",
                "combined_score",
                "expression_score",
                "seed_score",
                "n_studies",
                "sources",
            ]
        )

    combined = meta_expression.copy()
    combined = combined.dropna(subset=["disease_id", "gene_symbol", "meta_logfc", "meta_pvalue"])
    combined["expression_score_raw"] = _expression_score(combined["meta_logfc"], combined["meta_pvalue"])

    if seeds is not None and not seeds.empty:
        seed_table = seeds.copy()
        seed_table = seed_table.dropna(subset=["disease_id", "gene_symbol", "score"])
        seed_summary = (
            seed_table.groupby(["disease_id", "gene_symbol"], as_index=False)
            .agg(seed_score=("score", "mean"))
        )
        combined = combined.merge(seed_summary, on=["disease_id", "gene_symbol"], how="left")
    else:
        combined["seed_score"] = np.nan

    combined["expression_score"] = combined.groupby("disease_id")["expression_score_raw"].transform(_normalise)
    combined["seed_score_norm"] = combined.groupby("disease_id")["seed_score"].transform(_normalise)

    weight_sum = max(expression_weight + seed_weight, 1e-9)
    combined_score = (
        expression_weight * combined["expression_score"].fillna(0)
        + seed_weight * combined["seed_score_norm"].fillna(0)
    ) / weight_sum
    combined["combined_score"] = combined_score.clip(0, 1)

    def build_sources(row: pd.Series) -> List[str]:
        evidence: List[str] = []
        if not math.isnan(row.get("expression_score", np.nan)):
            evidence.append("expression")
        if not math.isnan(row.get("seed_score", np.nan)):
            evidence.append("seed")
        return evidence

    combined["sources"] = combined.apply(build_sources, axis=1)
    return combined[
        [
            "disease_id",
            "gene_symbol",
            "combined_score",
            "expression_score",
            "seed_score",
            "n_studies",
            "sources",
        ]
    ].sort_values(["disease_id", "combined_score"], ascending=[True, False])


def _expression_score(logfc: pd.Series, pvalue: pd.Series) -> pd.Series:
    with np.errstate(divide="ignore"):
        signed_logp = np.sign(logfc) * (-np.log10(pvalue.clip(lower=1e-12)))
    signed_logp = signed_logp.replace([np.inf, -np.inf], np.nan)
    signed_logp.fillna(0, inplace=True)
    return signed_logp


def _normalise(series: pd.Series) -> pd.Series:
    if series.isna().all():
        return series
    valid = series.dropna()
    if valid.empty:
        return series
    min_val = valid.min()
    max_val = valid.max()
    if math.isclose(min_val, max_val):
        return series.fillna(0)
    scaled = (series - min_val) / (max_val - min_val)
    return scaled.fillna(0)


def load_study_config(path: str | pathlib.Path) -> List[ExpressionStudy]:
    """Load expression studies from a JSON or CSV manifest."""

    path = pathlib.Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    if path.suffix.lower() in {".json", ".jsonl"}:
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        entries = payload if isinstance(payload, list) else payload.get("studies", [])
        records = pd.DataFrame(entries)
    else:
        records = pd.read_csv(path)

    required = {"study_id", "disease_id", "expression_path", "case_samples", "control_samples"}
    missing = required.difference(records.columns)
    if missing:
        raise ValueError(f"Study config missing columns: {sorted(missing)}")

    studies: List[ExpressionStudy] = []
    for row in records.to_dict(orient="records"):
        gene_column = row.get("gene_column", "gene_symbol")
        case_samples = _parse_list(row["case_samples"])
        control_samples = _parse_list(row["control_samples"])
        expression = _load_expression_table(row["expression_path"])
        studies.append(
            ExpressionStudy(
                study_id=str(row["study_id"]),
                disease_id=str(row["disease_id"]),
                expression=expression,
                case_samples=case_samples,
                control_samples=control_samples,
                gene_column=str(gene_column),
            )
        )
    return studies


def _parse_list(value: object) -> List[str]:
    if isinstance(value, (list, tuple)):
        return [str(item) for item in value]
    if isinstance(value, str):
        return [item.strip() for item in value.replace(";", ",").split(",") if item.strip()]
    raise TypeError(f"Cannot parse list from value: {value!r}")


def _load_expression_table(path: str | pathlib.Path) -> pd.DataFrame:
    file_path = pathlib.Path(path)
    if not file_path.exists():
        raise FileNotFoundError(file_path)
    suffix = file_path.suffix.lower()
    if suffix in {".parquet", ".pq"}:
        return pd.read_parquet(file_path)
    if suffix in {".tsv", ".tab"}:
        return pd.read_csv(file_path, sep="\t")
    return pd.read_csv(file_path)


def save_dataframe(df: pd.DataFrame, path: str | pathlib.Path) -> None:
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() in {".parquet", ".pq"}:
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, index=False)


def run_pipeline(
    study_manifest: str,
    output_expression: str,
    output_disease_gene: str,
    seed_path: Optional[str] = None,
    expression_weight: float = 0.6,
    seed_weight: float = 0.4,
) -> None:
    """Execute the disease module pipeline end-to-end."""

    studies = load_study_config(study_manifest)
    if not studies:
        raise ValueError("Study manifest produced no studies")

    tables = run_study_analysis(studies)
    meta = meta_analyse_studies(tables)
    save_dataframe(meta, output_expression)

    seeds = pd.read_csv(seed_path) if seed_path else None
    integrated = integrate_disease_evidence(meta, seeds, expression_weight, seed_weight)
    save_dataframe(integrated, output_disease_gene)


def build_arg_parser() -> "argparse.ArgumentParser":  # pragma: no cover - CLI wiring
    import argparse

    parser = argparse.ArgumentParser(description="Build disease gene modules from expression studies")
    parser.add_argument("--study-manifest", required=True, help="JSON/CSV manifest describing expression studies")
    parser.add_argument("--output-expression", required=True, help="Path to write meta-analysed expression table")
    parser.add_argument("--output-disease-gene", required=True, help="Path to write integrated disease-gene scores")
    parser.add_argument("--seed-path", help="Optional CSV containing seed evidence")
    parser.add_argument("--expression-weight", type=float, default=0.6)
    parser.add_argument("--seed-weight", type=float, default=0.4)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:  # pragma: no cover - CLI wiring
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    run_pipeline(
        study_manifest=args.study_manifest,
        output_expression=args.output_expression,
        output_disease_gene=args.output_disease_gene,
        seed_path=args.seed_path,
        expression_weight=args.expression_weight,
        seed_weight=args.seed_weight,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI wiring
    raise SystemExit(main())

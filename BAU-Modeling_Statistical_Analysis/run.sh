#!/usr/bin/env bash
set -euo pipefail

export BAU_CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$BAU_CODE_ROOT/settings.sh"
export PYTHONPATH="$BAU_CODE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1
cd "$BAU_CODE_ROOT"

stage="${1:-help}"
dry_run=false
if [[ "${2:-}" == --dry-run ]]; then
  dry_run=true
elif [[ $# -gt 1 ]]; then
  printf 'Unknown argument: %s\n' "$2" >&2
  exit 2
fi
if [[ $# -gt 2 ]]; then
  printf 'Usage: bash run.sh [stage] [--dry-run]\n' >&2
  exit 2
fi
model_options=()
if [[ "$BAU_USE_FROZEN_PARAMETERS" == 1 ]]; then
  model_options+=(--use-frozen-parameters)
fi

run() {
  if "$dry_run"; then
    printf '%q ' "$@"
    printf '\n'
  else
    "$@"
  fi
}

run_stage() {
  case "$1" in
    prepare)
      run "$BAU_PYTHON" analysis/prepare/ukb.py
      run "$BAU_PYTHON" analysis/prepare/supporting.py
      run "$BAU_PYTHON" analysis/prepare/organ_panels.py
      run Rscript analysis/prepare/external.R
      ;;
    models)
      for modality in nmr protein blood; do
        for fold in 1 2 3 4 5; do
          run "$BAU_PYTHON" analysis/models/ukb.py run-fold \
            --modality "$modality" --model ngboost --fold "$fold" \
            --threads "$BAU_THREADS" "${model_options[@]}"
        done
        run "$BAU_PYTHON" analysis/models/ukb.py aggregate \
          --modality "$modality" --model ngboost
      done
      for model in lightgbmlss xgboostlss; do
        for modality in nmr protein; do
          for fold in 1 2 3 4 5; do
            run "$BAU_PYTHON" analysis/models/ukb.py run-fold \
              --modality "$modality" --model "$model" --fold "$fold" \
              --threads "$BAU_THREADS" "${model_options[@]}"
          done
          run "$BAU_PYTHON" analysis/models/ukb.py aggregate \
            --modality "$modality" --model "$model"
        done
      done
      for cohort in nhanes charls chns; do
        run "$BAU_PYTHON" analysis/models/external.py "$cohort" \
          --threads "$BAU_THREADS" "${model_options[@]}"
      done
      ;;
    phenotypes)
      run Rscript analysis/phenotypes/model_diagnostics.R
      run Rscript analysis/phenotypes/matched_examples.R
      run "$BAU_PYTHON" analysis/phenotypes/split_domain_reference.py
      run Rscript analysis/phenotypes/domain_clocks.R
      run Rscript analysis/phenotypes/molecular_support.R
      run "$BAU_PYTHON" analysis/phenotypes/molecular_disagreement.py
      run Rscript analysis/phenotypes/disagreement_associations.R
      run "$BAU_PYTHON" analysis/phenotypes/score_repeat.py
      run Rscript analysis/phenotypes/repeatability.R
      run "$BAU_PYTHON" analysis/phenotypes/explainability.py
      run Rscript analysis/phenotypes/algorithm_robustness.R
      ;;
    outcomes)
      for model in ngboost lightgbmlss xgboostlss; do
        run Rscript analysis/outcomes/ukb_associations.R --model="$model"
      done
      run Rscript analysis/outcomes/discrimination.R
      run Rscript analysis/outcomes/longitudinal.R
      run Rscript analysis/outcomes/sensitivity.R
      run Rscript analysis/outcomes/subgroups.R
      run Rscript analysis/outcomes/external_mortality.R
      ;;
    gene-sets)
      run Rscript analysis/molecular/prepare_gene_sets.R
      ;;
    molecular)
      run Rscript analysis/molecular/cross_omic_associations.R
      run Rscript analysis/molecular/physiology.R
      run Rscript analysis/molecular/sensitivity.R
      run Rscript analysis/molecular/healthy_associations.R
      run Rscript analysis/molecular/enrichment.R --population=primary
      run Rscript analysis/molecular/enrichment.R --population=healthy
      ;;
    tables)
      run "$BAU_PYTHON" analysis/tables/baseline_characteristics.py
      run "$BAU_PYTHON" analysis/tables/analysis_tables.py
      ;;
    figures)
      run Rscript figures/schematics.R
      run Rscript figures/construct.R
      run Rscript figures/outcomes.R
      run Rscript figures/molecular.R
      ;;
    test)
      run "$BAU_PYTHON" -m unittest discover -s tests -v
      run Rscript tests/test_statistics.R
      ;;
    all)
      for next_stage in prepare models phenotypes outcomes molecular tables figures; do
        run_stage "$next_stage"
      done
      ;;
    help|-h|--help)
      printf '%s\n' 'Stages: prepare, models, phenotypes, outcomes, molecular, tables, figures, gene-sets, test, all'
      ;;
    *)
      printf 'Unknown stage: %s\n' "$1" >&2
      exit 2
      ;;
  esac
}

run_stage "$stage"

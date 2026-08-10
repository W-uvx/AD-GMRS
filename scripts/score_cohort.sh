#!/usr/bin/env bash
set -euo pipefail
source "${GMRS_ROOT}/scripts/lib.sh"

target="$1"; bfile="$(resolve_path "$2")"; score_root="${SCORE_ROOT:-$GMRS_ROOT/results/runs/training/scores}"
[[ -f "${bfile}.bed" && -f "${bfile}.bim" && -f "${bfile}.fam" ]] || die "incomplete PLINK files: $bfile"
mkdir -p "$score_root"

for layer in $(layers_for "$target"); do
  manifest="$(manifest_for "$layer")"; outdir="$score_root/$layer"; mkdir -p "$outdir/logs"
  printf 'feature_id\tweight_snps\tmatched_snps\tmatch_rate\n' > "$outdir/snp_qc.tsv"
  printf 'feature_id\tsha256\n' > "$outdir/weight_checksums.tsv"
  while IFS=$'\t' read -r feature chr _summary _n _; do
    [[ "$feature" == feature_id ]] && continue
    weight_prefix="$(weight_dir "$layer")/$feature"; out="$outdir/$feature"
    mapfile -t weights < <(compgen -G "${weight_prefix}_pst_eff_a1_b0.5_phiauto*.txt" || true)
    [[ ${#weights[@]} -gt 0 ]] || die "no PRS-CS weight for $layer/$feature"
    merged="$outdir/.${feature}.weights.noheader.tsv"
    tmp_weights="${merged}.tmp"
    if ! awk 'NF>=6 {signature=$1 FS $3 FS $4 FS $5 FS $6; if (($2 in signatures) && signatures[$2]!=signature) {print "conflicting duplicate SNP: " $2 > "/dev/stderr"; exit 2} signatures[$2]=signature; if(!seen[$2]++) print $1,$2,$3,$4,$5,$6}' OFS='\t' "${weights[@]}" > "$tmp_weights"; then die "conflicting or invalid weight rows for $layer/$feature"; fi
    sort -k1,1n -k3,3n -k2,2 "$tmp_weights" > "$merged"
    rm -f -- "$tmp_weights"
    checksum="$(sha256sum "$merged" | awk '{print $1}')"
    printf '%s\t%s\n' "$feature" "$checksum" >> "$outdir/weight_checksums.tsv"
    n_weight="$(awk '!seen[$2]++ {n++} END{print n+0}' "$merged")"
    n_match="$(awk 'NR==FNR {b[$2]=1;next} b[$2]&&!seen[$2]++ {n++} END{print n+0}' "${bfile}.bim" "$merged")"
    [[ "$n_weight" -gt 0 ]] || die "empty weight file for $layer/$feature"
    rate="$(awk -v a="$n_match" -v b="$n_weight" 'BEGIN{printf "%.6f",a/b}')"
    awk -v r="$rate" -v m="$MIN_SNP_MATCH_RATE" 'BEGIN{exit !(r+0>=m+0)}' || die "SNP ID match rate $rate below $MIN_SNP_MATCH_RATE for $layer/$feature"
    chr_args=(); [[ "$chr" != ALL ]] && chr_args=(--chr "$chr")
    keep_args=(); [[ -n "${KEEP_FILE:-}" ]] && keep_args=(--keep "$(resolve_path "$KEEP_FILE")")
    log "PLINK score: $layer/$feature (SNP match=$rate)"
    "$PLINK_BIN" --bfile "$bfile" "${keep_args[@]}" "${chr_args[@]}" --out "$out" \
      --score "$merged" 2 4 6 sum --threads "$THREADS" >"$outdir/logs/$feature.stdout.log" 2>&1
    [[ -s "$out.profile" ]] || die "PLINK did not produce $out.profile"
    valid="$(awk '/valid predictor.*loaded/ {for(i=1;i<=NF;i++) if($i=="valid") print $(i-1)}' "$out.log" | tail -1)"
    [[ "$valid" =~ ^[0-9]+$ ]] || die "cannot parse valid predictor count from $out.log"
    rate="$(awk -v a="$valid" -v b="$n_weight" 'BEGIN{printf "%.6f",a/b}')"
    printf '%s\t%s\t%s\t%s\n' "$feature" "$n_weight" "$valid" "$rate" >> "$outdir/snp_qc.tsv"
    awk -v r="$rate" -v m="$MIN_SNP_MATCH_RATE" 'BEGIN{exit !(r+0>=m+0)}' || die "PLINK valid-predictor rate $rate below $MIN_SNP_MATCH_RATE for $layer/$feature"
  done < "$manifest"
done


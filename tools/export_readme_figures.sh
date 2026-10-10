#!/bin/sh
# Export the technical note's vector figures for Markdown display.
set -eu

project_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
figures_dir="$project_root/docs/technical_note/figures"
assets_dir="$project_root/docs/assets"

if ! command -v pdftoppm >/dev/null 2>&1; then
    printf '%s\n' "Missing pdftoppm: install Poppler to export the README figures." >&2
    exit 1
fi

for figure in transport convergence relaxation shear_wave hilker; do
    if [ ! -f "$figures_dir/$figure.pdf" ]; then
        printf 'Missing figure: %s\n' "$figures_dir/$figure.pdf" >&2
        exit 1
    fi
done

mkdir -p "$assets_dir"
for figure in transport convergence relaxation shear_wave hilker; do
    pdftoppm -f 1 -singlefile -png -scale-to-x 2400 -scale-to-y -1 \
        "$figures_dir/$figure.pdf" "$assets_dir/$figure"
    printf 'Exported docs/assets/%s.png\n' "$figure"
done

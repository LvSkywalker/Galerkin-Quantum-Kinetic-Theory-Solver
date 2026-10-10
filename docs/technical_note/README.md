# Technical note

`main.tex` builds the standalone project note. It includes the model, nondimensionalization, matched trial spaces, collision quadrature, transport verification, homogeneous relaxation and optional spatial applications.

The source is split into `computational_section.tex`, `dynamics_section.tex` and `hilker_section.tex`. `macros.tex` holds notation, `bibliography_entries.tex` contains references, `data/` contains exported tables and plot data, and `figures/` contains both ready-made vector PDFs and PGFPlots sources.

Compile this directory's `main.tex` with `latexmk -pdf main.tex` or `tectonic main.tex`. Overleaf can compile the whole directory as an uploaded project, with `main.tex` selected as the main document. No numerical simulation is run by LaTeX.

To refresh data from C++ results and rebuild figures, follow [the reproduction guide](../REPRODUCIBILITY.md#technical-note). The note is independently readable; no separate background document is required.

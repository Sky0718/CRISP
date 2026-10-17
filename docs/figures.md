# Figure generation

From the code repository root, install the plotting dependencies and supply the
complete extracted `Supplementary_Data` collection:

```bash
python -m pip install -r requirements_figures.txt
python crisp.py analysis figures --data ../Supplementary_Data --output ../CRISP_figures
```

This command renders all five main figures and four supplementary figures. Inputs
are numerical records from `Main_Figure_Reproduction`, `Figure_Source_Data` and
`Allocation_and_Graph_Controls/Composition`. Research data remain outside the
repository. Each invocation requires a new output directory.

## Output map

Paths below are relative to the chosen output directory.

| Manuscript figure | PDF output | Plotting source under `src/crisp/analyses/` |
| --- | --- | --- |
| Figure 1 | `figures/Fig_Overview/Fig_Overview.pdf` | `Main_Figure_Reproduction/build_figures.py` |
| Figure 2 | `figures/Figure_2.pdf` | `_figure_2.py` |
| Figure 3 | `figures/Fig_Candidates/Fig_Candidates.pdf` | `Main_Figure_Reproduction/build_figures.py` |
| Figure 4 | `figures/Fig_Composition.pdf` | `Variant_and_Network_Context/Figures/build_composition_figure.py` |
| Figure 5 | `figures/Fig_4/Fig_4.pdf` | `Main_Figure_Reproduction/build_figures.py` |
| Supplementary Figure 1 | `figures/Supplementary_Figure_1.pdf` | `_figure_s1.py` |
| Supplementary Figure 2 | `figures/Supplementary_Figure_2.pdf` | `_figure_s2.py` |
| Supplementary Figure 3 | `figures/Supplementary_Figure_3.pdf` | `_figure_s3.py` |
| Supplementary Figure 4 | `figures/Supplementary_Figure_4.pdf` | `_figure_s4.py` |

Main Figures 1, 2, 3 and 5 and all supplementary figures also have SVG and PNG
exports. Figure 4 has a PNG preview under `composition/`. The remaining output
files contain numerical sources and layout information. Figure captions belong in
the manuscript and SI, separately from the figure artwork.

## Selected supplementary panels

After installing the package, a subset of the JSON-based figures can be selected:

```bash
python -m crisp.analyses.figure_panels --data ../Supplementary_Data --output ../CRISP_selected_figures --figures 2 S1 S4
```

The available labels are `2`, `S1`, `S2`, `S3` and `S4`; `2` denotes main Figure 2.
For this module, `--data` can also point directly to `Figure_Source_Data`.

Arial is preferred, with an installed sans-serif fallback when necessary.
PDF fonts are embedded and SVG labels remain editable text. Preserve the supplied
data values, category order, axis definitions and panel dimensions when adapting
the output to another document format.

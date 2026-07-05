# SeeProteins

**[Live app →](https://zacharyardern.github.io/SeeProteins/)**

Browser-based protein structure alignment and visualisation. All computation runs locally in the browser via WebAssembly — no data is sent to any server.

## Features

- Rigid alignment via TM-align (TM-score, RMSD, sequence alignment panel)
- Flexible alignment with automatic hinge detection via BICExact
- Three-panel 3D viewer: overlay + individual structures
- Mutation highlighting (yellow Cα spheres on differing residues)
- Arrow-key rotation per selected panel
- Accepts PDB and mmCIF input

## Dependencies

| Tool | Purpose | Paper |
|------|---------|-------|
| [tmalign-wasm](https://github.com/milot-mirdita/tmalign-wasm) | TM-align compiled to WebAssembly | [Zhang & Skolnick (2005)](https://doi.org/10.1093/nar/gki524) |
| [BICExact](https://github.com/KoyanoBunsho/BICExact) | Hinge detection for flexible alignment | [Koyano & Shibuya (2025)](https://doi.org/10.1089/cmb.2024.0731) |
| [3Dmol.js](https://github.com/3dmol/3Dmol.js) | 3D molecular visualisation | [Rego & Koes (2015)](https://doi.org/10.1093/bioinformatics/btu829) |
| [Eigen](https://gitlab.com/libeigen/eigen) | Linear algebra (used within BICExact) | — |

## References

- Zhang Y, Skolnick J. TM-align: a protein structure alignment algorithm based on the TM-score. *Nucleic Acids Res.* 2005;33(8):2302–2309. doi:[10.1093/nar/gki524](https://doi.org/10.1093/nar/gki524)
- Koyano B, Shibuya T. Faster and More Accurate Estimation of Protein Hinges Based on Information Criteria. *J Comput Biol.* 2025;32(5):498–519. doi:[10.1089/cmb.2024.0731](https://doi.org/10.1089/cmb.2024.0731)
- Rego N, Koes D. 3Dmol.js: molecular visualization with WebGL. *Bioinformatics.* 2015;31(8):1322–1324. doi:[10.1093/bioinformatics/btu829](https://doi.org/10.1093/bioinformatics/btu829)
- Kabsch W. A solution for the best rotation to relate two sets of vectors. *Acta Crystallogr A.* 1976;32(5):922. doi:[10.1107/S0567739476001873](https://doi.org/10.1107/S0567739476001873)

## Citation

```bibtex
@software{ardern2026seeproteins,
  author  = {Ardern, Zachary},
  title   = {SeeProteins: browser-based flexible protein structure alignment and visualisation},
  year    = {2026},
  url     = {https://github.com/ZacharyArdern/SeeProteins}
}
```

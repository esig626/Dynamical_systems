# Archived modelling workspace

This repository is retained as a historical record. It was originally used as a mixed working archive and combined several unrelated modelling projects. Those projects have now been separated so that the public research record is easier to navigate.

## Salivary secretion work

The salivary fluid-secretion and cell-volume code formerly stored under `Volume_Model/` has been consolidated into [Salivary_Secretion](https://github.com/esig626/Salivary_Secretion), under `model/archive/volume_model_2022/`.

That destination contains the runnable standalone implementation and its small parameter files, preserved with their original Git blob hashes. The main salivary repository also contains the publication-era single-cell, reconstructed-cell, and multicellular models.

## Birmingham mitochondrial calcium work

The exploratory 2022 closed-cell mitochondrial Ca²⁺ model formerly stored under `Ca2+ Model/` has been preserved in the author's later metabolic-modelling archive. It is kept as precursor material and is not presented as part of the final SDHB model.

## Files intentionally retained here

Only two historical files remain in the current tree.

- `legacy/metabolism/model.m` is the original 2022 Tennant Lab kinetic cell-metabolism model with glycolysis, fatty-acid metabolism, TCA-cycle reactions, malate-aspartate-shuttle transport, redox balance, electron transport, oxidative phosphorylation, and free-radical-related modelling. It is unrelated to salivary secretion and is retained here byte-for-byte as historical source.
- `legacy/salivary/P.mat` is the original large MATLAB binary associated with the historical geometry-aware salivary volume-model workspace. It is retained here byte-for-byte rather than regenerated or converted during repository consolidation.

## Provenance

The complete pre-consolidation repository remains recoverable from Git history. The last commit containing the original mixed working tree is `d95cacfad8890a35e9a533e0e847921b128f8693`.

This repository is not an active software release and should not be used as the primary entry point for the salivary-secretion research programme. For that work, use [Salivary_Secretion](https://github.com/esig626/Salivary_Secretion).

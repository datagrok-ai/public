# Sources of the control-growth priors

The priors in `engine_priors_v*.json` are aggregates estimated from the vehicle/control arms of these open-access papers (PubMed Central): every paper with a control experiment in a line's base pool, which feeds the MAP priors, the calibration, the leave-one-out checks and the Monitor validation, and the papers of the treated arms behind the therapy archetypes (a paper can enter some of these and not the MAP). Per-paper results are not distributed: the `studies` table uses anonymous ids. The one exception is the example study (`example/`, CC BY 4.0, cited in `example/README.md`), distributed animal by animal under its licence.

`band_calibration_v1.json` holds no data from these papers: factors of the drawn bands, from the validation on papers left out of the fit (its provenance inside).

| DOI | Journal | Year | Licence | Lines |
|---|---|---|---|---|
| [10.1038/s41467-018-08115-8](https://doi.org/10.1038/s41467-018-08115-8) | Nature communications | 2019 | cc by | PDX-PDAC |
| [10.1038/s41467-020-15290-0](https://doi.org/10.1038/s41467-020-15290-0) | Nature communications | 2020 | cc by | MDA-MB-231 |
| [10.1038/s41467-020-17802-4](https://doi.org/10.1038/s41467-020-17802-4) | Nature communications | 2020 | cc by | MDA-MB-231 |
| [10.1038/s41467-020-18399-4](https://doi.org/10.1038/s41467-020-18399-4) | Nature communications | 2020 | cc by | PDX-PDAC |
| [10.1038/s41467-020-19920-5](https://doi.org/10.1038/s41467-020-19920-5) | Nature communications | 2020 | cc by | 4T1 |
| [10.1038/s41467-021-22676-1](https://doi.org/10.1038/s41467-021-22676-1) | Nature communications | 2021 | cc by | PDX-NSCLC |
| [10.1038/s41467-021-24099-4](https://doi.org/10.1038/s41467-021-24099-4) | Nature communications | 2021 | cc by | A549 |
| [10.1038/s41467-021-24798-y](https://doi.org/10.1038/s41467-021-24798-y) | Nature communications | 2021 | cc by | PDX-PDAC |
| [10.1038/s42003-021-02004-5](https://doi.org/10.1038/s42003-021-02004-5) | Communications biology | 2021 | cc by | B16F10 |
| [10.1038/s42003-021-02191-1](https://doi.org/10.1038/s42003-021-02191-1) | Communications biology | 2021 | cc by | CT26 |
| [10.1038/s42003-021-02741-7](https://doi.org/10.1038/s42003-021-02741-7) | Communications biology | 2021 | cc by | MDA-MB-231 |
| [10.1038/s41467-021-27936-8](https://doi.org/10.1038/s41467-021-27936-8) | Nature communications | 2022 | cc by | LLC |
| [10.1038/s41467-022-29444-9](https://doi.org/10.1038/s41467-022-29444-9) | Nature communications | 2022 | cc by | PDX-NSCLC |
| [10.1038/s41467-022-29606-9](https://doi.org/10.1038/s41467-022-29606-9) | Nature communications | 2022 | cc by | 4T1 |
| [10.1038/s41467-022-29938-6](https://doi.org/10.1038/s41467-022-29938-6) | Nature communications | 2022 | cc by | 4T1 |
| [10.1038/s41467-022-30149-2](https://doi.org/10.1038/s41467-022-30149-2) | Nature communications | 2022 | cc by | A549 |
| [10.1038/s41467-022-33430-6](https://doi.org/10.1038/s41467-022-33430-6) | Nature communications | 2022 | cc by | A549 |
| [10.1038/s41467-022-34559-0](https://doi.org/10.1038/s41467-022-34559-0) | Nature communications | 2022 | cc by | MDA-MB-231 |
| [10.1038/s41467-022-34643-5](https://doi.org/10.1038/s41467-022-34643-5) | Nature communications | 2022 | cc by | 4T1 |
| [10.1038/s43018-022-00393-y](https://doi.org/10.1038/s43018-022-00393-y) | Nature cancer | 2022 | cc by | PDX-PDAC |
| [10.1038/s43018-022-00438-2](https://doi.org/10.1038/s43018-022-00438-2) | Nature cancer | 2022 | cc by | SKOV3 |
| [10.1038/s43018-022-00470-2](https://doi.org/10.1038/s43018-022-00470-2) | Nature cancer | 2022 | cc by | Panc02 |
| [10.1038/s41467-023-37515-8](https://doi.org/10.1038/s41467-023-37515-8) | Nature communications | 2023 | cc by | A549, LLC |
| [10.1038/s41467-023-37656-w](https://doi.org/10.1038/s41467-023-37656-w) | Nature communications | 2023 | cc by | Panc02 |
| [10.1038/s41467-023-38128-x](https://doi.org/10.1038/s41467-023-38128-x) | Nature communications | 2023 | cc by | MDA-MB-231, Panc02 |
| [10.1038/s41467-023-38605-3](https://doi.org/10.1038/s41467-023-38605-3) | Nature communications | 2023 | cc by | LLC |
| [10.1038/s41467-023-38651-x](https://doi.org/10.1038/s41467-023-38651-x) | Nature communications | 2023 | cc by | A549 |
| [10.1038/s41467-023-39035-x](https://doi.org/10.1038/s41467-023-39035-x) | Nature communications | 2023 | cc by | Panc02 |
| [10.1038/s41467-023-40015-4](https://doi.org/10.1038/s41467-023-40015-4) | Nature communications | 2023 | cc by | 4T1 |
| [10.1038/s41467-023-40207-y](https://doi.org/10.1038/s41467-023-40207-y) | Nature communications | 2023 | cc by | Panc02 |
| [10.1038/s41467-023-40312-y](https://doi.org/10.1038/s41467-023-40312-y) | Nature communications | 2023 | cc by | LLC |
| [10.1038/s41467-023-40814-9](https://doi.org/10.1038/s41467-023-40814-9) | Nature communications | 2023 | cc by | LLC |
| [10.1038/s41467-023-40826-5](https://doi.org/10.1038/s41467-023-40826-5) | Nature communications | 2023 | cc by | 4T1 |
| [10.1038/s41467-023-41335-1](https://doi.org/10.1038/s41467-023-41335-1) | Nature communications | 2023 | cc by | 4T1 |
| [10.1038/s41467-023-41782-w](https://doi.org/10.1038/s41467-023-41782-w) | Nature communications | 2023 | cc by | LLC |
| [10.1038/s41467-023-42248-9](https://doi.org/10.1038/s41467-023-42248-9) | Nature communications | 2023 | cc by | 4T1 |
| [10.1038/s41467-023-43350-8](https://doi.org/10.1038/s41467-023-43350-8) | Nature communications | 2023 | cc by | Panc02 |
| [10.1038/s42003-023-04521-x](https://doi.org/10.1038/s42003-023-04521-x) | Communications biology | 2023 | cc by | Hepa1-6 |
| [10.1038/s42003-023-04983-z](https://doi.org/10.1038/s42003-023-04983-z) | Communications biology | 2023 | cc by | 4T1, Hepa1-6 |
| [10.1038/s42003-023-04993-x](https://doi.org/10.1038/s42003-023-04993-x) | Communications biology | 2023 | cc by | SKOV3 |
| [10.1038/s42003-023-05167-5](https://doi.org/10.1038/s42003-023-05167-5) | Communications biology | 2023 | cc by | MDA-MB-231 |
| [10.1038/s43018-023-00522-1](https://doi.org/10.1038/s43018-023-00522-1) | Nature cancer | 2023 | cc by | 4T1, LLC |
| [10.1038/s41467-023-43743-9](https://doi.org/10.1038/s41467-023-43743-9) | Nature communications | 2024 | cc by | Hepa1-6 |
| [10.1038/s41467-024-44779-1](https://doi.org/10.1038/s41467-024-44779-1) | Nature communications | 2024 | cc by | Panc02 |
| [10.1038/s41467-024-45340-w](https://doi.org/10.1038/s41467-024-45340-w) | Nature communications | 2024 | cc by | LLC |
| [10.1038/s41467-024-45698-x](https://doi.org/10.1038/s41467-024-45698-x) | Nature communications | 2024 | cc by | SKOV3 |
| [10.1038/s41467-024-45712-2](https://doi.org/10.1038/s41467-024-45712-2) | Nature communications | 2024 | cc by | Hepa1-6 |
| [10.1038/s41467-024-46778-8](https://doi.org/10.1038/s41467-024-46778-8) | Nature communications | 2024 | cc by | LLC |
| [10.1038/s41467-024-46811-w](https://doi.org/10.1038/s41467-024-46811-w) | Nature communications | 2024 | cc by | PDX-PDAC |
| [10.1038/s41467-024-47799-z](https://doi.org/10.1038/s41467-024-47799-z) | Nature communications | 2024 | cc by | SKOV3 |
| [10.1038/s41467-024-48364-4](https://doi.org/10.1038/s41467-024-48364-4) | Nature communications | 2024 | cc by | PDX-NSCLC |
| [10.1038/s41467-024-49573-7](https://doi.org/10.1038/s41467-024-49573-7) | Nature communications | 2024 | cc by | Hepa1-6 |
| [10.1038/s41467-024-49932-4](https://doi.org/10.1038/s41467-024-49932-4) | Nature communications | 2024 | cc by | SKOV3 |
| [10.1038/s41467-024-52170-3](https://doi.org/10.1038/s41467-024-52170-3) | Nature communications | 2024 | cc by | Hepa1-6 |
| [10.1038/s41467-024-53015-9](https://doi.org/10.1038/s41467-024-53015-9) | Nature communications | 2024 | cc by | CT26, MC38 |
| [10.1038/s41467-024-55006-2](https://doi.org/10.1038/s41467-024-55006-2) | Nature communications | 2024 | cc by | B16F10, MC38 (also the example study, `example/`: its controls are in the MC38 model) |
| [10.1038/s41467-024-55083-3](https://doi.org/10.1038/s41467-024-55083-3) | Nature communications | 2024 | cc by | 4T1, B16F10 |
| [10.1038/s42003-024-05854-x](https://doi.org/10.1038/s42003-024-05854-x) | Communications biology | 2024 | cc by | A549 |
| [10.1038/s42003-024-06271-w](https://doi.org/10.1038/s42003-024-06271-w) | Communications biology | 2024 | cc by | B16F10, CT26, MC38 |
| [10.1038/s42003-024-06402-3](https://doi.org/10.1038/s42003-024-06402-3) | Communications biology | 2024 | cc by | LLC |
| [10.1038/s42003-024-06839-6](https://doi.org/10.1038/s42003-024-06839-6) | Communications biology | 2024 | cc by | A549 |
| [10.1038/s43018-024-00754-9](https://doi.org/10.1038/s43018-024-00754-9) | Nature cancer | 2024 | cc by | PDX-NSCLC |
| [10.1038/s43018-024-00771-8](https://doi.org/10.1038/s43018-024-00771-8) | Nature cancer | 2024 | cc by | Panc02 |
| [10.1038/s41467-024-55768-9](https://doi.org/10.1038/s41467-024-55768-9) | Nature communications | 2025 | cc by | A549, LLC |
| [10.1038/s41467-025-57196-9](https://doi.org/10.1038/s41467-025-57196-9) | Nature communications | 2025 | cc by | B16F10 |
| [10.1038/s41467-025-58510-1](https://doi.org/10.1038/s41467-025-58510-1) | Nature communications | 2025 | cc by | B16F10 |
| [10.1038/s41467-025-58552-5](https://doi.org/10.1038/s41467-025-58552-5) | Nature communications | 2025 | cc by | 4T1 |
| [10.1038/s41467-025-59840-w](https://doi.org/10.1038/s41467-025-59840-w) | Nature communications | 2025 | cc by | B16F10, MC38 |
| [10.1038/s41467-025-60223-4](https://doi.org/10.1038/s41467-025-60223-4) | Nature communications | 2025 | cc by | A549 |
| [10.1038/s41467-025-60745-x](https://doi.org/10.1038/s41467-025-60745-x) | Nature communications | 2025 | cc by | SKOV3 |
| [10.1038/s41467-025-60768-4](https://doi.org/10.1038/s41467-025-60768-4) | Nature communications | 2025 | cc by | 4T1, MC38 |
| [10.1038/s41467-025-62295-8](https://doi.org/10.1038/s41467-025-62295-8) | Nature communications | 2025 | cc by | MDA-MB-231, PDX-PDAC |
| [10.1038/s41467-025-62999-x](https://doi.org/10.1038/s41467-025-62999-x) | Nature communications | 2025 | cc by | CT26 |
| [10.1038/s41467-025-63198-4](https://doi.org/10.1038/s41467-025-63198-4) | Nature communications | 2025 | cc by | 4T1, SKOV3 |
| [10.1038/s41467-025-63882-5](https://doi.org/10.1038/s41467-025-63882-5) | Nature communications | 2025 | cc by | MC38 |
| [10.1038/s41467-025-63892-3](https://doi.org/10.1038/s41467-025-63892-3) | Nature communications | 2025 | cc by | PDX-NSCLC |
| [10.1038/s41467-025-64081-y](https://doi.org/10.1038/s41467-025-64081-y) | Nature communications | 2025 | cc by | CT26, MC38 |
| [10.1038/s41467-025-67097-6](https://doi.org/10.1038/s41467-025-67097-6) | Nature communications | 2025 | cc by | 4T1 |
| [10.1038/s41467-025-67523-9](https://doi.org/10.1038/s41467-025-67523-9) | Nature communications | 2025 | cc by | MDA-MB-231 |
| [10.1038/s42003-025-09050-3](https://doi.org/10.1038/s42003-025-09050-3) | Communications biology | 2025 | cc by | 4T1 |
| [10.1038/s43018-024-00889-9](https://doi.org/10.1038/s43018-024-00889-9) | Nature cancer | 2025 | cc by | 4T1 |
| [10.1038/s41467-025-65876-9](https://doi.org/10.1038/s41467-025-65876-9) | Nature communications | 2026 | cc by | CT26 |
| [10.1038/s41467-026-68741-5](https://doi.org/10.1038/s41467-026-68741-5) | Nature communications | 2026 | cc by | MC38 |
| [10.1038/s41467-026-70215-7](https://doi.org/10.1038/s41467-026-70215-7) | Nature communications | 2026 | cc by | Panc02 |
| [10.1038/s41467-026-70949-4](https://doi.org/10.1038/s41467-026-70949-4) | Nature communications | 2026 | cc by | B16F10, CT26 |
| [10.1038/s41467-026-73300-z](https://doi.org/10.1038/s41467-026-73300-z) | Nature communications | 2026 | cc by | 4T1, LLC |
| [10.1038/s41467-026-73339-y](https://doi.org/10.1038/s41467-026-73339-y) | Nature communications | 2026 | cc by | B16F10 |
| [10.1038/s41467-026-73379-4](https://doi.org/10.1038/s41467-026-73379-4) | Nature communications | 2026 | cc by | B16F10, LLC |
| [10.1038/s41467-026-73788-5](https://doi.org/10.1038/s41467-026-73788-5) | Nature communications | 2026 | cc by | 4T1 |
| [10.1038/s41467-026-74022-y](https://doi.org/10.1038/s41467-026-74022-y) | Nature communications | 2026 | cc by | A549 |
| [10.1038/s41467-026-75247-7](https://doi.org/10.1038/s41467-026-75247-7) | Nature communications | 2026 | cc by | 4T1 |
| [10.1038/s41467-026-75300-5](https://doi.org/10.1038/s41467-026-75300-5) | Nature communications | 2026 | cc by | A549 |
| [10.1038/s41467-026-75852-6](https://doi.org/10.1038/s41467-026-75852-6) | Nature communications | 2026 | cc by | MC38 |
| [10.1038/s42003-026-10329-2](https://doi.org/10.1038/s42003-026-10329-2) | Communications biology | 2026 | cc by | MC38 |
| [10.1038/s43018-026-01156-9](https://doi.org/10.1038/s43018-026-01156-9) | Nature cancer | 2026 | cc by | Panc02 |
| [10.1038/s41467-024-51275-z](https://doi.org/10.1038/s41467-024-51275-z) | Nature communications | 2024 | cc by-nc-nd | LLC |
| [10.1038/s41467-024-51498-0](https://doi.org/10.1038/s41467-024-51498-0) | Nature communications | 2024 | cc by-nc-nd | MC38 |
| [10.1038/s41467-024-51534-z](https://doi.org/10.1038/s41467-024-51534-z) | Nature communications | 2024 | cc by-nc-nd | Panc02 |
| [10.1038/s41467-024-51945-y](https://doi.org/10.1038/s41467-024-51945-y) | Nature communications | 2024 | cc by-nc-nd | 4T1, CT26 |
| [10.1038/s41467-024-52202-y](https://doi.org/10.1038/s41467-024-52202-y) | Nature communications | 2024 | cc by-nc-nd | Hepa1-6 |
| [10.1038/s41467-024-53084-w](https://doi.org/10.1038/s41467-024-53084-w) | Nature communications | 2024 | cc by-nc-nd | 4T1, CT26 |
| [10.1038/s41467-024-53135-2](https://doi.org/10.1038/s41467-024-53135-2) | Nature communications | 2024 | cc by-nc-nd | 4T1, CT26 |
| [10.1038/s41467-024-53347-6](https://doi.org/10.1038/s41467-024-53347-6) | Nature communications | 2024 | cc by-nc-nd | CT26, MC38 |
| [10.1038/s41467-024-53392-1](https://doi.org/10.1038/s41467-024-53392-1) | Nature communications | 2024 | cc by-nc-nd | 4T1 |
| [10.1038/s41467-024-53706-3](https://doi.org/10.1038/s41467-024-53706-3) | Nature communications | 2024 | cc by-nc-nd | LLC |
| [10.1038/s41467-024-53771-8](https://doi.org/10.1038/s41467-024-53771-8) | Nature communications | 2024 | cc by-nc-nd | 4T1 |
| [10.1038/s41467-024-54081-9](https://doi.org/10.1038/s41467-024-54081-9) | Nature communications | 2024 | cc by-nc-nd | 4T1 |
| [10.1038/s41467-024-54292-0](https://doi.org/10.1038/s41467-024-54292-0) | Nature communications | 2024 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-024-54402-y](https://doi.org/10.1038/s41467-024-54402-y) | Nature communications | 2024 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-024-54650-y](https://doi.org/10.1038/s41467-024-54650-y) | Nature communications | 2024 | cc by-nc-nd | B16F10, CT26, LLC |
| [10.1038/s41467-024-54874-y](https://doi.org/10.1038/s41467-024-54874-y) | Nature communications | 2024 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-024-54877-9](https://doi.org/10.1038/s41467-024-54877-9) | Nature communications | 2024 | cc by-nc-nd | MC38 |
| [10.1038/s41467-024-54883-x](https://doi.org/10.1038/s41467-024-54883-x) | Nature communications | 2024 | cc by-nc-nd | CT26 |
| [10.1038/s41467-024-55207-9](https://doi.org/10.1038/s41467-024-55207-9) | Nature communications | 2024 | cc by-nc-nd | MDA-MB-231 |
| [10.1038/s41467-024-54132-1](https://doi.org/10.1038/s41467-024-54132-1) | Nature communications | 2025 | cc by-nc-nd | CT26 |
| [10.1038/s41467-025-55876-0](https://doi.org/10.1038/s41467-025-55876-0) | Nature communications | 2025 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-025-56137-w](https://doi.org/10.1038/s41467-025-56137-w) | Nature communications | 2025 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-025-56507-4](https://doi.org/10.1038/s41467-025-56507-4) | Nature communications | 2025 | cc by-nc-nd | MDA-MB-231 |
| [10.1038/s41467-025-56648-6](https://doi.org/10.1038/s41467-025-56648-6) | Nature communications | 2025 | cc by-nc-nd | MC38 |
| [10.1038/s41467-025-57191-0](https://doi.org/10.1038/s41467-025-57191-0) | Nature communications | 2025 | cc by-nc-nd | 4T1 |
| [10.1038/s41467-025-57364-x](https://doi.org/10.1038/s41467-025-57364-x) | Nature communications | 2025 | cc by-nc-nd | LLC, PDX-NSCLC |
| [10.1038/s41467-025-57533-y](https://doi.org/10.1038/s41467-025-57533-y) | Nature communications | 2025 | cc by-nc-nd | CT26 |
| [10.1038/s41467-025-57603-1](https://doi.org/10.1038/s41467-025-57603-1) | Nature communications | 2025 | cc by-nc-nd | LLC |
| [10.1038/s41467-025-58740-3](https://doi.org/10.1038/s41467-025-58740-3) | Nature communications | 2025 | cc by-nc-nd | 4T1 |
| [10.1038/s41467-025-59193-4](https://doi.org/10.1038/s41467-025-59193-4) | Nature communications | 2025 | cc by-nc-nd | MC38 |
| [10.1038/s41467-025-59490-y](https://doi.org/10.1038/s41467-025-59490-y) | Nature communications | 2025 | cc by-nc-nd | Panc02 |
| [10.1038/s41467-025-59621-5](https://doi.org/10.1038/s41467-025-59621-5) | Nature communications | 2025 | cc by-nc-nd | LLC |
| [10.1038/s41467-025-59874-0](https://doi.org/10.1038/s41467-025-59874-0) | Nature communications | 2025 | cc by-nc-nd | CT26, MC38 |
| [10.1038/s41467-025-60463-4](https://doi.org/10.1038/s41467-025-60463-4) | Nature communications | 2025 | cc by-nc-nd | B16F10, MC38 |
| [10.1038/s41467-025-60769-3](https://doi.org/10.1038/s41467-025-60769-3) | Nature communications | 2025 | cc by-nc-nd | B16F10, MC38 |
| [10.1038/s41467-025-60927-7](https://doi.org/10.1038/s41467-025-60927-7) | Nature communications | 2025 | cc by-nc-nd | 4T1, CT26, MC38 |
| [10.1038/s41467-025-60958-0](https://doi.org/10.1038/s41467-025-60958-0) | Nature communications | 2025 | cc by-nc-nd | CT26, MC38 |
| [10.1038/s41467-025-61299-8](https://doi.org/10.1038/s41467-025-61299-8) | Nature communications | 2025 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-025-61439-0](https://doi.org/10.1038/s41467-025-61439-0) | Nature communications | 2025 | cc by-nc-nd | B16F10, CT26, MC38 |
| [10.1038/s41467-025-62132-y](https://doi.org/10.1038/s41467-025-62132-y) | Nature communications | 2025 | cc by-nc-nd | Hepa1-6 |
| [10.1038/s41467-025-62134-w](https://doi.org/10.1038/s41467-025-62134-w) | Nature communications | 2025 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-025-62558-4](https://doi.org/10.1038/s41467-025-62558-4) | Nature communications | 2025 | cc by-nc-nd | MC38 |
| [10.1038/s41467-025-62776-w](https://doi.org/10.1038/s41467-025-62776-w) | Nature communications | 2025 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-025-63311-7](https://doi.org/10.1038/s41467-025-63311-7) | Nature communications | 2025 | cc by-nc-nd | Hepa1-6 |
| [10.1038/s41467-025-63338-w](https://doi.org/10.1038/s41467-025-63338-w) | Nature communications | 2025 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-025-63659-w](https://doi.org/10.1038/s41467-025-63659-w) | Nature communications | 2025 | cc by-nc-nd | MC38 |
| [10.1038/s41467-025-63687-6](https://doi.org/10.1038/s41467-025-63687-6) | Nature communications | 2025 | cc by-nc-nd | A549, LLC |
| [10.1038/s41467-025-64149-9](https://doi.org/10.1038/s41467-025-64149-9) | Nature communications | 2025 | cc by-nc-nd | 4T1 |
| [10.1038/s41467-025-64572-y](https://doi.org/10.1038/s41467-025-64572-y) | Nature communications | 2025 | cc by-nc-nd | Hepa1-6 |
| [10.1038/s41467-025-64675-6](https://doi.org/10.1038/s41467-025-64675-6) | Nature communications | 2025 | cc by-nc-nd | PDX-NSCLC |
| [10.1038/s41467-025-65025-2](https://doi.org/10.1038/s41467-025-65025-2) | Nature communications | 2025 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-025-65048-9](https://doi.org/10.1038/s41467-025-65048-9) | Nature communications | 2025 | cc by-nc-nd | MC38 |
| [10.1038/s41467-025-65762-4](https://doi.org/10.1038/s41467-025-65762-4) | Nature communications | 2025 | cc by-nc-nd | PDX-PDAC |
| [10.1038/s41467-025-65927-1](https://doi.org/10.1038/s41467-025-65927-1) | Nature communications | 2025 | cc by-nc-nd | Hepa1-6 |
| [10.1038/s41467-025-65945-z](https://doi.org/10.1038/s41467-025-65945-z) | Nature communications | 2025 | cc by-nc-nd | Hepa1-6 |
| [10.1038/s41467-025-66089-w](https://doi.org/10.1038/s41467-025-66089-w) | Nature communications | 2025 | cc by-nc-nd | CT26 |
| [10.1038/s41467-025-66168-y](https://doi.org/10.1038/s41467-025-66168-y) | Nature communications | 2025 | cc by-nc-nd | Hepa1-6 |
| [10.1038/s41467-025-66266-x](https://doi.org/10.1038/s41467-025-66266-x) | Nature communications | 2025 | cc by-nc-nd | MDA-MB-231 |
| [10.1038/s41467-025-66630-x](https://doi.org/10.1038/s41467-025-66630-x) | Nature communications | 2025 | cc by-nc-nd | B16F10, MC38, MDA-MB-231 |
| [10.1038/s41467-025-67254-x](https://doi.org/10.1038/s41467-025-67254-x) | Nature communications | 2025 | cc by-nc-nd | Panc02 |
| [10.1038/s41467-025-67356-6](https://doi.org/10.1038/s41467-025-67356-6) | Nature communications | 2025 | cc by-nc-nd | MC38 |
| [10.1038/s41467-025-67365-5](https://doi.org/10.1038/s41467-025-67365-5) | Nature communications | 2025 | cc by-nc-nd | CT26, MC38 |
| [10.1038/s41467-025-67437-6](https://doi.org/10.1038/s41467-025-67437-6) | Nature communications | 2025 | cc by-nc-nd | MC38 |
| [10.1038/s41467-025-67674-9](https://doi.org/10.1038/s41467-025-67674-9) | Nature communications | 2025 | cc by-nc-nd | CT26 |
| [10.1038/s41467-025-67837-8](https://doi.org/10.1038/s41467-025-67837-8) | Nature communications | 2025 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-025-67967-z](https://doi.org/10.1038/s41467-025-67967-z) | Nature communications | 2025 | cc by-nc-nd | MC38 |
| [10.1038/s41586-025-08974-4](https://doi.org/10.1038/s41586-025-08974-4) | Nature | 2025 | cc by-nc-nd | 4T1 |
| [10.1038/s42003-025-07452-x](https://doi.org/10.1038/s42003-025-07452-x) | Communications biology | 2025 | cc by-nc-nd | A549 |
| [10.1038/s42003-025-07897-0](https://doi.org/10.1038/s42003-025-07897-0) | Communications biology | 2025 | cc by-nc-nd | Panc02 |
| [10.1038/s42003-025-08108-6](https://doi.org/10.1038/s42003-025-08108-6) | Communications biology | 2025 | cc by-nc-nd | LLC |
| [10.1038/s42003-025-08520-y](https://doi.org/10.1038/s42003-025-08520-y) | Communications biology | 2025 | cc by-nc-nd | A549 |
| [10.1038/s42003-025-09071-y](https://doi.org/10.1038/s42003-025-09071-y) | Communications biology | 2025 | cc by-nc-nd | A549 |
| [10.1038/s42003-025-09331-x](https://doi.org/10.1038/s42003-025-09331-x) | Communications biology | 2025 | cc by-nc-nd | SKOV3 |
| [10.1038/s41467-025-67984-y](https://doi.org/10.1038/s41467-025-67984-y) | Nature communications | 2026 | cc by-nc-nd | B16F10, MC38 |
| [10.1038/s41467-025-68152-y](https://doi.org/10.1038/s41467-025-68152-y) | Nature communications | 2026 | cc by-nc-nd | MDA-MB-231 |
| [10.1038/s41467-025-68279-y](https://doi.org/10.1038/s41467-025-68279-y) | Nature communications | 2026 | cc by-nc-nd | 4T1, B16F10, CT26 |
| [10.1038/s41467-026-68466-5](https://doi.org/10.1038/s41467-026-68466-5) | Nature communications | 2026 | cc by-nc-nd | CT26, MC38 |
| [10.1038/s41467-026-68490-5](https://doi.org/10.1038/s41467-026-68490-5) | Nature communications | 2026 | cc by-nc-nd | Panc02 |
| [10.1038/s41467-026-68502-4](https://doi.org/10.1038/s41467-026-68502-4) | Nature communications | 2026 | cc by-nc-nd | 4T1 |
| [10.1038/s41467-026-68525-x](https://doi.org/10.1038/s41467-026-68525-x) | Nature communications | 2026 | cc by-nc-nd | CT26, MC38 |
| [10.1038/s41467-026-68649-0](https://doi.org/10.1038/s41467-026-68649-0) | Nature communications | 2026 | cc by-nc-nd | B16F10, CT26, Hepa1-6, MC38 |
| [10.1038/s41467-026-68661-4](https://doi.org/10.1038/s41467-026-68661-4) | Nature communications | 2026 | cc by-nc-nd | B16F10, MC38 |
| [10.1038/s41467-026-68765-x](https://doi.org/10.1038/s41467-026-68765-x) | Nature communications | 2026 | cc by-nc-nd | LLC |
| [10.1038/s41467-026-68834-1](https://doi.org/10.1038/s41467-026-68834-1) | Nature communications | 2026 | cc by-nc-nd | MC38 |
| [10.1038/s41467-026-69305-3](https://doi.org/10.1038/s41467-026-69305-3) | Nature communications | 2026 | cc by-nc-nd | 4T1, CT26, LLC |
| [10.1038/s41467-026-69311-5](https://doi.org/10.1038/s41467-026-69311-5) | Nature communications | 2026 | cc by-nc-nd | Panc02 |
| [10.1038/s41467-026-69339-7](https://doi.org/10.1038/s41467-026-69339-7) | Nature communications | 2026 | cc by-nc-nd | MC38 |
| [10.1038/s41467-026-69456-3](https://doi.org/10.1038/s41467-026-69456-3) | Nature communications | 2026 | cc by-nc-nd | 4T1, CT26, MC38 |
| [10.1038/s41467-026-69462-5](https://doi.org/10.1038/s41467-026-69462-5) | Nature communications | 2026 | cc by-nc-nd | CT26, MC38 |
| [10.1038/s41467-026-69695-4](https://doi.org/10.1038/s41467-026-69695-4) | Nature communications | 2026 | cc by-nc-nd | PDX-PDAC |
| [10.1038/s41467-026-69796-0](https://doi.org/10.1038/s41467-026-69796-0) | Nature communications | 2026 | cc by-nc-nd | Panc02 |
| [10.1038/s41467-026-70657-z](https://doi.org/10.1038/s41467-026-70657-z) | Nature communications | 2026 | cc by-nc-nd | CT26, MC38 |
| [10.1038/s41467-026-70795-4](https://doi.org/10.1038/s41467-026-70795-4) | Nature communications | 2026 | cc by-nc-nd | 4T1 |
| [10.1038/s41467-026-70859-5](https://doi.org/10.1038/s41467-026-70859-5) | Nature communications | 2026 | cc by-nc-nd | 4T1, CT26 |
| [10.1038/s41467-026-70867-5](https://doi.org/10.1038/s41467-026-70867-5) | Nature communications | 2026 | cc by-nc-nd | MC38 |
| [10.1038/s41467-026-70876-4](https://doi.org/10.1038/s41467-026-70876-4) | Nature communications | 2026 | cc by-nc-nd | B16F10, MC38 |
| [10.1038/s41467-026-71363-6](https://doi.org/10.1038/s41467-026-71363-6) | Nature communications | 2026 | cc by-nc-nd | MC38 |
| [10.1038/s41467-026-71476-y](https://doi.org/10.1038/s41467-026-71476-y) | Nature communications | 2026 | cc by-nc-nd | 4T1 |
| [10.1038/s41467-026-71739-8](https://doi.org/10.1038/s41467-026-71739-8) | Nature communications | 2026 | cc by-nc-nd | LLC |
| [10.1038/s41467-026-72062-y](https://doi.org/10.1038/s41467-026-72062-y) | Nature communications | 2026 | cc by-nc-nd | A549 |
| [10.1038/s41467-026-72814-w](https://doi.org/10.1038/s41467-026-72814-w) | Nature communications | 2026 | cc by-nc-nd | Hepa1-6 |
| [10.1038/s41467-026-72973-w](https://doi.org/10.1038/s41467-026-72973-w) | Nature communications | 2026 | cc by-nc-nd | MC38 |
| [10.1038/s41467-026-72978-5](https://doi.org/10.1038/s41467-026-72978-5) | Nature communications | 2026 | cc by-nc-nd | 4T1, CT26, MDA-MB-231 |
| [10.1038/s41467-026-73109-w](https://doi.org/10.1038/s41467-026-73109-w) | Nature communications | 2026 | cc by-nc-nd | Panc02 |
| [10.1038/s41467-026-73228-4](https://doi.org/10.1038/s41467-026-73228-4) | Nature communications | 2026 | cc by-nc-nd | Hepa1-6 |
| [10.1038/s41467-026-73234-6](https://doi.org/10.1038/s41467-026-73234-6) | Nature communications | 2026 | cc by-nc-nd | MC38 |
| [10.1038/s41467-026-73380-x](https://doi.org/10.1038/s41467-026-73380-x) | Nature communications | 2026 | cc by-nc-nd | PDX-PDAC |
| [10.1038/s41467-026-73401-9](https://doi.org/10.1038/s41467-026-73401-9) | Nature communications | 2026 | cc by-nc-nd | 4T1 |
| [10.1038/s41467-026-73683-z](https://doi.org/10.1038/s41467-026-73683-z) | Nature communications | 2026 | cc by-nc-nd | CT26, LLC |
| [10.1038/s41467-026-73942-z](https://doi.org/10.1038/s41467-026-73942-z) | Nature communications | 2026 | cc by-nc-nd | CT26, LLC, MC38 |
| [10.1038/s41467-026-74162-1](https://doi.org/10.1038/s41467-026-74162-1) | Nature communications | 2026 | cc by-nc-nd | LLC |
| [10.1038/s41467-026-74352-x](https://doi.org/10.1038/s41467-026-74352-x) | Nature communications | 2026 | cc by-nc-nd | CT26 |
| [10.1038/s41467-026-74891-3](https://doi.org/10.1038/s41467-026-74891-3) | Nature communications | 2026 | cc by-nc-nd | B16F10 |
| [10.1038/s41467-026-75132-3](https://doi.org/10.1038/s41467-026-75132-3) | Nature communications | 2026 | cc by-nc-nd | PDX-PDAC |
| [10.1038/s41467-026-75208-0](https://doi.org/10.1038/s41467-026-75208-0) | Nature communications | 2026 | cc by-nc-nd | CT26, MC38 |
| [10.1038/s41467-026-76158-3](https://doi.org/10.1038/s41467-026-76158-3) | Nature communications | 2026 | cc by-nc-nd | 4T1, Panc02 |
| [10.1038/s41467-026-76239-3](https://doi.org/10.1038/s41467-026-76239-3) | Nature communications | 2026 | cc by-nc-nd | CT26 |
| [10.1038/s41467-026-76647-5](https://doi.org/10.1038/s41467-026-76647-5) | Nature communications | 2026 | cc by-nc-nd | 4T1, B16F10, LLC |
| [10.1038/s41467-026-77547-4](https://doi.org/10.1038/s41467-026-77547-4) | Nature communications | 2026 | cc by-nc-nd | 4T1, B16F10 |
| [10.1038/s41467-026-77904-3](https://doi.org/10.1038/s41467-026-77904-3) | Nature communications | 2026 | cc by-nc-nd | 4T1 |
| [10.1038/s41467-026-78085-9](https://doi.org/10.1038/s41467-026-78085-9) | Nature communications | 2026 | cc by-nc-nd | 4T1 |
| [10.1038/s42003-026-10070-w](https://doi.org/10.1038/s42003-026-10070-w) | Communications biology | 2026 | cc by-nc-nd | LLC |
| [10.1038/s42003-026-10424-4](https://doi.org/10.1038/s42003-026-10424-4) | Communications biology | 2026 | cc by-nc-nd | B16F10 |
| [10.1038/s42003-026-10433-3](https://doi.org/10.1038/s42003-026-10433-3) | Communications biology | 2026 | cc by-nc-nd | LLC |

# Sources of the tumor density

`tumour_density_v0.json` (the tumor-corrected body weight's density per line, mg per mm³ of caliper volume) is an aggregate of the excised tumors' weights over their last caliper volume in these open-access papers' per-animal source data (murine subcutaneous tumors of 200 mm³ or more). Per-animal values are not distributed.

| Paper | Licence | Lines |
|---|---|---|
| [10.1038/s41467-024-53135-2](https://doi.org/10.1038/s41467-024-53135-2) | cc by-nc-nd | CT26 |
| [10.1038/s41467-025-55876-0](https://doi.org/10.1038/s41467-025-55876-0) | cc by-nc-nd | B16F10 |
| [10.1038/s41467-025-59193-4](https://doi.org/10.1038/s41467-025-59193-4) | cc by-nc-nd | MC38, other murine |
| [10.1038/s41467-025-58552-5](https://doi.org/10.1038/s41467-025-58552-5) | cc by | other murine |
| [10.1038/s41467-025-60958-0](https://doi.org/10.1038/s41467-025-60958-0) | cc by-nc-nd | CT26, MC38 |
| [10.1038/s41467-025-62134-w](https://doi.org/10.1038/s41467-025-62134-w) | cc by-nc-nd | B16F10 |
| [10.1038/s41467-025-63338-w](https://doi.org/10.1038/s41467-025-63338-w) | cc by-nc-nd | B16F10 |
| [10.1038/s41467-025-67315-1](https://doi.org/10.1038/s41467-025-67315-1) | cc by-nc-nd | other murine |
| [10.1038/s41467-025-67356-6](https://doi.org/10.1038/s41467-025-67356-6) | cc by-nc-nd | MC38 |
| [10.1038/s41467-025-67437-6](https://doi.org/10.1038/s41467-025-67437-6) | cc by-nc-nd | MC38 |
| [10.1038/s41467-026-68525-x](https://doi.org/10.1038/s41467-026-68525-x) | cc by-nc-nd | CT26, MC38 |
| [10.1038/s41467-026-68765-x](https://doi.org/10.1038/s41467-026-68765-x) | cc by-nc-nd | MC38, other murine |
| [10.1038/s41467-026-68741-5](https://doi.org/10.1038/s41467-026-68741-5) | cc by | MC38 |
| [10.1038/s41467-026-69350-y](https://doi.org/10.1038/s41467-026-69350-y) | cc by-nc-nd | other murine |
| [10.1038/s41467-026-72847-1](https://doi.org/10.1038/s41467-026-72847-1) | cc by-nc-nd | B16F10 |
| [10.1038/s41467-026-73942-z](https://doi.org/10.1038/s41467-026-73942-z) | cc by-nc-nd | CT26 |
| [PMC13650124](https://pmc.ncbi.nlm.nih.gov/articles/PMC13650124/) | cc by-nc-nd | B16F10, other murine |
| [10.1038/s41467-026-77495-w](https://doi.org/10.1038/s41467-026-77495-w) | cc by-nc-nd | CT26, MC38 |

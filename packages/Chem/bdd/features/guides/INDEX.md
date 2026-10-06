# Questions answered as guides

One line per question: the question as asked, the feature that answers it, when, and the help
page it illustrates (its `@help:` tag) if any. `grok-bdd guide <feature>` films it; the `/bdd-answer`
skill is the recipe.

| Question | Feature | Date | Help page |
|---|---|---|---|
| How do I split my series into the common core and its substituents, and see which combinations of substituents I have? | r-groups-from-mcs.feature | 2026-09-24 | datagrok/solutions/domains/chem |
| How do I run a SAR Matrix analysis? | sar-matrix.feature | 2026-10-06 | datagrok/solutions/domains/chem |
| My compounds already come in parts — can the SAR Matrix use them? | sar-matrix-protac.feature | 2026-10-06 | datagrok/solutions/domains/chem |

## The two SAR Matrix films

`sar-matrix.feature` hands the analysis structures only, so it cuts them itself and the Summary's
answers are per series; `sar-matrix-protac.feature` hands it a core column and two R-group columns,
so nothing is cut, one fit ranks all three parts and each part's measured pairs are pooled
separately.

Both carry the caption strip and no title card. The video is 1600×800; filming the page at 1600×736
leaves exactly the 64 px the strip needs, so the composed frame lands on the output size without
being resampled.

```
cd packages/Chem/bdd
BDD_GUIDE_VIEWPORT=1600x736 BDD_TASKBAR_MS=1800000 \
  npx grok-bdd guide features/guides/sar-matrix.feature --timeout=1800000 \
  --render "--no-title --size 1600x800 --hold 1.6 --travel 0.6"
```

`BDD_TASKBAR_MS` and `--timeout` are only needed while the matrix build on the 10,000-row demo set
takes minutes. The GIF goes in the `help` submodule, which is a checkout of its own.

The film reads the demo set on `hERG_pIC50`, higher-is-better, scaling none. The dialog defaults to
`CYP3A4`, which with no scaling is read as a raw assay value — and a raw column has no log scale, so
every fitted effect on the Summary refuses to rank rather than claiming a fold it cannot support.
The film sets the activity and the direction so the walkthrough shows findings rather than that
refusal.

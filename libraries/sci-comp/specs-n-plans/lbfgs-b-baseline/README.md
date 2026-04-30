# L-BFGS-B Remediation: Baseline Snapshot

Raw stdout captures of the three benchmark runners taken **before** any of the
remediation commits in [`../lbfgs-b-remediation-plan.md`](../lbfgs-b-remediation-plan.md)
were applied. Use these as the reference point for spotting regressions and
confirming improvements as commits 1–8 land.

## Files

| File | Source command |
|---|---|
| `unconstrained.txt` | `npx tsx src/optimization/single-objective/benchmarks/unconstrained-benchmarks.ts` |
| `multistart.txt`    | `npx tsx src/optimization/single-objective/benchmarks/multistart-benchmarks.ts` |
| `bounded.txt`       | `npx tsx src/optimization/single-objective/benchmarks/bounded-benchmarks.ts` |

Captured from `master` at the commit immediately after
`GROK-20018: Sci Comp: Add L-BFGS-B remediation plan` (i.e., L-BFGS-B as
implemented in `GROK-20018: Sci Comp: Implement L-BFGS-B optimizer`, plain).

## Which commit checks which file

Per the plan's commit-to-runner mapping:

| Commit | Runner(s) to re-check after the commit |
|---|---|
| 1 | `unconstrained` (NaN-cap should not affect healthy problems) |
| 2 | `unconstrained`, `multistart` (different α₀ → different traces) |
| 3 | all three (last BFGS update skipped → −1 iter on converged tasks) |
| 4 | all three (curvature gate redraws `nskip`) |
| 5 | all three (snapshot is a no-op functionally) |
| 6 | `bounded` (SD fallback removed; previously-passing cases may flip) |
| 7 | `bounded` (false success at `stpMax = 0` corrected) |
| 8 | `bounded` — main target; expected significant improvement |

## Headline numbers from the baseline

`bounded.txt` — Success summary (feasible + accurate):

```
L-BFGS-B              7/7
L-BFGS                3/7
PSO                   3/7
Adam                  2/7
GradientDescent       2/7
Nelder-Mead           2/7
```

The L-BFGS-B 7/7 result is the bar to maintain or improve. Note that some of
those wins may be coming through the false-success paths (`stpMax = 0`,
projected-SD fallback) that commits 6 and 7 close — if 7/7 dips after those
commits, that is the masked-bug signal the plan warns about, and is expected
to be recovered by commit 8.

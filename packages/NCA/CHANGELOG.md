# Nca changelog

## v.next

* Hub demo studies: five simulated studies of one discovery program (Series X, lead DGK-0007), installed with the new `ncaSeedDemoStudies()` function, which any user can run — **once per server, and again by the demo's owner or an administrator after every sci-comp update** (it re-creates out-of-date demo studies). See README, *Demo data → Demo studies*.
* Studies and analyses are stored in the platform database (EMS). The analysed table is kept as a snapshot with each analysis.
* Studies saved before this version are not carried over.
* Studies are private to their creator until shared. Sharing grants View access; a View-only user can open and explore a study but not change it.
* Deleting a study or an analysis is final: its data snapshot is removed and it cannot be restored.
* Results, λz point overrides, review decisions, pipeline flags and pipeline runs are stored as records of their own instead of inside the analysis. Each results row keeps Mean, %CV and n per parameter as numbers, so results can be sorted and searched across studies. The audit trail shows each λz point and each review decision as a separate record with its analyst and time. A superseded review decision stays stored next to the current one. Each λz point and each review decision records the analyst as a platform user, not only by name, so a renamed user or two users with the same name stay distinguishable.
* The AUC method, the four BLQ rules and the sparse sampling design of each analysis are stored as columns, so a change of BLQ rule shows in the audit trail as that rule.
* Opening an analysis no longer rewrites it when its results did not change: the analysis keeps its version and the audit trail gets no record. The "Last computed" date of a study's analyses now shows when the results last changed.
* The study overview reads its results without the concentration-time curves and loads the curves just after it opens.
* Fixed: in Sparse NCA a composite row could show AUCinf smaller than AUClast (mouse sparse F: 2.173 < 2.312). AUCinf, %AUCextrap, CL, Vz, MRT and Vss of a composite now rest on the same linear-trapezoid AUClast that the row shows and that its Bailer / Holder SE describes, as PKNCA computes them on the composite mean. These values change on every sparse analysis, so review decisions recorded on its composites are reset the next time it is computed. On an IV-bolus composite the derived values still start from the back-extrapolated C0 while AUClast starts from zero; the row now says so and gives both AUClast values.
* Fixed: after a λz point edit, Save stored the profile's warnings from before the edit (the study grid's Warnings and the analysis status were off until the analysis was reopened).
* Fixed: after the first Save of a new analysis the address bar read `/apps/Nca/apps/Nca/analyses/…`, and reloading the page opened the Hub instead of the analysis.
* Fixed: an analysis grouped by a results column such as CL basis had no concentration-time curves in its study (the context panel showed "No measurable concentration–time points to plot"). Each group now has its own curve; an analysis computed before the fix gets its curves the next time it is opened.
* CL, Vz and Vss are reported per kg of body weight everywhere: CL in mL/min/kg, Vz and Vss in L/kg, with the unit in the column header. The absolute values are in the new CL abs. (L/h), Vz abs. (L) and Vss abs. (L) columns, with BW (kg), in the Full profile. Without a body weight, a mg/kg dose still gives per-kg values (no absolute ones), and an absolute dose gives only the absolute ones.
* Sparse composites report CL, Vz and Vss per kg (dose per kg ÷ AUC of the mean curve), also when no body weight is mapped.
* The projected exposure margin reads the per-kg CL and uses the cohort's own mean body weight as the rat reference weight when the data carry one (shown in the new Rat ref. BW columns).
* Studies computed before this version show "stale (refreshes on open)", and their CL, Vz and Vss cells stay empty until the analysis is reopened.
* The confirmation before computing a mg/kg dose without a body weight was removed.
* `ncaCompute` output contract: `CL` / `VZ` / `VSS` now mean per kg (mL/min/kg, L/kg) on every row, and the new `CL_ABS` / `VZ_ABS` / `VSS_ABS` columns carry the absolute values. `CL_ML_MIN_KG` and `VZ_L_KG` are kept for this release as deprecated aliases of `CL` and `VZ` (column tag `nca.deprecated`) and will be removed in the next one. A `RESULTS_BASIS` column and the `nca.resultsBasis` frame tag (`per-kg-1`) mark the new contract. `BW_KG` is now always in kg (a body weight mapped in grams was divided 1000× too little), and the CL unit tag is correct on a table that mixes subjects with and without a body weight.

## 0.0.1 (2026-05-12)

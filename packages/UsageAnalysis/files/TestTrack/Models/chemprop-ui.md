---
feature: models
target_layer: manual-only
coverage_type: regression
manual_only_reason: |
  Chemprop training runs multi-minute compute in the chem-chemprop Docker
  container whose startup state is not reliable in CI, and the final
  prediction-vs-actual check is a visual correlation judgment on a scatter
  plot.
related_bugs: []
---

# Chemprop — manual checks

## Docker container lifecycle

1. Go to **Browse > Platform > Dockers** and locate the container
   named **chem-chemprop**.
   **Verify:** the container card surfaces with its current
   status (Running or Stopped).
2. Right-click the **chem-chemprop** container and select
   **Stop**.
   **Verify:** the container transitions to Stopped (the card
   status updates; no error surface).
3. Once the container has stopped, right-click it again and select
   **Run** to restart it.
   **Verify:** the container transitions back to Running within
   ~2 minutes and the Dockers view lists it as Running.
4. Cleanup: if the run aborts mid-scenario, make sure the shared
   `chem-chemprop` container is restored to Running.

## Train, apply, and verify prediction–activity correlation

Precondition: the `chem-chemprop` Docker container is running.

1. Close all and open `System:AppData/Chem/mol1K.sdf` — the table renders with
   the molecule column and a numeric `pIC50_HIV_Integrase` column.
2. From the top menu, run **ML > Models > Train Model...** — the Train Model
   view opens and the **Chemprop** trainer is selectable for the molecule
   column.
3. Set the molecule column as **features** and `pIC50_HIV_Integrase` as the
   **predict** target; choose the Chemprop trainer. No console errors.
4. Run training (can take minutes), naming the model `test_chemprop_manual`;
   once the model is built, apply it back to the same dataset — a new
   prediction column is appended to the table, no console errors.
5. Add a scatter plot comparing the prediction column against
   `pIC50_HIV_Integrase` — the points cluster around the y = x diagonal
   (prediction nearly equal to actual). This is the end-to-end correctness
   check of the train + apply round-trip.
6. Cleanup: delete the `test_chemprop_manual` model and close the view.

---
{
  "order": 1,
  "datasets": ["System:AppData/Chem/mol1K.sdf"]
}

# Chemprop Docker container lifecycle

Verifies that the `chem-chemprop` Docker container backing Chemprop
training/apply can be stopped and restarted from
**Browse > Platform > Dockers** without breaking the model
workflow — training and applying a Chemprop model still works once
the container comes back up.

## Block 3: Manage the chem-chemprop Docker container

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
   **Verify:** the container transitions back to Running; the
   restart completes within a bounded wait; downstream Chemprop
   training / apply continue to function (sanity-check by
   navigating back to **Browse > Platform > Predictive models**
   and observing that the saved `test_chemprop` entity is still
   present).
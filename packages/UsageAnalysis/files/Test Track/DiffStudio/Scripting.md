### Verify Scripting and Model Interaction in Diff Studio

1. Open **Diff Studio**. Open JavaScript Script:
* Click on the JS button on the ribbon. A view with a JavaScript script should open.
2. Run the Script:
* Execute the JavaScript script.
* Move the slider for the Final at input.
* Observe Changes: Verify that the table and line chart are modified in real-time as you move the slider.
3. Save the Script:
* Add "//tags: model" to JS body. Save the script after confirming that it runs correctly.
4. Access Model in Model Catalog:
* Go to Apps > Run Model catalog > Open the saved model from the catalog.
5. Interact with Model:
* Move the slider for the Final at input in the Model catalog.

Expected Results:

Step 2: The table and line chart should update automatically without any delays or errors as the slider is moved.

Step 3: The script should save successfully.

Step 5: The saved model should function correctly, and moving the slider should again modify the table and line chart in real-time.


---
{
  "order": 2
}
### Similatity search

1. Open the smiles.csv Dataset:
- Press the 'star' icon in TestTrack (with tooltip "Open test data"). The smiles.csv dataset opens.
2. Initiate Similarity Search:
- Go to Top menu > Chem > Search > Similatity search. 
- The search results should display containing molecules.
3. Access Similarity Search Properties:
- Click the gear icon on the Chem Similarity Search panel to open the properties settings.
4. Test Property Modifications (no errors should appear):
- Fingerprint: Change and test different fingerprints.
- Limit: Increase or decrease the limit of results.
- Distance Metric: Change and test different distance metrics.
- Size: Test with all available sizes (small, normal, large).
- Molecule Properties: Add a few properties from the list; these should be added to the "Most similar structures" information.
- Cutoff: Set the cutoff to 1; only one molecule should remain when the cutoff equals 1.

Verify that changing the properties works without any errors or crashes.

---
{
  "order": 1,
  "datasets": ["System:DemoFiles/chem/smiles.csv"]  
}
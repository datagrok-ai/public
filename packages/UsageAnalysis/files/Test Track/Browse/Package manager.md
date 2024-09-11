### Package Manager - Version Management and Uninstall Functionality

1. Navigate to Plugins:
* Go to Browse > Platform > Plugins.
2. Check Version Change:
* Versions should be switchable.
* Verify that when "Latest version" is selected, the system uses the latest stable version of the package, excluding versions marked as *.X.

After implementing features from ticket [GROK-16545](https://reddata.atlassian.net/browse/GROK-16545) check the next behavior:

3. Confirm that the "Delete" option is no longer present for any package. 

4. Select a package and use the "Uninstall" option. Verify that the package remains in the list but is shown as 'greyed out,' indicating it has been disabled but not deleted.

5. Install/Uninstall Specific Versions:
* Use the interface to install a specific version of the package (e.g., a *.X version) from the repository.
* Confirm the installation by checking that the correct version is marked as active.
* Uninstall the specific version using the context pane.
* Ensure that the uninstalled version behaves as expected (i.e., it is greyed out, not deleted).

6. Check Context Pane Bug Fix: 
* Install two or more versions of a package consecutively. Verify that only one version is shown as “green” (active) in the context pane at any time.

7. Install/Uninstall from Context Pane and Files View:
* Use the context pane to install a new package / reinstall a previously uninstalled package / svitch a version of a package. Verify that the action is successful.

---
{
  "order": 5
}

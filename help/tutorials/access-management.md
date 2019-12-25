<!-- TITLE: Use Cases: Access Management -->
<!-- SUBTITLE: -->

# Use Cases: Access Management

Owner: Alexander

Features: group-based privileges, nested groups, group admins, role management system, 
membership request system, audit, membership editor

TODO: write


## Creating and sharing a project
1. Open demog.csv
2. Add and customize a scatter plot
3. Click ‘Share’
4. Add recipients (test@datagrok.com) and click ‘OK’

Expected:
* Recipient should receive an email with the link to the project
* Recipient should receive a notification with the link to the project
* The view (along with the customized scatter plot) should be the same
* There should be an audit record 

## Creating and sharing a script
1. Open demog.csv
2. Tools | Scripting | New Script
3. Edit script (TODO)
4. Run. Expected: 
   * Executes normally and returns a result
   * Script appears in the script browser
   * Audit record appears that links script execution with the dataset (available under both dataset’s and script’s Activity pane in properties)
5. Click ‘Share’, enter recipients, OK
   * Expected: Recipient should receieve email and/or notification


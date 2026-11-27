@guide @tutorial
Feature: How-to · Files
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). A folder of Browse > Files > My files shared with a colleague:
  Share folder... from the folder's context menu, the colleague picked (the library's sharing user, bddsecond), the
  access levels shown and View and use kept, Send notifications off (nothing is sent), OK. The share makes a
  connection for the folder, "My files/Assays", which shows under Files after a refresh; Share folder... again lists
  the colleague with View and use. Files are shared by folder, never one by one (help: access/files, "File sharing
  and access control"). What the colleague then sees is left out: public#4268 (the second account gets "Credentials
  required" on opening the shared folder, sharing/sharing-folder.feature's known failure).

  The folder "Assays" holds one file (files/solubility-assay.csv of this folder), written through the API. The
  connection the share makes (MyFilesAssays on the server) is removed now and at feature end, its grant with it; the
  file and the folder are deleted when the feature ends (and a leftover folder when it starts).

  Scenario: Share a folder
    Share a folder of your files with a colleague, and choose what they may do with it.

    Given user is logged in
    And simple mode is off
    And the user's own connection "MyFilesAssays" is removed now and at feature end
    And no file "Assays" is in the user's files, now or when the feature ends
    And the "features/howto/14-files/files/solubility-assay.csv" file of the project is in the user's files as "Assays/solubility-assay.csv"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When user clicks on browse tab
    And user expands Files tree node inside browse tree
    And user clicks on Files---My-files tree node inside browse tree
    And user expands Files---My-files tree node inside browse tree
    Then Files---My-files---Assays tree node inside browse tree should be visible
    When the video says "Files are shared by folder: right-click the folder"
    And user picks "Share folder..." from the context menu of Files---My-files---Assays tree node inside browse tree
    Then "Share folder" dialog should be visible
    When the video says "Pick a colleague or a group"
    And user picks the sharing user from the suggestions of "User, group, or email" input in "Share folder" dialog
    And the video says "View and use: they can read the files, not change them"
    And user clicks on share access selector
    Then "View and use" tree node inside privilege tree should be checked
    When user presses Escape
    Then share access selector should contain text "View and use"
    # caption: Notifications are optional; off here
    When user unchecks "Send notifications" input in "Share folder" dialog
    And user clicks on OK button in "Share folder" dialog
    Then the "Share folder" dialog should close
    And an info balloon containing "Shared" should have been shown
    And 1 connection named "MyFilesAssays" should be on the server
    When user refreshes the browse tree
    Then "My files/Assays" tree node inside browse tree should be visible
    When the video says "The shared folder now shows under Files too"
    And user hovers over "My files/Assays" tree node inside browse tree
    And user expands Files---My-files tree node inside browse tree
    And user picks "Share folder..." from the context menu of Files---My-files---Assays tree node inside browse tree
    Then "Share folder" dialog should be visible
    And the access level of the sharing user in "Share folder" dialog should be "View and use"
    When the video says "The folder's sharing lists who has access"
    And user clicks on CANCEL button in "Share folder" dialog
    Then the "Share folder" dialog should close
    And no errors should have been logged

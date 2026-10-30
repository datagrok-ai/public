@journey @browse
Feature: Files: a folder and a file created, renamed and deleted from the Browse tree
  The File Manager "allows you to manage connections, browse and preview file content, and perform standard
  actions such as opening, downloading, deleting, and renaming. Right-click on an object to see available
  context actions" (`help/access/files/files.md`, "File Manager"; "create, delete, rename" in its
  introduction). A folder's menu offers Create folder... and Create file..., each asking for a name; a folder
  and a file offer Rename..., which shows the current name, and Delete..., which asks "Are you sure?". Each
  is claimed by the node the tree shows after OK, and the old name gone.

  Left out on purpose, one line each: deleting a top-level folder of My files (and its connection) is
  `sharing/sharing-folder.feature`'s; Open, Download and Download Zip are `browse-files.feature`'s; the
  content of a file made by Create file... and what the server stores are API tests (ApiTests, Dapi: files);
  creating, renaming and deleting a file share (a connection) are the connection features'.

  Everything happens inside a folder of the account's My files made for the run, "bdd-files-{time}", which
  holds one file written through the API; the folder and all it holds are deleted when the feature ends (and
  an older run's folder of the same family when it starts). Every scenario deletes what it made.

  Background:
    Given user is logged in
    And the "fixtures/cars-small.csv" file of the project is in the user's files as "bdd-files-{time}/seed.csv"
    And the browse panel is open
    When user refreshes the browse tree
    Given Files tree node inside browse tree is expanded
    And Files---My-files tree node inside browse tree is expanded
    And Files---My-files---bdd-files-{time} tree node inside browse tree is expanded

  Scenario: Create folder... asks for a name and the tree shows the new folder
    Then Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree should be visible
    When user picks "Create folder..." from the context menu of Files---My-files---bdd-files-{time} tree node inside browse tree
    Then "Create folder" dialog should be visible
    When user enters "sub-a" into Name input in "Create folder" dialog
    And user clicks on OK button in "Create folder" dialog
    Then the "Create folder" dialog should close
    And Files---My-files---bdd-files-{time}---sub-a tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: Rename... shows the folder's name and the tree shows the new one
    When user picks "Rename..." from the context menu of Files---My-files---bdd-files-{time}---sub-a tree node inside browse tree
    Then "Folder name" input in Rename dialog should have value "sub-a"
    When user enters "sub-b" into "Folder name" input in Rename dialog
    And user clicks on OK button in Rename dialog
    Then the "Rename" dialog should close
    And Files---My-files---bdd-files-{time}---sub-b tree node inside browse tree should be visible
    And Files---My-files---bdd-files-{time}---sub-a tree node inside browse tree should be absent
    And no errors should have been logged

  Scenario: Delete... asks for confirmation and the folder leaves the tree
    When user picks "Delete..." from the context menu of Files---My-files---bdd-files-{time}---sub-b tree node inside browse tree
    Then "Are you sure?" dialog should contain text "bdd-files-{time}/sub-b"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And Files---My-files---bdd-files-{time}---sub-b tree node inside browse tree should be absent
    And Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: Create file... asks for a name and the tree shows the new file
    When user picks "Create file..." from the context menu of Files---My-files---bdd-files-{time} tree node inside browse tree
    Then "Create file" dialog should be visible
    When user enters "notes-a.txt" into Name input in "Create file" dialog
    And user clicks on OK button in "Create file" dialog
    Then the "Create file" dialog should close
    And Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: Rename... shows the file's name and the tree shows the new one
    When user picks "Rename..." from the context menu of Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree
    Then "File name" input in Rename dialog should have value "notes-a.txt"
    When user enters "notes-b.txt" into "File name" input in Rename dialog
    And user clicks on OK button in Rename dialog
    Then the "Rename" dialog should close
    And Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree should be visible
    And Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree should be absent
    And no errors should have been logged

  Scenario: Delete... asks for confirmation and the file leaves the tree
    When user picks "Delete..." from the context menu of Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree
    Then "Are you sure?" dialog should contain text "notes-b.txt"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree should be absent
    And Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree should be visible
    And no errors should have been logged

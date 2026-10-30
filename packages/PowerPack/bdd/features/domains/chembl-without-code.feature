@guide @full-stand
Feature: Connect ChEMBL and get an application over it without code
  A guide for the pharma demo: "how fast can our own database become something our people use?". The database is a
  customer's ChEMBL 37 (71 tables, 2.9M compounds, 24.5M activities) in their own Postgres, reached with a login that
  may write. Create Domain Schema reads its tables, keys and foreign keys, recognizes the structures in its text
  columns, and registers the tables as domain tables that query ChEMBL live: nothing is copied. A synonym added, edited
  and deleted in Datagrok is read back from ChEMBL itself after each save. The connection, the schema and the synonym
  are removed at feature end.

  Scenario: Connect, name three tables, browse ChEMBL and change it from Datagrok
    Given user is logged in
    And simple mode is off
    And the context panel is closed
    And the user's own connection "ChEMBL37" is removed now and at feature end
    And the domain schema "chembl37" is removed now and at feature end
    And the "enableDomainDatabases" shell setting is put back at feature end
    And the "enableDomainDatabases" shell setting is on
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    # caption: The customer's ChEMBL 37: a Postgres database Datagrok has never seen
    When user picks "New connection..." from the context menu of Databases---Postgres tree node inside browse tree
    Then "Add new connection" dialog should be visible
    When user enters "ChEMBL37" into Name input in "Add new connection" dialog
    And user enters "host.docker.internal" into Server input in "Add new connection" dialog
    And user enters "5441" into Port input in "Add new connection" dialog
    And user enters "chembl" into Db input in "Add new connection" dialog
    And user enters "chembl_app" into Login input in "Add new connection" dialog
    And user enters "chembl_app" into Password input in "Add new connection" dialog
    And user clicks on OK button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And the synonym 9000001 is removed from the "ChEMBL37" database now and at feature end
    # caption: Right-click its schema › Create Domain Schema…: every table, key and link is read from the database
    When user expands Databases---Postgres tree node inside browse tree
    And user expands Databases---Postgres---ChEMBL37 tree node inside browse tree
    And user expands Databases---Postgres---ChEMBL37---Schemas tree node inside browse tree
    And user picks "Create Domain Schema..." from the context menu of Databases---Postgres---ChEMBL37---Schemas---public tree node inside browse tree
    Then "Create domain schema" dialog should be visible
    And "Create domain schema" dialog should contain text "71 of 71 bindable tables"
    # caption: Name it, let people write, and give three tables the names people use
    When user enters "ChEMBL" into "Friendly name" input in "Create domain schema" dialog
    And user enters "chembl37" into Identifier input in "Create domain schema" dialog
    And user checks Writable input in "Create domain schema" dialog
    And user clicks on the "molecule_dictionary" table in schema design
    And user enters "Compounds" into "Friendly name" input in "Create domain schema" dialog
    And user clicks on the "target_dictionary" table in schema design
    And user enters "Targets" into "Friendly name" input in "Create domain schema" dialog
    And user clicks on "Next" button in "Create domain schema" dialog
    And user clicks on "Validate" button in "Create domain schema" dialog
    Then "Create domain schema" dialog should contain text "No findings"
    When user clicks on "Create" button in "Create domain schema" dialog
    Then the "Create domain schema" dialog should close
    And the "chembl37.molecule_dictionary" domain table should hold 2921148 rows on the server
    # caption: The schema's page: its 71 tables, registered under the names given them — nothing copied
    And domain tables should contain text "Compounds"
    # caption: 2.9 million compounds, searched live in ChEMBL
    When user double-clicks on the "Compounds" card of domain tables
    Then domain status count should contain text "2,921,148"
    When user enters "IBRUTINIB" into domain search box
    And user opens the "IBRUTINIB" row of domain rows
    # caption: Its structure comes from compound_structures, found through ChEMBL's own keys
    Then domain header should contain text "IBRUTINIB"
    And domain header should show a structure
    And domain children should contain text "Activities (1685)"
    When user clicks on "Molecule synonyms (20)" tab in domain children
    Then domain children should show their rows
    # caption: Every link is a page: its mechanism of action, and the target it acts on
    When user clicks on "Drug mechanism (1)" tab in domain children
    Then domain children should show their rows
    When user clicks on the "cell 1 of mec_id" area of grid in domain children
    Then domain breadcrumbs should contain text "Drug mechanism"
    When user follows the "Tyrosine-protein kinase BTK" link in domain form
    Then domain breadcrumbs should contain text "Tyrosine-protein kinase BTK"
    And domain children should contain text "Assays (1662)"
    # caption: A curator adds a trade name: a new row in ChEMBL's molecule_synonyms
    When user opens the address "/domains/chembl37/molecule_synonyms?entity=new"
    And user picks "IBRUTINIB" in the "Molregno" link of domain form
    And user enters "TRADE_NAME" into "Syn type" input in domain form
    And user enters "9000001" into "Molsyn id" input in domain form
    And user enters "IMBRUVICA (DEMO)" into Synonyms input in domain form
    And user clicks on "Save" button in domain ribbon
    Then the synonym 9000001 in the "ChEMBL37" database should read "IMBRUVICA (DEMO)"
    # caption: …and it is in the customer's database
    When user reads "select molsyn_id, syn_type, synonyms from molecule_synonyms where molsyn_id = 9000001" from the "ChEMBL37" database as "In ChEMBL"
    Then "In ChEMBL" dialog should be visible
    When user clicks on OK button in "In ChEMBL" dialog
    # caption: An edit, saved: ChEMBL holds the new text
    And user enters "IMBRUVICA 560 MG (DEMO)" into Synonyms input in domain form
    And user clicks on "Save" button in domain ribbon
    Then the synonym 9000001 in the "ChEMBL37" database should read "IMBRUVICA 560 MG (DEMO)"
    When user reads "select molsyn_id, syn_type, synonyms from molecule_synonyms where molsyn_id = 9000001" from the "ChEMBL37" database as "In ChEMBL after the edit"
    Then "In ChEMBL after the edit" dialog should be visible
    When user clicks on OK button in "In ChEMBL after the edit" dialog
    # caption: A delete, saved: the row is gone from ChEMBL
    And user picks "Delete" from domain actions
    And user clicks on "Save" button in domain ribbon
    Then the synonym 9000001 should be gone from the "ChEMBL37" database
    When user reads "select molsyn_id, syn_type, synonyms from molecule_synonyms where molsyn_id = 9000001" from the "ChEMBL37" database as "In ChEMBL after the delete"
    Then "In ChEMBL after the delete" dialog should be visible
    When user clicks on OK button in "In ChEMBL after the delete" dialog

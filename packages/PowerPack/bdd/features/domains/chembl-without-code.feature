@guide @full-stand
Feature: Connect ChEMBL and get an application over it without code
  A guide for the pharma demo: "how fast can our own database become something our scientists use?". The database is
  the real ChEMBL 37 on db.datagrok.ai (71 tables, 2.4M compounds, 24.5M activities), read through a read-only login.
  A Postgres connection reaches it; Create Domain Schema reads its tables, primary and foreign keys, recognizes the
  structures in its text columns, and registers the tables as domain tables that query ChEMBL live. Nothing is copied.
  The connection and the schema are removed at feature end; ChEMBL itself is never changed.

  Scenario: Connect, bind every table, and walk ChEMBL along its own links
    Given user is logged in
    And simple mode is off
    And the user's own connection "ChEMBL37" is removed now and at feature end
    And the domain schema "chembl37" is removed now and at feature end
    And the "enableDomainDatabases" shell setting is put back at feature end
    And the "enableDomainDatabases" shell setting is on
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    # caption: A Postgres connection to ChEMBL, with a read-only login
    When user picks "New connection..." from the context menu of Databases---Postgres tree node inside browse tree
    Then "Add new connection" dialog should be visible
    When user enters "ChEMBL37" into Name input in "Add new connection" dialog
    And user enters "db.datagrok.ai" into Server input in "Add new connection" dialog
    And user enters "54325" into Port input in "Add new connection" dialog
    And user enters "chembl" into Db input in "Add new connection" dialog
    And user enters "datagrok" into Login input in "Add new connection" dialog
    And user enters "datagrok" into Password input in "Add new connection" dialog
    And user clicks on OK button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And 1 connection named "ChEMBL37" should be on the server
    # caption: Right-click its schema › Create Domain Schema…
    When user expands Databases---Postgres tree node inside browse tree
    And user expands Databases---Postgres---ChEMBL37 tree node inside browse tree
    And user expands Databases---Postgres---ChEMBL37---Schemas tree node inside browse tree
    And user picks "Create Domain Schema..." from the context menu of Databases---Postgres---ChEMBL37---Schemas---public tree node inside browse tree
    Then "Create domain schema" dialog should be visible
    # caption: All 71 tables are bindable: keys and links come from the database, structures are recognized
    And "Create domain schema" dialog should contain text "71 of 71 bindable tables"
    # caption: Name the schema, and give the tables people will open a readable name and a caption column
    When user enters "ChEMBL" into "Friendly name" input in "Create domain schema" dialog
    And user enters "chembl37" into Identifier input in "Create domain schema" dialog
    And user leaves the "chembl_id_lookup" table out of schema design
    And user clicks on the "target_dictionary" table in schema design
    And user enters "Targets" into "Friendly name" input in "Create domain schema" dialog
    And user selects "pref_name" in "Name column" input in "Create domain schema" dialog
    And user clicks on the "molecule_dictionary" table in schema design
    And user enters "Compounds" into "Friendly name" input in "Create domain schema" dialog
    And user selects "chembl_id" in "Name column" input in "Create domain schema" dialog
    And user selects "pref_name" in "Searchable" input in "Create domain schema" dialog
    And user clicks on the "assays" table in schema design
    And user enters "Assays" into "Friendly name" input in "Create domain schema" dialog
    And user selects "chembl_id" in "Name column" input in "Create domain schema" dialog
    And user clicks on the "docs" table in schema design
    And user enters "Documents" into "Friendly name" input in "Create domain schema" dialog
    And user selects "chembl_id" in "Name column" input in "Create domain schema" dialog
    And user clicks on "Next" button in "Create domain schema" dialog
    # caption: Validate checks every table against the live database before anything is registered
    And user clicks on "Validate" button in "Create domain schema" dialog
    Then "Create domain schema" dialog should contain text "No findings"
    When user clicks on "Create" button in "Create domain schema" dialog
    Then the "Create domain schema" dialog should close
    And the "chembl37.target_dictionary" domain table should hold 18552 rows on the server
    # caption: Every ChEMBL target, read live, with the condition builder beside them
    When user opens the address "/domains/chembl37/target_dictionary"
    Then domain status count should contain text "18,552"
    When user enters "Tyrosine-protein kinase BTK" into domain search box
    And user opens the "Tyrosine-protein kinase BTK" row of domain rows
    # caption: BTK, and the assays measured on it, found through ChEMBL's own foreign keys
    Then domain breadcrumbs should contain text "Tyrosine-protein kinase BTK"
    And domain children should contain text "Assays (1662)"
    # caption: 2.9 million compounds, searched live by name
    When user opens the address "/domains/chembl37/molecule_dictionary"
    And user enters "IBRUTINIB" into domain search box
    And user opens the "CHEMBL1873475" row of domain rows
    Then domain breadcrumbs should contain text "CHEMBL1873475"
    # caption: Its structure lives in another table; Create Domain Schema recognized it as a molecule
    And the "canonical_smiles" column of the "chembl37.compound_structures" domain table should be recognized as "Molecule"
    When user clicks on "Compound structures (1)" tab in domain children
    # caption: Its 1,685 activities, and one of them opened
    And user clicks on "Activities (1685)" tab in domain children
    And user clicks on the "cell 1 of activity_id" area of grid in domain children
    Then domain breadcrumbs should contain text "6975312"
    # caption: Every reference is a link: the assay it was measured in
    When user follows the "CHEMBL1794553" link in domain form
    Then domain breadcrumbs should contain text "CHEMBL1794553"
    And domain form should contain text "qHTS for Inhibitors of TGF-b"
    # caption: …and back to the compound, one click away
    When user goes back
    And user follows the "CHEMBL1201862" link in domain form
    Then domain breadcrumbs should contain text "CHEMBL1201862"
    And domain form should contain text "PubChem BioAssay data set"
    # caption: The document behind it: 7.4 million activities, counted live in ChEMBL
    And domain children should contain text "Activities (7434992)"
    And domain children should show their rows

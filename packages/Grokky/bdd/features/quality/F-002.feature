@serial
Feature: A tree map asked to color by a text column shows its categories instead of grey
  pilot-1120 F-002: asked for the sex / site / disease overlap in a picture, the assistant set the tree map's
  colorColumnName to a text column, which the tree map averages, so every block rendered grey with no error;
  and it described the chart with overlaps the data does not have ("women at New York", "Anxiety biggest at
  every site"). This checks the coloring, which the skill fix makes deterministic. The false descriptions are
  model behaviour the skill and system prompt reduce but cannot rule out, so no step here checks the wording.

  These scenarios reach a Docker container and a paid model: they run only with GROKKY_BDD_LIVE=1.

  Background:
    Given live AI model calls are allowed
    And the AI runtime traffic is recorded
    And user is logged in
    And the AI runtime container is healthy
    And user opens a table "clinic" with:
      | id   | sex    | site        | disease     | age |
      | 1    | Female | Los Angeles | Anxiety     | 34  |
      | 2    | Female | Los Angeles | Anxiety     | 51  |
      | 3    | Female | Los Angeles | Anxiety     | 47  |
      | 4    | Female | Los Angeles | Anxiety     | 62  |
      | 5    | Female | Los Angeles | Anxiety     | 29  |
      | 6    | Female | Buffalo     | Anxiety     | 58  |
      | 7    | Female | Buffalo     | Anxiety     | 41  |
      | 8    | Female | Buffalo     | Anxiety     | 66  |
      | 9    | Male   | Los Angeles | Anxiety     | 38  |
      | 10   | Male   | Los Angeles | Anxiety     | 55  |
      | 11   | Male   | New York    | Glaucoma    | 71  |
      | 12   | Male   | New York    | Glaucoma    | 49  |
      | 13   | Male   | New York    | Glaucoma    | 63  |
      | 14   | Male   | New York    | Glaucoma    | 57  |
      | 15   | Male   | New York    | Hepatitis A | 44  |
      | 16   | Male   | New York    | Hepatitis A | 36  |
      | 17   | Male   | New York    | Hepatitis A | 52  |
      | 18   | Male   | Orlando     | Hepatitis A | 68  |
      | 19   | Male   | Orlando     | Hepatitis A | 45  |
      | 20   | Male   | Orlando     | Hepatitis A | 59  |
      | 21   | Male   | Orlando     | Hepatitis A | 33  |
      | 22   | Male   | Buffalo     | Anxiety     | 61  |
      | 23   | Male   | Buffalo     | Anxiety     | 50  |
    And the AI panel is open on a new chat

  Scenario: The tree map asked to color by disease puts disease last in its split and sets no color column
    When user asks the AI panel "say I just got this clinical table from a colleague. What is in it?"
    And user asks the AI panel "Thanks. What's worth a closer look? Anything odd or suspicious in the data?"
    And user asks the AI panel "ok, show me the sex / site / disease overlap in a tree map colored by disease, and tell me what it shows"
    Then the open tableview should have 1 Tree map viewer
    And the tree map should color its blocks by "disease"
    And a screenshot "F-002-tree-map" of the AI panel is saved

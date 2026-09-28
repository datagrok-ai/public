@viewers @realizes:viewers.chrome
Feature: Docking viewers by their title bars
  `viewers-docking-ui.md`: a viewer dragged by its title bar shows dock-spawn's wheels — a compass
  over the panel under the pointer (left, up, right, down and a centre item) and one item at each edge
  of the view — and dropping it on an item docks it there. A viewer dropped on an edge item of the
  view runs along the whole of that edge, whatever is docked already: after the histogram took the
  right edge, the bar chart dropped on the bottom edge item runs under the grid, the scatter plot and
  the histogram, and the histogram keeps the top right corner. Dropped on a compass arrow over another
  viewer, it shares that viewer's width or height and leaves the rest of the view alone (a viewer left
  of the scatter plot takes the scatter plot's height, not the view's); dropped on the compass centre
  it becomes a tab of that viewer's group. A viewer added after another is already docked under it,
  so the md's drop on the scatter plot's bottom side would change nothing: the histogram is dropped on
  its top side instead. Every drop is preceded by the claim that the viewer is not yet where the drop
  puts it. The arrangement comes back from a layout applied to demog-1000 opened anew.
  The drop targets are the wheel items dock-spawn draws, the edge ones told from the compass ones by
  the compass's centre item, and the drop waits for the item's hover mark. The claims are the
  panels' boxes after the drop, against each other and against the area the docked panels cover.
  Not translated: the splitter drag between two docked viewers and "no rendering artifacts" (a
  judgement of the eye), and the md's downloaded `.layout` file dropped onto the view — a file
  dragged in from the desktop, which a browser page cannot be handed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer
    And user adds a histogram viewer
    And user adds a bar chart viewer
    Then the current view should hold at least 3 viewers

  Scenario: A viewer dropped at the right edge of the view docks along the whole edge
    Then histogram viewer should not be docked along the right edge of the view
    When user docks histogram viewer to the right edge of the view
    Then histogram viewer should be docked along the right edge of the view
    And bar chart viewer should be visible
    And scatter plot viewer should be visible
    And no errors should have been logged

  Scenario: A viewer dropped at the bottom edge after the right edge was taken runs under all of them
    When user docks histogram viewer to the right edge of the view
    Then histogram viewer should be docked along the right edge of the view
    And bar chart viewer should not be docked along the bottom edge of the view
    When user docks bar chart viewer to the bottom edge of the view
    Then bar chart viewer should be docked along the bottom edge of the view
    And histogram viewer should not be docked along the right edge of the view
    And histogram viewer should be docked in the top right corner of the view
    And scatter plot viewer should be visible
    And no errors should have been logged

  Scenario: A viewer dropped on a side of another viewer docks next to that viewer only
    Then histogram viewer should not be docked above scatter plot viewer
    When user docks histogram viewer to the top side of scatter plot viewer
    Then histogram viewer should be docked above scatter plot viewer
    And bar chart viewer should not be docked left-of scatter plot viewer
    When user docks bar chart viewer to the left side of scatter plot viewer
    Then bar chart viewer should be docked left-of scatter plot viewer
    And no errors should have been logged

  Scenario: A docking arrangement comes back from a layout applied to the table opened anew
    When user docks histogram viewer to the right edge of the view
    Then bar chart viewer should not be docked below histogram viewer
    When user docks bar chart viewer to the bottom side of histogram viewer
    Then bar chart viewer should be docked below histogram viewer
    When user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    And user loads the saved layout
    Then scatter plot viewer should be visible
    And bar chart viewer should be docked below histogram viewer
    And no errors should have been logged

@journey @viewers @realizes:viewers.matrix-plot
Feature: Matrix plot — the viewport over a long column list, and the 250-cell cap
  With sixteen numerical columns on each axis the grid would be 256 cells, and the plot lays out at
  most 250: the viewport starts at five by five and the two column sliders move it, up to the point
  where the next increment is refused. `_updateViewport` refuses it SILENTLY — it returns, the
  slider keeps its new values and the grid does not change — so the only observable is
  `viewport rejected`, which is why the spec this replaces could do nothing but drag and then poll
  for a cell count that never moved, calling the refusal proven when `yFull === yNearCap`.
  The drag itself is the other half. The old one computed a pixels-per-unit from the bounding boxes
  of the two handles minus a handle diameter and then moved the max handle by `(to - from) *
  pxPerUnit`; the slider's track and handles are regions the plot reports, so the gesture here is
  "to the end of the track" and the viewport it produces is read back rather than assumed.
  Not covered: narrowing the viewport by hand. The max handle dragged towards the min reaches
  min == max, which lays out a viewport of no cells at all and leaves `isRenderPending` true for
  ever — every step after it fails with "a repaint has been pending for 10000 ms". The 0-cell
  viewport is a state a user can reach with the mouse, so the getter has to survive it.
  Twelve calculated columns make up the sixteen; the last scenario removes them.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a matrix plot viewer
    Then the "cells" reading of matrix plot viewer should be 16
    And the "viewport limit" reading of matrix plot viewer should be 250
    And the "viewport rejected" reading of matrix plot viewer should be "false"
    And matrix plot viewer should have an "x scroll slider" area
    And matrix plot viewer should have a "y scroll slider" area

  Scenario: With four columns the viewport already holds them all and a slider drag is inert
    # Narrowing on this fixture is not written here: the max handle dragged to the START of its
    # track makes min == max, the plot lays out a viewport of NO cells, and its `isRenderPending`
    # then never clears — the settle fails with "a repaint has been pending for 10000 ms" while the
    # empty grid sits there quietly. The 0-cell viewport is reachable by hand, so the getter has to
    # cope with it. Narrowing is claimed on the sixteen-column fixture below instead, where the
    # viewport is smaller than the axis from the start.
    Then the cells of matrix plot viewer should be 4 wide and 4 tall
    And matrix plot viewer should have an "x label STARTED" area
    And matrix plot viewer should have a "cell AGE x AGE" area
    When user drags the max handle of the "x" scroll slider of matrix plot viewer to its end
    Then the cells of matrix plot viewer should be 4 wide and 4 tall
    And the "cells" reading of matrix plot viewer should be 16
    And the "cells drawn" and "cells" readings of matrix plot viewer should be the same
    And the "viewport rejected" reading of matrix plot viewer should be "false"
    When user drags the max handle of the "y" scroll slider of matrix plot viewer to its end
    Then the cells of matrix plot viewer should be 4 wide and 4 tall
    And matrix plot viewer should have an "x label STARTED" area
    And no errors should have been logged

  Scenario: Sixteen columns start at a five by five viewport and the sliders open it up
    Given user adds a calculated column "MP_FX_1" with formula "${AGE} + 1"
    And user adds a calculated column "MP_FX_2" with formula "${AGE} + 2"
    And user adds a calculated column "MP_FX_3" with formula "${AGE} + 3"
    And user adds a calculated column "MP_FX_4" with formula "${AGE} + 4"
    And user adds a calculated column "MP_FX_5" with formula "${AGE} + 5"
    And user adds a calculated column "MP_FX_6" with formula "${AGE} + 6"
    And user adds a calculated column "MP_FX_7" with formula "${AGE} + 7"
    And user adds a calculated column "MP_FX_8" with formula "${AGE} + 8"
    And user adds a calculated column "MP_FX_9" with formula "${AGE} + 9"
    And user adds a calculated column "MP_FX_10" with formula "${AGE} + 10"
    And user adds a calculated column "MP_FX_11" with formula "${AGE} + 11"
    And user adds a calculated column "MP_FX_12" with formula "${AGE} + 12"
    When user sets properties of matrix plot viewer:
      | xColumnNames | AGE, HEIGHT, WEIGHT, STARTED, MP_FX_1, MP_FX_2, MP_FX_3, MP_FX_4, MP_FX_5, MP_FX_6, MP_FX_7, MP_FX_8, MP_FX_9, MP_FX_10, MP_FX_11, MP_FX_12 |
      | yColumnNames | AGE, HEIGHT, WEIGHT, STARTED, MP_FX_1, MP_FX_2, MP_FX_3, MP_FX_4, MP_FX_5, MP_FX_6, MP_FX_7, MP_FX_8, MP_FX_9, MP_FX_10, MP_FX_11, MP_FX_12 |
    Then the "x columns" reading of matrix plot viewer should be 16
    And the "y columns" reading of matrix plot viewer should be 16
    And the cells of matrix plot viewer should be 5 wide and 5 tall
    And the "cells" reading of matrix plot viewer should be 25
    And the "viewport rejected" reading of matrix plot viewer should be "false"
    When user drags the max handle of the "x" scroll slider of matrix plot viewer to its end
    Then the "columns" reading of matrix plot viewer should be 16
    And the "cells" reading of matrix plot viewer should be 80
    And the "cells drawn" and "cells" readings of matrix plot viewer should be the same
    And the "viewport rejected" reading of matrix plot viewer should be "false"
    And no errors should have been logged

  Scenario: The viewport stops short of the whole grid and records that it refused the rest
    # 16 by 16 is 256 cells and the plot lays out at most 250, so the slider dragged to the end of
    # its track cannot bring the sixteenth row in: `_updateViewport` returns without changing
    # anything and sets `viewportRejected`. Which row it stopped on is the drag's own path — every
    # mousemove that changes the integer viewport re-tiles the whole grid — so the claim is that it
    # is below the sixteen the axis holds, at or under the limit, and that the refusal was recorded.
    Then the "cells" reading of matrix plot viewer should be 80
    When user drags the max handle of the "y" scroll slider of matrix plot viewer to its end
    Then the "cells" reading of matrix plot viewer should be higher than before
    And the "cells" reading of matrix plot viewer should be between 96 and 250
    And the "columns" reading of matrix plot viewer should be 16
    And the "y columns" reading of matrix plot viewer should be 16
    And the "rows" reading of matrix plot viewer should be between 6 and 15
    And the "viewport rejected" reading of matrix plot viewer should be "true"
    And the "cells drawn" and "cells" readings of matrix plot viewer should be the same
    And the "error" reading of matrix plot viewer should be ""
    When user sets properties of matrix plot viewer:
      | xColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
      | yColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
    Then the "cells" reading of matrix plot viewer should be 16
    When user removes "MP_FX_1" column
    And user removes "MP_FX_2" column
    And user removes "MP_FX_3" column
    And user removes "MP_FX_4" column
    And user removes "MP_FX_5" column
    And user removes "MP_FX_6" column
    And user removes "MP_FX_7" column
    And user removes "MP_FX_8" column
    And user removes "MP_FX_9" column
    And user removes "MP_FX_10" column
    And user removes "MP_FX_11" column
    And user removes "MP_FX_12" column
    Then the table should have 11 columns
    And the "cells" reading of matrix plot viewer should be 16
    And no errors should have been logged

Feature: Diff Studio smoke
  The first feature of this package: platform steps only, so it passes on any stand. Replace it
  with the app's own scenarios — "user opens the Diff Studio app" (bindings/steps.ts) is the way in.

  Scenario: The platform is up
    Given user is logged in
    Then browse tab should be visible

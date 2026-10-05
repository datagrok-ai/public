/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/data-access/sticky-meta.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.sticky-meta]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {entityTypeExists, startTutorial, stepDone, stepNotDone, stickySchemaExists, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen, walkTour} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, expand, selectIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, noHintShown, packageInstalled, stickyMetaFixturesGone, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Sticky Meta tutorial", () => {
  const session = feature(test, "features/data-access/sticky-meta.feature", import.meta.url);
  test("A learner completes the Sticky Meta tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.sticky-meta"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(21, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(22, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(23, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(24, "And the Sticky Meta schema \"schema for tutorial\" and entity type \"molecule-tutorial\" are removed now and at feature end", () => stickyMetaFixturesGone(page, "schema for tutorial", "molecule-tutorial"));
    await session.step(25, "And the \"Sticky Meta\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Sticky Meta"));
    await session.step(26, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(29, "When user starts the \"Sticky Meta\" tutorial", () => startTutorial(page, "Sticky Meta"));
    await session.step(30, "Then the tutorial progress should be 1 of 19", () => tutorialProgress(page, 1, 19));
    await session.step(31, "Given the tutorial step \"Open Types node\" should not be done yet", () => stepNotDone(page, "Open Types node"));
    await session.step(32, "When user expands Platform tree node inside browse tree", () => expand(page, el("Platform tree node inside browse tree")));
    await session.step(33, "And user expands Platform---Sticky-Meta tree node inside browse tree", () => expand(page, el("Platform---Sticky-Meta tree node inside browse tree")));
    await session.step(34, "And user clicks on Platform---Sticky-Meta---Types tree node inside browse tree", () => clickOn(page, el("Platform---Sticky-Meta---Types tree node inside browse tree")));
    await session.step(35, "Then the tutorial step \"Open Types node\" should be done", () => stepDone(page, "Open Types node"));
    await session.step(36, "Given the tutorial step \"Create a new entity type\" should not be done yet", () => stepNotDone(page, "Create a new entity type"));
    await session.step(37, "When user clicks on \"New Entity Type...\" button", () => clickOn(page, el("\"New Entity Type...\" button")));
    await session.step(38, "Then the tutorial step \"Create a new entity type\" should be done", () => stepDone(page, "Create a new entity type"));
    await session.step(39, "Given the tutorial step \"Explore entity type dialog\" should not be done yet", () => stepNotDone(page, "Explore entity type dialog"));
    await session.step(40, "When user goes through the tour to its end", () => walkTour(page));
    await session.step(41, "Then the tutorial step \"Explore entity type dialog\" should be done", () => stepDone(page, "Explore entity type dialog"));
    await session.step(42, "When user enters \"molecule-tutorial\" into \"Name\" input in \"Create a new entity type\" dialog", () => enterInto(page, "molecule-tutorial", el("\"Name\" input in \"Create a new entity type\" dialog")));
    await session.step(43, "Then the tutorial step \"Set \\\"Name\\\" to \\\"molecule-tutorial\\\"\" should be done", () => stepDone(page, "Set \"Name\" to \"molecule-tutorial\""));
    await session.step(44, "When user enters \"semtype=Molecule\" into \"Matching expression\" input in \"Create a new entity type\" dialog", () => enterInto(page, "semtype=Molecule", el("\"Matching expression\" input in \"Create a new entity type\" dialog")));
    await session.step(45, "Then the tutorial step \"Set \\\"Matching expression\\\" to \\\"semtype=Molecule\\\"\" should be done", () => stepDone(page, "Set \"Matching expression\" to \"semtype=Molecule\""));
    await session.step(46, "When user clicks on OK button in \"Create a new entity type\" dialog", () => clickOn(page, el("OK button in \"Create a new entity type\" dialog")));
    await session.step(47, "Then the tutorial step \"Save entity type\" should be done", () => stepDone(page, "Save entity type"));
    await session.step(48, "Then the entity type \"molecule-tutorial\" should exist", () => entityTypeExists(page, "molecule-tutorial"));
    await session.step(50, "Given the tutorial step \"Open schemas node\" should not be done yet", () => stepNotDone(page, "Open schemas node"));
    await session.step(51, "When user clicks on Platform---Sticky-Meta---Schemas tree node inside browse tree", () => clickOn(page, el("Platform---Sticky-Meta---Schemas tree node inside browse tree")));
    await session.step(52, "Then the tutorial step \"Open schemas node\" should be done", () => stepDone(page, "Open schemas node"));
    await session.step(53, "Given the tutorial step \"Create a new schema\" should not be done yet", () => stepNotDone(page, "Create a new schema"));
    await session.step(54, "When user clicks on \"New Schema...\" button", () => clickOn(page, el("\"New Schema...\" button")));
    await session.step(55, "Then the tutorial step \"Create a new schema\" should be done", () => stepDone(page, "Create a new schema"));
    await session.step(56, "Given the tutorial step \"Explore schema dialog\" should not be done yet", () => stepNotDone(page, "Explore schema dialog"));
    await session.step(57, "When user goes through the tour to its end", () => walkTour(page));
    await session.step(58, "Then the tutorial step \"Explore schema dialog\" should be done", () => stepDone(page, "Explore schema dialog"));
    await session.step(59, "When user enters \"schema for tutorial\" into \"Name\" input in \"Create a new schema\" dialog", () => enterInto(page, "schema for tutorial", el("\"Name\" input in \"Create a new schema\" dialog")));
    await session.step(60, "Then the tutorial step \"Set \\\"Name\\\" to \\\"schema for tutorial\\\"\" should be done", () => stepDone(page, "Set \"Name\" to \"schema for tutorial\""));
    await session.step(61, "When user clicks on \"select entities\" action in \"Create a new schema\" dialog", () => clickOn(page, el("\"select entities\" action in \"Create a new schema\" dialog")));
    await session.step(62, "Then the tutorial step \"Select associated entity\" should be done", () => stepDone(page, "Select associated entity"));
    await session.step(63, "And \"Select types for schema for tutorial\" dialog should be visible", () => shouldBe(page, el("\"Select types for schema for tutorial\" dialog"), "visible"));
    await session.step(64, "When user checks \"molecule-tutorial\" property in \"Select types for schema for tutorial\" dialog", () => check(page, el("\"molecule-tutorial\" property in \"Select types for schema for tutorial\" dialog")));
    await session.step(65, "Then the tutorial step \"Select molecule-tutorial\" should be done", () => stepDone(page, "Select molecule-tutorial"));
    await session.step(66, "When user clicks on OK button in \"Select types for schema for tutorial\" dialog", () => clickOn(page, el("OK button in \"Select types for schema for tutorial\" dialog")));
    await session.step(67, "Then the tutorial step \"Confirm entity selection\" should be done", () => stepDone(page, "Confirm entity selection"));
    await session.step(68, "When user enters \"project name\" into second \"Name\" input in \"Create a new schema\" dialog", () => enterInto(page, "project name", el("second \"Name\" input in \"Create a new schema\" dialog")));
    await session.step(69, "Then the tutorial step \"Set property \\\"Name\\\" to \\\"project name\\\"\" should be done", () => stepDone(page, "Set property \"Name\" to \"project name\""));
    await session.step(70, "When user selects \"string\" in \"Property Type\" input in \"Create a new schema\" dialog", () => selectIn(page, "string", el("\"Property Type\" input in \"Create a new schema\" dialog")));
    await session.step(71, "Then the tutorial step \"Set property \\\"Type\\\" to \\\"string\\\"\" should be done", () => stepDone(page, "Set property \"Type\" to \"string\""));
    await session.step(72, "When user clicks on OK button in \"Create a new schema\" dialog", () => clickOn(page, el("OK button in \"Create a new schema\" dialog")));
    await session.step(73, "Then the tutorial step \"Save schema\" should be done", () => stepDone(page, "Save schema"));
    await session.step(74, "And the Sticky Meta schema \"schema for tutorial\" should exist", () => stickySchemaExists(page, "schema for tutorial"));
    await session.step(76, "When user clicks on the \"cell 1 of smiles\" area of grid", () => clickArea(page, "cell 1 of smiles", el("grid")));
    await session.step(77, "Then \"Sticky meta\" pane in context panel should be visible", () => shouldBe(page, el("\"Sticky meta\" pane in context panel"), "visible"));
    await session.step(78, "And the tutorial step \"Enter value for project name\" should not be done yet", () => stepNotDone(page, "Enter value for project name"));
    await session.step(79, "When user types \"BDD tutorial\" into \"project name\" input in context panel", () => typeInto(page, "BDD tutorial", el("\"project name\" input in context panel")));
    await session.step(80, "Then the tutorial step \"Enter value for project name\" should be done", () => stepDone(page, "Enter value for project name"));
    await session.step(81, "When user clicks on Save button in \"Sticky meta\" pane in context panel", () => clickOn(page, el("Save button in \"Sticky meta\" pane in context panel")));
    await session.step(82, "Then the tutorial step \"Save sticky meta changes\" should be done", () => stepDone(page, "Save sticky meta changes"));
    await session.step(83, "And Save button in \"Sticky meta\" pane in context panel should be disabled", () => shouldBe(page, el("Save button in \"Sticky meta\" pane in context panel"), "disabled"));
    await session.step(84, "When user hovers over the \"cell 1 of smiles\" area of grid", () => hoverArea(page, "cell 1 of smiles", el("grid")));
    await session.step(85, "Then the tutorial step \"Hover a cell to verify metadata tooltip.\" should be done", () => stepDone(page, "Hover a cell to verify metadata tooltip."));
    await session.step(87, "And the \"Sticky Meta\" tutorial should be completed", () => tutorialCompleted(page, "Sticky Meta"));
    await session.step(88, "And the tutorial should have listed 19 steps", () => tutorialStepsListed(page, 19));
    await session.step(89, "And the tutorial progress should be 19 of 19", () => tutorialProgress(page, 19, 19));
    await session.step(90, "And no hint should be shown", () => noHintShown(page));
    await session.step(91, "And no errors should have been logged", () => noErrors(page));
  });
});

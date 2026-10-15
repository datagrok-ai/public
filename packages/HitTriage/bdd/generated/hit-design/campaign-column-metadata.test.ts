/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/hit-design/campaign-column-metadata.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addTag, campaignOnServer} from '../../bindings/hit-design.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, pressKeyIn, shouldOfferChoice, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openApp, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea, infoBalloonText, noBalloons, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A Hit Design campaign keeps its columns' tags and formulas when it is saved and reopened", () => {
  const session = feature(test, "features/hit-design/campaign-column-metadata.feature", import.meta.url);
  test("A Hit Design campaign keeps its columns' tags and formulas when it is saved and reopened", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the Hit Design campaign \"BDD Columns\" with a table of:", () => campaignOnServer(page, "BDD Columns", [["Molecule","Stage","V-iD","Score","Notes"],["CCO","Design","","10",""],["c1ccccc1","Design","","4",""]]), [["Molecule","Stage","V-iD","Score","Notes"],["CCO","Design","","10",""],["c1ccccc1","Design","","4",""]]);
    await session.step(18, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
    await session.step(19, "When user clicks on \"BDD Columns\" link", () => clickOn(page, el("\"BDD Columns\" link")));
    await session.step(20, "Then the \"BDD Columns\" view should be current", () => viewIsCurrent(page, "BDD Columns"));
    await session.step(21, "And the table should have 2 rows", () => rowCount(page, 2));
    await run.scenario("A tag set on an empty column is there when the campaign is reopened", async () => {
      await session.step(24, "When user picks \"Column Properties...\" from the context menu of the \"header Notes\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header Notes", el("grid")));
      await session.step(25, "And user adds the \".choices\" tag with value '[\"low\",\"high\"]' to \"Tags\" input in \"Notes\" dialog", () => addTag(page, ".choices", "[\"low\",\"high\"]", el("\"Tags\" input in \"Notes\" dialog")));
      await session.step(26, "And user clicks on OK button in \"Notes\" dialog", () => clickOn(page, el("OK button in \"Notes\" dialog")));
      await session.step(27, "And user double-clicks on the \"cell 1 of Notes\" area of grid", () => doubleClickArea(page, "cell 1 of Notes", el("grid")));
      await session.step(28, "Then cell editor should offer the choice \"low\"", () => shouldOfferChoice(page, el("cell editor"), "low"));
      await session.step(29, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(30, "And user clicks on \"Submit...\" button", () => clickOn(page, el("\"Submit...\" button")));
      await session.step(31, "And user clicks on Save button in \"Submit\" dialog", () => clickOn(page, el("Save button in \"Submit\" dialog")));
      await session.step(32, "Then an info balloon containing \"Campaign saved successfully\" should have been shown", () => infoBalloonText(page, "Campaign saved successfully"));
      await session.step(33, "When user clicks on \"Hit Design\" link", () => clickOn(page, el("\"Hit Design\" link")));
      await session.step(34, "And user clicks on \"BDD Columns\" link", () => clickOn(page, el("\"BDD Columns\" link")));
      await session.step(35, "Then the \"BDD Columns\" view should be current", () => viewIsCurrent(page, "BDD Columns"));
      await session.step(36, "And the table should have the columns \"Molecule, Stage, V-iD, Score, Notes\"", () => columnsExactly(page, "Molecule, Stage, V-iD, Score, Notes"));
      await session.step(37, "And \"Notes\" column should have tag \".choices\" equal to '[\"low\",\"high\"]'", () => columnTag(page, "Notes", ".choices", "[\"low\",\"high\"]"));
      await session.step(38, "When user double-clicks on the \"cell 1 of Notes\" area of grid", () => doubleClickArea(page, "cell 1 of Notes", el("grid")));
      await session.step(39, "Then cell editor should offer the choice \"low\"", () => shouldOfferChoice(page, el("cell editor"), "low"));
      await session.step(40, "And cell editor should offer the choice \"high\"", () => shouldOfferChoice(page, el("cell editor"), "high"));
      await session.step(41, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(42, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A calculated column keeps its formula and follows its source when the campaign is reopened", async () => {
      await session.step(45, "When user adds a calculated column \"Score x2\" with formula \"${Score} * 2\"", () => addCalculated(page, "Score x2", "${Score} * 2"));
      await session.step(46, "Then the value of \"Score x2\" column in row 1 should be \"20\"", () => valueInRow(page, "Score x2", 1, "20"));
      await session.step(47, "When user clicks on \"Submit...\" button", () => clickOn(page, el("\"Submit...\" button")));
      await session.step(48, "And user clicks on Save button in \"Submit\" dialog", () => clickOn(page, el("Save button in \"Submit\" dialog")));
      await session.step(49, "Then an info balloon containing \"Campaign saved successfully\" should have been shown", () => infoBalloonText(page, "Campaign saved successfully"));
      await session.step(50, "When user clicks on \"Hit Design\" link", () => clickOn(page, el("\"Hit Design\" link")));
      await session.step(51, "And user clicks on \"BDD Columns\" link", () => clickOn(page, el("\"BDD Columns\" link")));
      await session.step(52, "Then the \"BDD Columns\" view should be current", () => viewIsCurrent(page, "BDD Columns"));
      await session.step(53, "And the table should have the columns \"Molecule, Stage, V-iD, Score, Notes, Score x2\"", () => columnsExactly(page, "Molecule, Stage, V-iD, Score, Notes, Score x2"));
      await session.step(54, "And \"Score x2\" column should have tag \"formula\" equal to \"${Score} * 2\"", () => columnTag(page, "Score x2", "formula", "${Score} * 2"));
      await session.step(55, "When user double-clicks on the \"cell 1 of Score\" area of grid", () => doubleClickArea(page, "cell 1 of Score", el("grid")));
      await session.step(56, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(57, "And user types \"7\" into cell editor", () => typeInto(page, "7", el("cell editor")));
      await session.step(58, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(59, "Then the value of \"Score\" column in row 1 should be \"7\"", () => valueInRow(page, "Score", 1, "7"));
      await session.step(60, "And the value of \"Score x2\" column in row 1 should be \"14\"", () => valueInRow(page, "Score x2", 1, "14"));
      await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});

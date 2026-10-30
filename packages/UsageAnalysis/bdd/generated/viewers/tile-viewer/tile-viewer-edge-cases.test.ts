/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer, core#1217]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {addTileViewerWith, deleteValueField, designerValueFields, dragCardIntoLane} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, selectNoRows, selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pointerAway, readingIs, readingReads, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer edge cases: nothing selected under Row Source Selected, one row", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-edge-cases.feature", import.meta.url);
  test("Row Source Selected with nothing selected shows no card and no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
    await session.step(15, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(16, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(19, "Given user selects no rows", () => selectNoRows(page));
    await session.step(20, "When user sets \"Row Source\" property of tile viewer to \"Selected\"", () => setProperty(page, "Row Source", el("tile viewer"), "Selected"));
    await session.step(21, "Then tile viewer should show 0 rows", () => showsRows(page, el("tile viewer"), 0));
    await session.step(22, "And the \"tiles\" reading of tile viewer should be 0", () => readingIs(page, "tiles", el("tile viewer"), 0));
    await session.step(23, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
    await session.step(24, "Then tile viewer should show 15 rows", () => showsRows(page, el("tile viewer"), 15));
    await session.step(25, "And the \"lane of row 10\" reading of tile viewer should be \"Asian\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Asian"));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
  });
  test("A one-row table shows its one card", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
    await session.step(15, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(16, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(29, "Given user opens a table \"one-row\" with:", () => openTableOf(page, "one-row", [["name","group","x"],["a","g1","42"]]), [["name","group","x"],["a","g1","42"]]);
    await session.step(32, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","group"]]), [["Lanes Column Name","group"]]);
    await session.step(34, "Then tile viewer should show 1 rows", () => showsRows(page, el("tile viewer"), 1));
    await session.step(35, "And the \"tiles\" reading of tile viewer should be 1", () => readingIs(page, "tiles", el("tile viewer"), 1));
    await session.step(36, "And the \"lane of row 1\" reading of tile viewer should be \"g1\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "g1"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("Escape in the form designer, a view and not a dialog, leaves it open with the field deleted", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
    await session.step(15, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(16, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(40, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
    await session.step(41, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
    await session.step(42, "When user deletes the \"AGE\" value field in the form designer", () => deleteValueField(page, "AGE"));
    await session.step(43, "Then the form designer should show value fields \"CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
    await session.step(44, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(45, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
    await session.step(46, "And the form designer should show value fields \"CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT\"", () => designerValueFields(page, "CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
  test("A card dropped into the (empty) lane gets an empty value", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tile-viewer", "@known-failure:62", "@realizes:core#1217"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
    await session.step(15, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(16, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
    await knownFailure(async () => {
      await session.step(53, "Given user opens a table \"with-empty\" with:", () => openTableOf(page, "with-empty", [["site","arm"],["A","x"],["A",""],["B","y"]]), [["site","arm"],["A","x"],["A",""],["B","y"]]);
      await session.step(58, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","arm"]]), [["Lanes Column Name","arm"]]);
      await session.step(60, "Then the \"lane of row 1\" reading of tile viewer should be \"x\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "x"));
      await session.step(61, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
      await session.step(62, "And user drags the card of row 1 of tile viewer into lane \"(empty)\"", () => dragCardIntoLane(page, 1, el("tile viewer"), "(empty)"));
      await session.step(63, "Then the value of \"arm\" column in row 1 should be \"\"", () => valueInRow(page, "arm", 1, ""));
      await session.step(64, "And the \"lane of row 1\" reading of tile viewer should be \"(empty)\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "(empty)"));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    }, 62);
  });
  test("A card dragged between the lanes of a calculated column writes the lane's value", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
    await session.step(15, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(16, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(69, "When user adds a calculated column \"AGE_GROUP\" with formula \"if(${AGE} > 40, \\\"old\\\", \\\"young\\\")\"", () => addCalculated(page, "AGE_GROUP", "if(${AGE} > 40, \"old\", \"young\")"));
    await session.step(70, "And user sets \"lanesColumnName\" property of tile viewer to \"AGE_GROUP\"", () => setProperty(page, "lanesColumnName", el("tile viewer"), "AGE_GROUP"));
    await session.step(71, "Then the \"lane of row 1\" reading of tile viewer should be \"young\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "young"));
    await session.step(72, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
    await session.step(73, "And user drags the card of row 1 of tile viewer into lane \"old\"", () => dragCardIntoLane(page, 1, el("tile viewer"), "old"));
    await session.step(74, "Then the \"lane of row 1\" reading of tile viewer should be \"old\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "old"));
    await session.step(75, "And the value of \"AGE_GROUP\" column in row 1 should be \"old\"", () => valueInRow(page, "AGE_GROUP", 1, "old"));
    await session.step(76, "And the value of \"AGE\" column in row 1 should be \"26\"", () => valueInRow(page, "AGE", 1, "26"));
    await session.step(77, "And no errors should have been logged", () => noErrors(page));
  });
});

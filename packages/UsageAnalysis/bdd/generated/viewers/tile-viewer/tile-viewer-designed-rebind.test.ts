/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-designed-rebind.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer]
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
import {addTileViewerWith, deleteValueField, designerValueFields} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, noErrors, readingIs, readingReads, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu, readingContains, readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer — a card designed on one table, then rebound to another", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-designed-rebind.feature", import.meta.url);
  test("The designed spgi card gives way to demog's own card when the table changes", {tag: ["@viewers", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(20, "And user switches to the \"spgi-100\" table view", () => switchTableView(page, "spgi-100"));
    await session.step(21, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","Stereo Category"]]), [["Lanes Column Name","Stereo Category"]]);
    await session.step(23, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
    await session.step(24, "And the \"lane names\" reading of tile viewer should be \"R_ONE, S_ABS, S_ACHIR, S_PART, S_UNKN\"", () => readingReads(page, "lane names", el("tile viewer"), "R_ONE, S_ABS, S_ACHIR, S_PART, S_UNKN"));
    await session.step(25, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
    await session.step(26, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
    await session.step(27, "And the form designer should show value fields \"Idea-Created, Primary-Series-Name, Scaffold-Names\"", () => designerValueFields(page, "Idea-Created, Primary-Series-Name, Scaffold-Names"));
    await session.step(28, "When user deletes the \"Scaffold-Names\" value field in the form designer", () => deleteValueField(page, "Scaffold-Names"));
    await session.step(29, "And user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
    await session.step(30, "Then form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
    await session.step(31, "And the \"form designed\" reading of tile viewer should be \"true\"", () => readingReads(page, "form designed", el("tile viewer"), "true"));
    await session.step(32, "And the \"fields shown\" reading of tile viewer should be 2", () => readingIs(page, "fields shown", el("tile viewer"), 2));
    await session.step(33, "And the \"fields\" reading of tile viewer should contain \"Idea Created\"", () => readingContains(page, "fields", el("tile viewer"), "Idea Created"));
    await session.step(34, "And the \"fields\" reading of tile viewer should not contain \"Scaffold Names\"", () => readingNotContains(page, "fields", el("tile viewer"), "Scaffold Names"));
    await session.step(35, "When user sets \"Table\" property of tile viewer to \"demog-1000\"", () => setProperty(page, "Table", el("tile viewer"), "demog-1000"));
    await session.step(36, "Then tile viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("tile viewer"), "demog-1000"));
    await session.step(37, "And the \"table\" reading of tile viewer should be \"demog-1000\"", () => readingReads(page, "table", el("tile viewer"), "demog-1000"));
    await session.step(38, "And tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(39, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
    await session.step(40, "And the \"fields\" reading of tile viewer should contain \"SEX\"", () => readingContains(page, "fields", el("tile viewer"), "SEX"));
    await session.step(41, "And the \"fields\" reading of tile viewer should not contain \"Primary Series Name\"", () => readingNotContains(page, "fields", el("tile viewer"), "Primary Series Name"));
    await session.step(42, "And the \"fields\" reading of tile viewer should not contain \"Idea Created\"", () => readingNotContains(page, "fields", el("tile viewer"), "Idea Created"));
    await session.step(43, "And the \"SEX of row 1\" reading of tile viewer should be \"F\"", () => readingReads(page, "SEX of row 1", el("tile viewer"), "F"));
    await session.step(44, "And the \"AGE of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "26"));
    await session.step(45, "And the \"RACE of row 1\" reading of tile viewer should be \"Caucasian\"", () => readingReads(page, "RACE of row 1", el("tile viewer"), "Caucasian"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
});

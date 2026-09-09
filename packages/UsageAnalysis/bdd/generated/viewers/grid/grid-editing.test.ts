/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-editing.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clipboardContains, clipboardHas, pressKey, pressKeyIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnIncomplete, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, noneSelected, rowCount, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, doubleClickArea, eventFired, listenFor, noErrors, propertyShouldBe, setProperty, showsRows, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid cell editing and the clipboard", () => {
  const session = feature(test, "features/viewers/grid/grid-editing.feature", import.meta.url);
  test("Grid cell editing and the clipboard", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(14, "And \"Allow Edit\" property of grid should be \"true\"", () => propertyShouldBe(page, "Allow Edit", el("grid"), "true"));
    await session.step(15, "And cell editor should be hidden", () => shouldBe(page, el("cell editor"), "hidden"));
    await run.scenario("Enter commits the edited value and Escape cancels the next edit", async () => {
      await session.step(18, "Given user listens for \"d4-grid-cell-value-edited\" event on grid", () => listenFor(page, "d4-grid-cell-value-edited", el("grid")));
      await session.step(19, "When user double-clicks on the \"cell 4 of AGE\" area of grid", () => doubleClickArea(page, "cell 4 of AGE", el("grid")));
      await session.step(20, "Then cell editor should be visible", () => shouldBe(page, el("cell editor"), "visible"));
      await session.step(21, "When user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(22, "And user types \"99\" into cell editor", () => typeInto(page, "99", el("cell editor")));
      await session.step(23, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(24, "Then \"d4-grid-cell-value-edited\" event should have fired on grid", () => eventFired(page, "d4-grid-cell-value-edited", el("grid")));
      await session.step(25, "And the value of \"AGE\" column in row 4 should be \"99\"", () => valueInRow(page, "AGE", 4, "99"));
      await session.step(26, "And cell editor should be hidden", () => shouldBe(page, el("cell editor"), "hidden"));
      await session.step(27, "When user double-clicks on the \"cell 4 of AGE\" area of grid", () => doubleClickArea(page, "cell 4 of AGE", el("grid")));
      await session.step(28, "Then cell editor should be visible", () => shouldBe(page, el("cell editor"), "visible"));
      await session.step(29, "When user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(30, "And user types \"77\" into cell editor", () => typeInto(page, "77", el("cell editor")));
      await session.step(31, "And user presses Escape", () => pressKey(page, "Escape"));
      await session.step(32, "Then the value of \"AGE\" column in row 4 should be \"99\"", () => valueInRow(page, "AGE", 4, "99"));
      await session.step(33, "And cell editor should be hidden", () => shouldBe(page, el("cell editor"), "hidden"));
      await session.step(34, "When user double-clicks on the \"cell 4 of AGE\" area of grid", () => doubleClickArea(page, "cell 4 of AGE", el("grid")));
      await session.step(35, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(36, "And user types \"45\" into cell editor", () => typeInto(page, "45", el("cell editor")));
      await session.step(37, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(38, "Then the value of \"AGE\" column in row 4 should be \"45\"", () => valueInRow(page, "AGE", 4, "45"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Delete clears the current cell without opening the editor", async () => {
      await session.step(42, "When user clicks on the \"cell 5 of AGE\" area of grid", () => clickArea(page, "cell 5 of AGE", el("grid")));
      await session.step(43, "And user presses Delete", () => pressKey(page, "Delete"));
      await session.step(44, "Then \"AGE\" column should have missing values", () => columnIncomplete(page, "AGE"));
      await session.step(45, "And cell editor should be hidden", () => shouldBe(page, el("cell editor"), "hidden"));
      await session.step(46, "When user double-clicks on the \"cell 5 of AGE\" area of grid", () => doubleClickArea(page, "cell 5 of AGE", el("grid")));
      await session.step(47, "And user types \"51\" into cell editor", () => typeInto(page, "51", el("cell editor")));
      await session.step(48, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(49, "Then the value of \"AGE\" column in row 5 should be \"51\"", () => valueInRow(page, "AGE", 5, "51"));
      await session.step(50, "And \"AGE\" column should have no missing values", () => columnComplete(page, "AGE"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A read-only grid refuses the edit and warns", async () => {
      await session.step(54, "When user sets \"Allow Edit\" property of grid to \"false\"", () => setProperty(page, "Allow Edit", el("grid"), "false"));
      await session.step(55, "And user double-clicks on the \"cell 6 of AGE\" area of grid", () => doubleClickArea(page, "cell 6 of AGE", el("grid")));
      await session.step(56, "Then a warning balloon containing \"read-only\" should have been shown", () => warningBalloonText(page, "read-only"));
      await session.step(57, "And cell editor should be hidden", () => shouldBe(page, el("cell editor"), "hidden"));
      await session.step(58, "When user sets \"Allow Edit\" property of grid to \"true\"", () => setProperty(page, "Allow Edit", el("grid"), "true"));
      await session.step(59, "And user clicks on the \"cell 6 of AGE\" area of grid", () => clickArea(page, "cell 6 of AGE", el("grid")));
      await session.step(60, "And user presses 7", () => pressKey(page, "7"));
      await session.step(61, "Then cell editor should be visible", () => shouldBe(page, el("cell editor"), "visible"));
      await session.step(62, "When user presses Enter", () => pressKey(page, "Enter"));
      await session.step(63, "Then the value of \"AGE\" column in row 6 should be \"7\"", () => valueInRow(page, "AGE", 6, "7"));
      await session.step(64, "When user double-clicks on the \"cell 6 of AGE\" area of grid", () => doubleClickArea(page, "cell 6 of AGE", el("grid")));
      await session.step(65, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(66, "And user types \"49\" into cell editor", () => typeInto(page, "49", el("cell editor")));
      await session.step(67, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(68, "Then the value of \"AGE\" column in row 6 should be \"49\"", () => valueInRow(page, "AGE", 6, "49"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control+Shift+C copies the current cell", async () => {
      await session.step(72, "When user clicks on the \"cell 4 of AGE\" area of grid", () => clickArea(page, "cell 4 of AGE", el("grid")));
      await session.step(73, "And user presses Control+Shift+C", () => pressKey(page, "Control+Shift+C"));
      await session.step(74, "Then the clipboard should have the text \"45\"", () => clipboardHas(page, "45"));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control+C copies the selected rows", async () => {
      await session.step(78, "When user clicks on the \"row header 1\" area of grid", () => clickArea(page, "row header 1", el("grid")));
      await session.step(79, "And user clicks on the \"row header 5\" area of grid holding Shift", () => clickAreaHolding(page, "row header 5", el("grid"), "Shift"));
      await session.step(80, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
      await session.step(81, "When user presses Control+C", () => pressKey(page, "Control+C"));
      await session.step(82, "Then the clipboard should contain the text \"X0273T21000300003\"", () => clipboardContains(page, "X0273T21000300003"));
      await session.step(83, "And the clipboard should contain the text \"X0273T21000500006\"", () => clipboardContains(page, "X0273T21000500006"));
      await session.step(84, "When user clears the row selection", () => clearSelection(page));
      await session.step(85, "Then no rows should be selected", () => noneSelected(page));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Select all, copy and paste leaves the table as it was", async () => {
      await session.step(89, "When user clicks on the \"cell 1 of USUBJID\" area of grid", () => clickArea(page, "cell 1 of USUBJID", el("grid")));
      await session.step(90, "And user presses Control+A", () => pressKey(page, "Control+A"));
      await session.step(91, "And user presses Control+C", () => pressKey(page, "Control+C"));
      await session.step(92, "And user presses Control+V", () => pressKey(page, "Control+V"));
      await session.step(93, "Then the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(94, "And the value of \"USUBJID\" column in row 1 should be \"X0273T21000300003\"", () => valueInRow(page, "USUBJID", 1, "X0273T21000300003"));
      await session.step(95, "And the value of \"AGE\" column in row 1 should be \"26\"", () => valueInRow(page, "AGE", 1, "26"));
      await session.step(96, "When user clears the row selection", () => clearSelection(page));
      await session.step(97, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control+V pastes a copied value into another cell", async () => {
      await session.step(100, "When user clicks on the \"cell 2 of AGE\" area of grid", () => clickArea(page, "cell 2 of AGE", el("grid")));
      await session.step(101, "And user presses Control+C", () => pressKey(page, "Control+C"));
      await session.step(102, "And user clicks on the \"cell 11 of AGE\" area of grid", () => clickArea(page, "cell 11 of AGE", el("grid")));
      await session.step(103, "And user presses Control+V", () => pressKey(page, "Control+V"));
      await session.step(104, "Then the value of \"AGE\" column in row 11 should be \"30\"", () => valueInRow(page, "AGE", 11, "30"));
      await session.step(105, "When user double-clicks on the \"cell 11 of AGE\" area of grid", () => doubleClickArea(page, "cell 11 of AGE", el("grid")));
      await session.step(106, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(107, "And user types \"46\" into cell editor", () => typeInto(page, "46", el("cell editor")));
      await session.step(108, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(109, "Then the value of \"AGE\" column in row 11 should be \"46\"", () => valueInRow(page, "AGE", 11, "46"));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shift+Delete removes the selected rows and Control+Z brings them back", async () => {
      await session.step(113, "When user clicks on the \"row header 4\" area of grid holding Control", () => clickAreaHolding(page, "row header 4", el("grid"), "Control"));
      await session.step(114, "And user clicks on the \"row header 6\" area of grid holding Control", () => clickAreaHolding(page, "row header 6", el("grid"), "Control"));
      await session.step(115, "And user clicks on the \"row header 8\" area of grid holding Control", () => clickAreaHolding(page, "row header 8", el("grid"), "Control"));
      await session.step(116, "Then 3 rows should be selected", () => selectedRowCount(page, 3));
      await session.step(117, "When user presses Shift+Delete", () => pressKey(page, "Shift+Delete"));
      await session.step(118, "Then the table should have 997 rows", () => rowCount(page, 997));
      await session.step(119, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(120, "Then the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(121, "When user clears the row selection", () => clearSelection(page));
      await session.step(122, "Then no rows should be selected", () => noneSelected(page));
      await session.step(123, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});

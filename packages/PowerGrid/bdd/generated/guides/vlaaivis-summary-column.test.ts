/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/vlaaivis-summary-column.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/vlaaivis.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragTo, hoverOver, selectIn, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, contextPanelOpen, dialogCloses, openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hasArea, pickColorSwatch, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragColumnToEdge, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Score compounds with a VlaaiVis summary column", () => {
  const session = feature(test, "features/guides/vlaaivis-summary-column.feature", import.meta.url);
  test("Add a VlaaiVis column, group its columns into sectors and tune their scores", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(22, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(24, "When user picks \"Add > Summary Columns > VlaaiVis\" from the context menu of the \"cell 2 of Id\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > VlaaiVis", "cell 2 of Id", el("grid")), undefined, "Right-click any cell and pick Add > Summary Columns > VlaaiVis");
    await session.step(26, "Then grid should have a \"header VlaaiVis\" area", () => hasArea(page, el("grid"), "header VlaaiVis"), undefined, "The VlaaiVis column is added after the last column");
    await session.step(28, "When user drags the \"header VlaaiVis\" area of grid to the first position", () => dragColumnToEdge(page, "header VlaaiVis", el("grid"), "first"), undefined, "Drag its header onto the arrow that shows left of the grid: it becomes the first column");
    await session.step(30, "And user clicks on the \"header VlaaiVis\" area of grid", () => clickArea(page, "header VlaaiVis", el("grid")), undefined, "Click the VlaaiVis header: the Context Panel shows its settings");
    await session.step(32, "Then \"CAST Idea ID\" VlaaiVis sector should be visible", () => shouldBe(page, el("\"CAST Idea ID\" VlaaiVis sector"), "visible"), undefined, "It starts with a sector for each of the first three numeric columns");
    await session.step(34, "When user clicks on Columns input in context panel", () => clickOn(page, el("Columns input in context panel")), undefined, "Click Columns to choose the columns it scores");
    await session.step(36, "And user toggles the \"CAST Idea ID\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "CAST Idea ID", el("\"Select columns...\" dialog")), undefined, "Uncheck CAST Idea ID");
    await session.step(38, "And user toggles the \"Chemical Space X\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Chemical Space X", el("\"Select columns...\" dialog")), undefined, "Uncheck Chemical Space X");
    await session.step(40, "And user toggles the \"Chemical Space Y\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Chemical Space Y", el("\"Select columns...\" dialog")), undefined, "Uncheck Chemical Space Y");
    await session.step(42, "And user toggles the \"Board ID\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Board ID", el("\"Select columns...\" dialog")), undefined, "Uncheck Board ID");
    await session.step(44, "And user toggles the \"Competition assay\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Competition assay", el("\"Select columns...\" dialog")), undefined, "Check Competition assay");
    await session.step(46, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")), undefined, "Click OK");
    await session.step(48, "Then Columns input in context panel should contain text \"(7) Average Mass, TPSA\"", () => shouldContainText(page, el("Columns input in context panel"), "(7) Average Mass, TPSA"), undefined, "Seven columns are scored now; the ones no sector holds are listed under Unassigned columns");
    await session.step(50, "When user hovers over header of \"CAST Idea ID\" VlaaiVis sector", () => hoverOver(page, el("header of \"CAST Idea ID\" VlaaiVis sector")), undefined, "Point at the CAST Idea ID sector to show its buttons");
    await session.step(52, "And user clicks on \"Rename sector\" icon in \"CAST Idea ID\" VlaaiVis sector", () => clickOn(page, el("\"Rename sector\" icon in \"CAST Idea ID\" VlaaiVis sector")), undefined, "Click the pencil to rename the sector");
    await session.step(54, "And user types \"Physchem\" into Name input in \"Rename sector\" dialog", () => typeInto(page, "Physchem", el("Name input in \"Rename sector\" dialog")), undefined, "Name it Physchem");
    await session.step(56, "And user clicks on OK button in \"Rename sector\" dialog", () => clickOn(page, el("OK button in \"Rename sector\" dialog")), undefined, "Click OK");
    await session.step(58, "And user drags \"TPSA\" text in \"Unassigned columns\" VlaaiVis sector to \"Physchem\" VlaaiVis sector", () => dragTo(page, el("\"TPSA\" text in \"Unassigned columns\" VlaaiVis sector"), el("\"Physchem\" VlaaiVis sector")), undefined, "Drag TPSA from Unassigned columns onto Physchem");
    await session.step(60, "And user drags \"NIBR logP\" text in \"Unassigned columns\" VlaaiVis sector to \"Physchem\" VlaaiVis sector", () => dragTo(page, el("\"NIBR logP\" text in \"Unassigned columns\" VlaaiVis sector"), el("\"Physchem\" VlaaiVis sector")), undefined, "Drag NIBR logP onto Physchem too");
    await session.step(62, "Then \"NIBR logP\" VlaaiVis column in \"Physchem\" VlaaiVis sector should be visible", () => shouldBe(page, el("\"NIBR logP\" VlaaiVis column in \"Physchem\" VlaaiVis sector"), "visible"), undefined, "Both are wedges of the Physchem sector now");
    await session.step(64, "When user hovers over header of \"Chemical Space X\" VlaaiVis sector", () => hoverOver(page, el("header of \"Chemical Space X\" VlaaiVis sector")), undefined, "Point at the Chemical Space X sector");
    await session.step(66, "And user clicks on \"Rename sector\" icon in \"Chemical Space X\" VlaaiVis sector", () => clickOn(page, el("\"Rename sector\" icon in \"Chemical Space X\" VlaaiVis sector")), undefined, "Click the pencil");
    await session.step(68, "And user types \"Potency\" into Name input in \"Rename sector\" dialog", () => typeInto(page, "Potency", el("Name input in \"Rename sector\" dialog")), undefined, "Name it Potency");
    await session.step(70, "And user clicks on OK button in \"Rename sector\" dialog", () => clickOn(page, el("OK button in \"Rename sector\" dialog")), undefined, "Click OK");
    await session.step(72, "And user drags \"Competition assay\" text in \"Unassigned columns\" VlaaiVis sector to \"Potency\" VlaaiVis sector", () => dragTo(page, el("\"Competition assay\" text in \"Unassigned columns\" VlaaiVis sector"), el("\"Potency\" VlaaiVis sector")), undefined, "Drag Competition assay onto Potency");
    await session.step(74, "Then \"Competition assay\" VlaaiVis column in \"Potency\" VlaaiVis sector should be visible", () => shouldBe(page, el("\"Competition assay\" VlaaiVis column in \"Potency\" VlaaiVis sector"), "visible"), undefined, "Competition assay is Potency's wedge now");
    await session.step(76, "When user hovers over header of \"Chemical Space Y\" VlaaiVis sector", () => hoverOver(page, el("header of \"Chemical Space Y\" VlaaiVis sector")), undefined, "Point at the empty Chemical Space Y sector");
    await session.step(78, "And user clicks on \"Delete sector\" icon in \"Chemical Space Y\" VlaaiVis sector", () => clickOn(page, el("\"Delete sector\" icon in \"Chemical Space Y\" VlaaiVis sector")), undefined, "Click the trash can to delete it");
    await session.step(80, "Then \"Chemical Space Y\" VlaaiVis sector should be absent", () => shouldBe(page, el("\"Chemical Space Y\" VlaaiVis sector"), "absent"), undefined, "The sector is gone");
    await session.step(82, "When user clicks on \"+ Add sector\" link", () => clickOn(page, el("\"+ Add sector\" link")), undefined, "Click + Add sector for a new one");
    await session.step(84, "And user types \"Size\" into Name input in \"New sector\" dialog", () => typeInto(page, "Size", el("Name input in \"New sector\" dialog")), undefined, "Name it Size");
    await session.step(86, "And user clicks on OK button in \"New sector\" dialog", () => clickOn(page, el("OK button in \"New sector\" dialog")), undefined, "Click OK");
    await session.step(88, "And user drags \"Num Heavy Atoms\" text in \"Unassigned columns\" VlaaiVis sector to \"Size\" VlaaiVis sector", () => dragTo(page, el("\"Num Heavy Atoms\" text in \"Unassigned columns\" VlaaiVis sector"), el("\"Size\" VlaaiVis sector")), undefined, "Drag Num Heavy Atoms onto Size");
    await session.step(90, "Then \"Num Heavy Atoms\" VlaaiVis column in \"Size\" VlaaiVis sector should be visible", () => shouldBe(page, el("\"Num Heavy Atoms\" VlaaiVis column in \"Size\" VlaaiVis sector"), "visible"), undefined, "Size has its first column");
    await session.step(92, "When user clicks on color of \"Potency\" VlaaiVis sector", () => clickOn(page, el("color of \"Potency\" VlaaiVis sector")), undefined, "Click Potency's color square to change its color");
    await session.step(94, "And user picks the color \"#d62728\" in the color picker", () => pickColorSwatch(page, "#d62728"), undefined, "Pick red");
    await session.step(96, "And user clicks on OK button", () => clickOn(page, el("OK button")), undefined, "Click OK: Potency's wedges turn red");
    await session.step(98, "And user clicks on \"Settings\" icon in \"TPSA\" VlaaiVis column", () => clickOn(page, el("\"Settings\" icon in \"TPSA\" VlaaiVis column")), undefined, "Click TPSA's gear to edit its score function");
    await session.step(100, "And user selects \"gaussian\" in Mode input in \"TPSA\" dialog", () => selectIn(page, "gaussian", el("Mode input in \"TPSA\" dialog")), undefined, "Switch the curve to gaussian");
    await session.step(102, "And user types \"75\" into Mean input in \"TPSA\" dialog", () => typeInto(page, "75", el("Mean input in \"TPSA\" dialog")), undefined, "Make 75 the most desirable TPSA");
    await session.step(104, "And user types \"12\" into Sigma input in \"TPSA\" dialog", () => typeInto(page, "12", el("Sigma input in \"TPSA\" dialog")), undefined, "Make the bell narrower");
    await session.step(106, "And user clicks on OK button in \"TPSA\" dialog", () => clickOn(page, el("OK button in \"TPSA\" dialog")), undefined, "Click OK");
    await session.step(107, "Then the \"TPSA\" dialog should close", () => dialogCloses(page, "TPSA"));
    await session.step(109, "When user types \"0.5\" into \"NIBR logP\" VlaaiVis column", () => typeInto(page, "0.5", el("\"NIBR logP\" VlaaiVis column")), undefined, "Halve NIBR logP's weight: its wedge gets narrower");
    await session.step(111, "And user types \"0.6\" into \"Target min\" input in context panel", () => typeInto(page, "0.6", el("\"Target min\" input in context panel")), undefined, "Set Target min, the inner edge of the shaded ring");
    await session.step(113, "And user types \"0.8\" into \"Target max\" input in context panel", () => typeInto(page, "0.8", el("\"Target max\" input in context panel")), undefined, "and Target max, its outer edge");
    await session.step(115, "And user clicks on the \"cell 3 of VlaaiVis\" area of grid", () => clickArea(page, "cell 3 of VlaaiVis", el("grid")), undefined, "Click a VlaaiVis cell: the Context Panel lists the values behind it");
    await session.step(117, "Then context panel should contain text \"Physchem\"", () => shouldContainText(page, el("context panel"), "Physchem"), undefined, "The values are grouped by sector, with the unassigned columns last");
  });
});

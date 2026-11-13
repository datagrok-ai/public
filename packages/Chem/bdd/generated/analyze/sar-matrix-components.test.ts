/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/sar-matrix-components.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.sar-matrix]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {analysisFinished, matricesBuilt, readingIsTextOf} from '../../bindings/sar-matrix.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingAsRemembered, readingIs, readingNotAsRemembered, readingReads, rememberReading, setProperties, setProperty, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix Summary on component columns", () => {
  const session = feature(test, "features/analyze/sar-matrix-components.feature", import.meta.url);
  test("SAR Matrix Summary on component columns", {tag: ["@journey", "@realizes:chem.cp.sar-matrix"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And user opens the \"SarMatrixProtacDemo\" project", () => openProject(page, "SarMatrixProtacDemo"));
    await session.step(19, "And SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
    await run.scenario("A component's ranking selects and opens what each row counts", async () => {
      await session.step(22, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(23, "And user clicks on value of Warhead component row", () => clickOn(page, el("value of Warhead component row")));
      await session.step(24, "Then Summary evidence should contain text \"Warhead — effect vs the average compound\"", () => shouldContainText(page, el("Summary evidence"), "Warhead — effect vs the average compound"));
      await session.step(25, "And value of first evidence row should have text \"+1.65\"", () => shouldHaveText(page, el("value of first evidence row"), "+1.65"));
      await session.step(26, "And count badge of first evidence row should have text \"8 cpd\"", () => shouldHaveText(page, el("count badge of first evidence row"), "8 cpd"));
      await session.step(27, "When user clicks on count badge of first evidence row", () => clickOn(page, el("count badge of first evidence row")));
      await session.step(28, "Then 8 rows should be selected", () => selectedRowCount(page, 8));
      await session.step(29, "When user clicks on open icon of first evidence row", () => clickOn(page, el("open icon of first evidence row")));
      await session.step(30, "Then the \"tab\" reading of SAR Matrix Viewer viewer should be \"SAR Matrix\"", () => readingReads(page, "tab", el("SAR Matrix Viewer viewer"), "SAR Matrix"));
      await session.step(31, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(32, "Then the \"series\" reading of SAR Matrix Viewer viewer should be the text of name of first evidence row", () => readingIsTextOf(page, "series", el("SAR Matrix Viewer viewer"), el("name of first evidence row")));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A component's best swap selects the compounds with the new value", async () => {
      await session.step(36, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(37, "And user clicks on swap of Warhead component row", () => clickOn(page, el("swap of Warhead component row")));
      await session.step(38, "Then Summary evidence should contain text \"Swapping Warhead — what it bought\"", () => shouldContainText(page, el("Summary evidence"), "Swapping Warhead — what it bought"));
      await session.step(39, "And first evidence row should contain text \"· 4 pairs\"", () => shouldContainText(page, el("first evidence row"), "· 4 pairs"));
      await session.step(40, "And value of first evidence row should have text \"+4.14\"", () => shouldHaveText(page, el("value of first evidence row"), "+4.14"));
      await session.step(41, "When user clicks on open icon of first evidence row", () => clickOn(page, el("open icon of first evidence row")));
      await session.step(42, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Lower is better re-reads every ranking without a rebuild", async () => {
      await session.step(46, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(47, "And user remembers the \"key series\" reading of SAR Matrix Viewer viewer", () => rememberReading(page, "key series", el("SAR Matrix Viewer viewer")));
      await session.step(48, "And user sets \"Activity direction\" property of SAR Matrix Viewer viewer to \"Lower is better\"", () => setProperty(page, "Activity direction", el("SAR Matrix Viewer viewer"), "Lower is better"));
      await session.step(49, "Then value of Warhead component row should contain text \"−2.45\"", () => shouldContainText(page, el("value of Warhead component row"), "−2.45"));
      await session.step(50, "And the \"key series\" reading of SAR Matrix Viewer viewer should not be as remembered", () => readingNotAsRemembered(page, "key series", el("SAR Matrix Viewer viewer")));
      await session.step(51, "And the \"restored\" reading of SAR Matrix Viewer viewer should be \"true\"", () => readingReads(page, "restored", el("SAR Matrix Viewer viewer"), "true"));
      await session.step(52, "When user sets \"Activity direction\" property of SAR Matrix Viewer viewer to \"Higher is better\"", () => setProperty(page, "Activity direction", el("SAR Matrix Viewer viewer"), "Higher is better"));
      await session.step(53, "Then value of Warhead component row should contain text \"+1.65\"", () => shouldContainText(page, el("value of Warhead component row"), "+1.65"));
      await session.step(54, "And the \"key series\" reading of SAR Matrix Viewer viewer should be as remembered", () => readingAsRemembered(page, "key series", el("SAR Matrix Viewer viewer")));
      await session.step(55, "And the \"restored\" reading of SAR Matrix Viewer viewer should be \"true\"", () => readingReads(page, "restored", el("SAR Matrix Viewer viewer"), "true"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The dialog offers -lg again after an activity it cannot scale", async () => {
      await session.step(59, "When user picks \"Chem > Analyze > SAR Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > SAR Matrix..."));
      await session.step(60, "Then \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("\"SAR Matrix\" dialog"), "visible"));
      await session.step(61, "When user selects \"Solubility logS (pred)\" in Activity input in \"SAR Matrix\" dialog", () => selectIn(page, "Solubility logS (pred)", el("Activity input in \"SAR Matrix\" dialog")));
      await session.step(62, "Then Scaling input in \"SAR Matrix\" dialog should be disabled", () => shouldBe(page, el("Scaling input in \"SAR Matrix\" dialog"), "disabled"));
      await session.step(63, "And Scaling input in \"SAR Matrix\" dialog should have value \"none\"", () => shouldHaveValue(page, el("Scaling input in \"SAR Matrix\" dialog"), "none"));
      await session.step(64, "When user selects \"Average Mass\" in Activity input in \"SAR Matrix\" dialog", () => selectIn(page, "Average Mass", el("Activity input in \"SAR Matrix\" dialog")));
      await session.step(65, "Then Scaling input in \"SAR Matrix\" dialog should be enabled", () => shouldBe(page, el("Scaling input in \"SAR Matrix\" dialog"), "enabled"));
      await session.step(66, "And Scaling input in \"SAR Matrix\" dialog should have value \"-lg\"", () => shouldHaveValue(page, el("Scaling input in \"SAR Matrix\" dialog"), "-lg"));
      await session.step(67, "When user clicks on CANCEL button in \"SAR Matrix\" dialog", () => clickOn(page, el("CANCEL button in \"SAR Matrix\" dialog")));
      await session.step(68, "Then \"SAR Matrix\" dialog should be absent", () => shouldBe(page, el("\"SAR Matrix\" dialog"), "absent"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Another component across the columns rebuilds, a scaling it cannot take leaves none, and both go back at once", async () => {
      await session.step(72, "When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer", () => clickOn(page, el("SAR Matrix tab in SAR Matrix Viewer viewer")));
      await session.step(73, "And user selects \"E3 Ligand\" in Columns input in SAR Matrix Viewer viewer", () => selectIn(page, "E3 Ligand", el("Columns input in SAR Matrix Viewer viewer")));
      await session.step(74, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
      await session.step(75, "And the \"axis\" reading of SAR Matrix Viewer viewer should be \"E3 Ligand\"", () => readingReads(page, "axis", el("SAR Matrix Viewer viewer"), "E3 Ligand"));
      await session.step(76, "And the \"matrices\" reading of SAR Matrix Viewer viewer should be 71", () => readingIs(page, "matrices", el("SAR Matrix Viewer viewer"), 71));
      await session.step(77, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(78, "Then setup line should contain text \"columns: E3 Ligand\"", () => shouldContainText(page, el("setup line"), "columns: E3 Ligand"));
      await session.step(79, "And Warhead component row should be visible", () => shouldBe(page, el("Warhead component row"), "visible"));
      await session.step(80, "When user sets \"Scaling\" property of SAR Matrix Viewer viewer to \"lg\"", () => setProperty(page, "Scaling", el("SAR Matrix Viewer viewer"), "lg"));
      await session.step(81, "Then SAR Matrix Viewer viewer should have finished building", () => analysisFinished(page, el("SAR Matrix Viewer viewer")));
      await session.step(82, "And a warning balloon containing \"cannot be scaled by lg\" should have been shown", () => warningBalloonText(page, "cannot be scaled by lg"));
      await session.step(83, "And the \"matrices\" reading of SAR Matrix Viewer viewer should be 0", () => readingIs(page, "matrices", el("SAR Matrix Viewer viewer"), 0));
      await session.step(84, "And SAR Matrix note should have text \"No SAR matrices found: no activity value can be scaled by lg. Set Scaling to none.\"", () => shouldHaveText(page, el("SAR Matrix note"), "No SAR matrices found: no activity value can be scaled by lg. Set Scaling to none."));
      await session.step(85, "When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer", () => clickOn(page, el("SAR Matrix tab in SAR Matrix Viewer viewer")));
      await session.step(86, "Then SAR Matrix note should have text \"No SAR matrices found: no activity value can be scaled by lg. Set Scaling to none.\"", () => shouldHaveText(page, el("SAR Matrix note"), "No SAR matrices found: no activity value can be scaled by lg. Set Scaling to none."));
      await session.step(87, "When user clicks on SAR Transfer tab in SAR Matrix Viewer viewer", () => clickOn(page, el("SAR Transfer tab in SAR Matrix Viewer viewer")));
      await session.step(88, "Then SAR Matrix note should contain text \"No SAR matrices to compare\"", () => shouldContainText(page, el("SAR Matrix note"), "No SAR matrices to compare"));
      await session.step(89, "When user sets properties of SAR Matrix Viewer viewer:", () => setProperties(page, el("SAR Matrix Viewer viewer"), [["Scaling","none"],["Matrix columns","Warhead"]]), [["Scaling","none"],["Matrix columns","Warhead"]]);
      await session.step(92, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
      await session.step(93, "And the \"matrices\" reading of SAR Matrix Viewer viewer should be 70", () => readingIs(page, "matrices", el("SAR Matrix Viewer viewer"), 70));
      await session.step(94, "And the \"axis\" reading of SAR Matrix Viewer viewer should be \"Warhead\"", () => readingReads(page, "axis", el("SAR Matrix Viewer viewer"), "Warhead"));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});

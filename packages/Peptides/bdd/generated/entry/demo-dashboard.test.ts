/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/entry/demo-dashboard.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {analysisScaling} from '../../bindings/activity.js';
import {demoRegistration} from '../../bindings/demo.js';
import {openDemo, peptidesInitialized, sarReady, sarSetting} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {columnSemType, columnTag, columnUnits, hasColumn, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount, tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {areaAtLeastTall, areaColors, noBalloons, noErrors, painted, readingIs, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Peptide SAR demo dashboard", () => {
  const session = feature(test, "features/entry/demo-dashboard.feature", import.meta.url);
  test("The demo builds a working dashboard on Simple peptides", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(12, "Then the Peptide SAR demo should be registered as a Bioinformatics dashboard", () => demoRegistration(page));
    await session.step(13, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(14, "When user opens the Peptide SAR demo dashboard", () => openDemo(page));
    await session.step(15, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(16, "And table \"Simple peptides\" should be open", () => tableOpen(page, "Simple peptides"));
    await session.step(17, "And the table should have 647 rows", () => rowCount(page, 647));
    await session.step(18, "And \"AlignedSequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AlignedSequence", "Macromolecule"));
    await session.step(19, "And \"AlignedSequence\" column should have units \"fasta\"", () => columnUnits(page, "AlignedSequence", "fasta"));
    await session.step(20, "And \"AlignedSequence\" column should have tag \"alphabet\" equal to \"PT\"", () => columnTag(page, "AlignedSequence", "alphabet", "PT"));
    await session.step(21, "And \"AlignedSequence\" column should have tag \"aligned\" equal to \"SEQ.MSA\"", () => columnTag(page, "AlignedSequence", "aligned", "SEQ.MSA"));
    await session.step(22, "And the table should have a column \"15\"", () => hasColumn(page, "15"));
    await session.step(23, "And the table should not have a column \"16\"", () => hasNoColumn(page, "16"));
    await session.step(24, "And the \"header 2\" area of grid should be at least 100 pixels tall", () => areaAtLeastTall(page, "header 2", el("grid"), 100));
    await session.step(25, "And the \"header 2\" area of grid should be painted in at least 2 colors", () => areaColors(page, "header 2", el("grid"), 2));
    await session.step(26, "And the SAR setting \"activityScaling\" should be \"-lg\"", () => sarSetting(page, "activityScaling", "-lg"));
    await session.step(27, "And the SAR activity column should use \"-lg\" scaling", () => analysisScaling(page, "-lg"));
    await session.step(28, "And the SAR setting \"mclSettings.threshold\" should be \"94\"", () => sarSetting(page, "mclSettings.threshold", "94"));
    await session.step(29, "And the \"completed threshold\" reading of MCL viewer should be 94", () => readingIs(page, "completed threshold", el("MCL viewer"), 94));
    await session.step(30, "And the table should have a column \"Cluster (MCL)\"", () => hasColumn(page, "Cluster (MCL)"));
    await session.step(31, "And the open tableview should have 1 Sequence Variability Map viewer", () => viewerCount(page, 1, "Sequence Variability Map"));
    await session.step(32, "And the open tableview should have 1 Most Potent Residues viewer", () => viewerCount(page, 1, "Most Potent Residues"));
    await session.step(33, "And the open tableview should have 1 MCL viewer", () => viewerCount(page, 1, "MCL"));
    await session.step(34, "And the open tableview should have 1 Logo Summary Table viewer", () => viewerCount(page, 1, "Logo Summary Table"));
    await session.step(35, "And the \"positions\" reading of Sequence Variability Map viewer should be 15", () => readingIs(page, "positions", el("Sequence Variability Map viewer"), 15));
    await session.step(36, "And the \"activity scaling\" reading of Sequence Variability Map viewer should be \"-lg\"", () => readingReads(page, "activity scaling", el("Sequence Variability Map viewer"), "-lg"));
    await session.step(37, "And Sequence Variability Map viewer should be painted", () => painted(page, el("Sequence Variability Map viewer")));
    await session.step(38, "And Most Potent Residues viewer should be painted", () => painted(page, el("Most Potent Residues viewer")));
    await session.step(39, "And the \"positions\" reading of Most Potent Residues viewer should be 15", () => readingIs(page, "positions", el("Most Potent Residues viewer"), 15));
    await session.step(40, "And Logo Summary Table viewer should be painted", () => painted(page, el("Logo Summary Table viewer")));
    await session.step(41, "And the \"members total\" reading of Logo Summary Table viewer should be 647", () => readingIs(page, "members total", el("Logo Summary Table viewer"), 647));
    await session.step(42, "And scatter plot viewer in MCL viewer should be painted", () => painted(page, el("scatter plot viewer in MCL viewer")));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
    await session.step(44, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});

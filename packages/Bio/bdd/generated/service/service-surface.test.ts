/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/service/service-surface.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.service-surface]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {helmMonomersMatch} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {call, callWith, resultHasMethods, resultIsList, resultProperty} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The service surface other packages call", () => {
  const session = feature(test, "features/service/service-surface.feature", import.meta.url);
  test("The service surface other packages call", {tag: ["@journey", "@realizes:bio.service-surface"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await run.scenario("The singletons resolve with the methods their consumers use", async () => {
      await session.step(13, "When user calls \"Bio:getSeqHelper\" function", () => call(page, "Bio:getSeqHelper"));
      await session.step(14, "Then the result should have methods \"getSeqHandler, getSeqMonomers, helmToAtomicLevel, setUnitsToFastaColumn\"", () => resultHasMethods(page, "getSeqHandler, getSeqMonomers, helmToAtomicLevel, setUnitsToFastaColumn"));
      await session.step(15, "When user calls \"Bio:getMonomerLibHelper\" function", () => call(page, "Bio:getMonomerLibHelper"));
      await session.step(16, "Then the result should have methods \"getMonomerLib, awaitLoaded\"", () => resultHasMethods(page, "getMonomerLib, awaitLoaded"));
      await session.step(17, "When user calls \"Bio:getBioLib\" function", () => call(page, "Bio:getBioLib"));
      await session.step(18, "Then the result should have methods \"getMonomer, getMonomerSymbolsByType, getPolymerTypes\"", () => resultHasMethods(page, "getMonomer, getMonomerSymbolsByType, getPolymerTypes"));
      await session.step(19, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A sequence handler is per column and reports the column's notation", async () => {
      await session.step(22, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
      await session.step(23, "Then \"HELM string\" column should have units \"helm\"", () => columnUnits(page, "HELM string", "helm"));
      await session.step(24, "When user calls \"Bio:getSeqHandler\" function with:", () => callWith(page, "Bio:getSeqHandler", [["sequence","column:HELM string"]]), [["sequence","column:HELM string"]]);
      await session.step(26, "Then the result should have methods \"getSplitter, getRegion, convert\"", () => resultHasMethods(page, "getSplitter, getRegion, convert"));
      await session.step(27, "And the result should have a \"notation\" of \"helm\"", () => resultProperty(page, "notation", "helm"));
      await session.step(28, "Given user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
      await session.step(29, "When user calls \"Bio:getSeqHandler\" function with:", () => callWith(page, "Bio:getSeqHandler", [["sequence","column:fasta"]]), [["sequence","column:fasta"]]);
      await session.step(31, "Then the result should have a \"notation\" of \"fasta\"", () => resultProperty(page, "notation", "fasta"));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The HELM monomer list is exactly the monomers of the column", async () => {
      await session.step(35, "When user switches to the \"filter_HELM\" table view", () => switchTableView(page, "filter_HELM"));
      await session.step(36, "And user calls \"Bio:getHelmMonomers\" function with:", () => callWith(page, "Bio:getHelmMonomers", [["sequence","column:HELM string"]]), [["sequence","column:HELM string"]]);
      await session.step(38, "Then the result should be a list of 1 or more items", () => resultIsList(page, 1));
      await session.step(39, "And the result should be exactly the monomers of \"HELM string\" column", () => helmMonomersMatch(page, "HELM string"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});

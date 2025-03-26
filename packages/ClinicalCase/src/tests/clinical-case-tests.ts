import { category, test } from "@datagrok-libraries/utils/src/test";
import { _testBaselineEndpoint, _testCumulativeEnrollment, _testHysLaw, _testLabChangesBetweenMinMax, _testLabChangesComparedToBaseline, _testLabChangesRelatedToRef, _testLabValuesByVisit, _testOpenApp, _testPatientVisit, _testStudyVisit, _testSurvivalDataframe, _testValidation } from "./utils";

category('clinicalCase', () => {

    test('clinicalCase.appOpensWithountExceptions', async () => {
        await _testOpenApp();
    })

    test('clinicalCase.cumulativeEnrollment', async () => {
        await _testCumulativeEnrollment();
    })

    test('clinicalCase.validation', async () => {
        await _testValidation();
    })

    test('clinicalCase.hysLaw', async () => {
        await _testHysLaw();
    })

    test('clinicalCase.baselineEndpoint', async () => {
        await _testBaselineEndpoint();
    })

    test('clinicalCase.labValuesByVisit', async () => {
        await _testLabValuesByVisit();
    })

    test('clinicalCase.survivalDataframe', async () => {
        await _testSurvivalDataframe();
    })

    test('clinicalCase.labChangesComparedToBaseline', async () => {
        await _testLabChangesComparedToBaseline();
    })

    test('clinicalCase.labChangesBetweenMinMax', async () => {
        await _testLabChangesBetweenMinMax();
    })

    test('clinicalCase.labChangesRelatedToRef', async () => {
        await _testLabChangesRelatedToRef();
    })

    test('clinicalCase.studyVisit', async () => {
        await _testStudyVisit();
    })

    test('clinicalCase.patientVisit', async () => {
        await _testPatientVisit();
    })

});
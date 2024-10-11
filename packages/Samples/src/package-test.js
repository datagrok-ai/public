var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
import * as DG from 'datagrok-api/dg';
import './test/demo-project-test';
// import './test/projects-test';
import { runTests, tests, initAutoTests as initTests } from '@datagrok-libraries/utils/src/test';
export const _package = new DG.Package();
export { tests };
//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export function test(category, test, testContext) {
    return __awaiter(this, void 0, void 0, function* () {
        const data = yield runTests({ category, test, testContext });
        return DG.DataFrame.fromObjects(data);
    });
}
//name: initAutoTests
export function initAutoTests() {
    return __awaiter(this, void 0, void 0, function* () {
        yield initTests(_package, _package.getModule('package-test.js'));
    });
}
//# sourceMappingURL=package-test.js.map
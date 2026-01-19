import type * as _DG from 'datagrok-api/dg';
import { Observable } from 'rxjs';
export declare const tests: {
    [key: string]: Category;
};
export declare let currentCategory: string;
export declare namespace assure {
    function notNull(value: any, name?: string): void;
}
export interface TestOptions {
    timeout?: number;
    benchmarkWarnTimeout?: number;
    benchmarkTimeout?: number;
    unhandledExceptionTimeout?: number;
    skipReason?: string;
    isAggregated?: boolean;
    benchmark?: boolean;
    stressTest?: boolean;
    owner?: string;
    tags?: string[];
}
export interface TestResult {
    date: string;
    category: string;
    name: string;
    success: boolean;
    result: any;
    ms: number;
    skipped: boolean;
    logs: string;
    owner: string;
    package: string;
    flaking: boolean;
}
export interface TestResultExtended extends TestResult {
    widgetsDifference: number;
}
export interface CategoryOptions {
    clear?: boolean;
    timeout?: number;
    benchmarks?: boolean;
    stressTests?: boolean;
    owner?: string;
}
export declare class TestContext {
    stressTest?: boolean;
    catchUnhandled: boolean;
    report: boolean;
    returnOnFail: boolean;
    constructor(catchUnhandled?: boolean, report?: boolean, returnOnFail?: boolean);
}
export declare class Test {
    test: () => Promise<any>;
    name: string;
    category: string;
    options?: TestOptions;
    constructor(category: string, name: string, test: () => Promise<any>, options?: TestOptions);
}
export declare class Category {
    tests?: Test[];
    before?: () => Promise<void>;
    after?: () => Promise<void>;
    beforeStatus?: string;
    afterStatus?: string;
    clear?: boolean;
    timeout?: number;
    benchmarks?: boolean;
    benchmarkTimeout?: number;
    stressTests?: boolean;
    owner?: string;
}
export declare class NodeTestExecutionOptions {
    package: _DG.Package;
}
export declare class TestExecutionOptions {
    category?: string;
    test?: string;
    testContext?: TestContext;
    exclude?: string[];
    verbose?: boolean;
    stressTest?: boolean;
    tags?: string[];
    nodeOptions?: NodeTestExecutionOptions;
    skipToCategory?: string;
    skipToTest?: string;
    returnOnFail?: boolean;
}
export declare function testEvent<T>(event: Observable<T>, handler: (args: T) => void, trigger: () => void, ms?: number, reason?: string): Promise<any>;
export declare function testEventAsync<T>(event: Observable<T>, handler: (args: T) => Promise<void>, trigger: () => void, ms?: number, reason?: string): Promise<any>;
export declare function test(name: string, test: () => Promise<any>, options?: TestOptions): void;
export declare function expect(actual: any, expected?: any, error?: string): void;
export declare function expectFloat(actual: number, expected: number, tolerance?: number, error?: string): void;
export declare function expectTable(actual: _DG.DataFrame, expected: _DG.DataFrame, error?: string): void;
export declare function expectObject(actual: {
    [key: string]: any;
}, expected: {
    [key: string]: any;
}): void;
export declare function expectArray(actual: ArrayLike<any>, expected: ArrayLike<any>): void;
export declare function category(category: string, tests_: () => void, options?: CategoryOptions): void;
export declare function before(before: () => Promise<void>): void;
export declare function after(after: () => Promise<void>): void;
export declare function initAutoTests(package_: _DG.Package, module?: any): Promise<void>;
export declare function runTests(options?: TestExecutionOptions): Promise<TestResultExtended[]>;
export declare function shuffle(array: any[]): any[];
export declare function delay(ms: number): Promise<void>;
export declare function awaitCheck(checkHandler: () => boolean, error?: string, wait?: number, interval?: number): Promise<any>;
export declare function timeout(func: () => Promise<any>, testTimeout: number, timeoutReason?: string): Promise<any>;
export declare function isDialogPresent(dialogTitle: string): boolean;
/** Expects an asynchronous {@link action} to throw an exception. Use {@link check} to perform
 * deeper inspection of the exception if necessary.
 * @param  {function(): Promise<void>} action
 * @param  {function(any): boolean} check
 * @return {Promise<void>}
 */
export declare function expectExceptionAsync(action: () => Promise<void>, check?: (exception: any) => boolean): Promise<void>;
/**
 * Universal test for viewers. It search viewers in DOM by tags: canvas, svg, img, input, h1, a
 * @param  {string} v Viewer name
 * @param  {_DG.DataFrame} df Dataframe to use. Should have at least 3 rows
 * @param  {boolean} options.detectSemanticTypes Specify whether to detect semantic types or not
 * @param  {boolean} options.readOnly If set to true, the dataframe will not be modified during the test
 * @param  {boolean} options.arbitraryDfTest If set to false, test on arbitrary dataframe
 * (one categorical column) will not be performed
 * @param  {object} options List of options (optional)
 * @return {Promise<void>} The test is considered successful if it completes without errors
 */
export declare function testViewer(v: string, df: _DG.DataFrame, options?: {
    detectSemanticTypes?: boolean;
    readOnly?: boolean;
    arbitraryDfTest?: boolean;
    packageName?: string;
    awaitViewer?: (viewer: _DG.Viewer) => Promise<void>;
}): Promise<void>;
//# sourceMappingURL=test.d.ts.map
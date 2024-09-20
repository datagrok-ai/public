export function setRandomOrder(tests: any[], workersAmount: number, testRepeats: number): any[][] {
    if (workersAmount > tests.length)
        workersAmount = tests.length;

    const repeatedTests = repeatTests(tests, testRepeats);
    const orderedTests = shuffle(repeatedTests);
    return splitArray(orderedTests, workersAmount);
}

export function setPackageRandomOrder(tests: any[], workersAmount: number, testRepeats: number): any[][] {
    const packages = splitTestsByPackages(tests);

    if (workersAmount > packages.size)
        workersAmount = packages.size;

    for (let packageName of packages.keys()) {
        packages.set(packageName, repeatTests(packages.get(packageName), testRepeats));
        packages.set(packageName, shuffle(packages.get(packageName)));
    }

    const splittedPackagesArray = splitArray([...packages.values()], workersAmount);
    const splittedTestsArray: any[] = [];
    
    for (let splitArray of splittedPackagesArray) {
        let testsArray: any = [];
        for (let packageTests of splitArray) {
            testsArray = testsArray.concat(packageTests);
        }
        splittedTestsArray.push(testsArray);
    } 
    
    return splittedTestsArray;
}

export function setAlphabeticalOrder(tests: any[], workersAmount: number, testRepeats: number): any[][] {
    if (workersAmount > tests.length)
        workersAmount = tests.length;

    const repeatedTests = repeatTests(tests, testRepeats);
    const orderedTests = repeatedTests.sort((test1: any, test2: any) => { return (`${test1.category}: ${test1.name}s`.localeCompare(`${test2.category}: ${test2.name}s`)) });

    return splitArray(orderedTests, workersAmount);
}

export function setPackageAlphabeticalOrder(tests: any[], workersAmount: number, testRepeats: number): any[][] {
    const packages = splitTestsByPackages(tests);

    if (workersAmount > packages.size)
        workersAmount = packages.size;

    for (let packageName of packages.keys()) {
        packages.set(packageName, repeatTests(packages.get(packageName), testRepeats));
        packages.set(packageName, packages.get(packageName).sort(alphabeticalSortFunction));
    }

    const splittedPackagesArray = splitArray([...packages.values()], workersAmount);
    const splittedTestsArray: any[] = [];

    for (let splitArray of splittedPackagesArray) {
        let testsArray: any = [];
        for (let packageTests of splitArray) {
            testsArray = testsArray.concat(packageTests);
        }
        splittedTestsArray.push(testsArray);
    } 

    return splittedTestsArray;
}

function alphabeticalSortFunction(test1: any, test2: any): number {
    return (`${test1.category}: ${test1.name}s`.localeCompare(`${test2.category}: ${test2.name}s`));
}

function splitTestsByPackages(tests: any[]): Map<string, any> {
    let resultMap = new Map<string, any>();

    for (let test of tests) {
        if (resultMap.has(test.packageName))
            resultMap.get(test.packageName).push(test);
        else
            resultMap.set(test.packageName, [test])
    }

    return resultMap;
}

function repeatTests(tests: Object[], testRepeats: number): Object[] {
    let repeatedTests: any[] = []
    for (let test of tests) {
        for (let i = 0; i < testRepeats; i++) {
            repeatedTests.push({ ...test })
        }
    }
    return repeatedTests;
}

function shuffle<T>(array: T[]): T[] {
    const newArr = array.slice();
    newArr.sort(() => Math.random() - 0.5);
    return newArr;
};

function splitArray<T>(arr: T[], countOfParts: number): T[][] {
    const result: T[][] = [];
    const partSize = Math.ceil(arr.length / countOfParts);

    for (let i = 0; i < arr.length; i += partSize) {
        result.push(arr.slice(i, i + partSize));
    }

    return result;
}
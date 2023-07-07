import * as DG from 'datagrok-api/dg';
import {DEMO_APP_HIERARCHY} from './const';


type Category = {
	name: string;
	children?: Category[];
};

function findCategoryNumber(categoryName: string, obj: Category[] = DEMO_APP_HIERARCHY.children): number {
	let result = -1;
	for (let i = 0; i < obj.length; i++) {
		if (obj[i].name === categoryName)
			return i;
		if (obj[i].children)
		  result = findCategoryNumber(categoryName, obj[i].children)
		if (result !== -1)
			return result;
	}
	return -1;
}

export function sortFunctionsByHierarchy(func1: DG.Func, func2: DG.Func): number {
	const func1DemoPath = getDemoPath(func1);
	const func2DemoPath = getDemoPath(func2);

	const length = Math.min(func1DemoPath.length, func2DemoPath.length) - 1;
	let func1CategoryNumber = 0;
	let func2CategoryNumber = 0;

	for (let i = 0; i < length; ++i) {
		func1CategoryNumber = findCategoryNumber(func1DemoPath[i]);
		func2CategoryNumber = findCategoryNumber(func2DemoPath[i]);

		if (func1CategoryNumber === func2CategoryNumber)
			continue;
		return func1CategoryNumber > func2CategoryNumber ? 1: -1;
	}

	if (func1DemoPath.length !== func2DemoPath.length) {
		return func1DemoPath.length > func2DemoPath.length ? -1 : 1;
	}

	func1CategoryNumber = findCategoryNumber(func1DemoPath[func1DemoPath.length - 1]);
	func2CategoryNumber = findCategoryNumber(func2DemoPath[func2DemoPath.length - 1]);
	return func1CategoryNumber > func2CategoryNumber ? 1: func1CategoryNumber < func2CategoryNumber ? -1 : 0;
}

function getDemoPath(func: DG.Func): string[] {
	return func.options[DG.FUNC_OPTIONS.DEMO_PATH].split('|').map((elem: string) => elem.trim());
}

export function getParentCategoryName(categoryName: string, parentName: string = '',
	obj: Category[] = DEMO_APP_HIERARCHY.children): string {
	let result = '';
	for (let i = 0; i < obj.length; i++) {
		if (obj[i].name.toLowerCase() === categoryName.toLowerCase())
			return parentName;
		if (obj[i].children)
			result = getParentCategoryName(categoryName, obj[i].name, obj[i].children);
		if (result !== '')
			return result;
	}
	return '';
}

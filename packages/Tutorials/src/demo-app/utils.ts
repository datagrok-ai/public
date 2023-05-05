import * as DG from 'datagrok-api/dg';


export function sortFunctionsByDemoPath(func1: DG.Func, func2: DG.Func): number {
	const func1Name = getDemoPath(func1);
	const func2Name = getDemoPath(func2);

	if (func1Name.length !== func2Name.length)
		return func1Name.length > func2Name.length ? 1 : -1;
	
	for (let i = 0; i < func1Name.length; i++) {
		if (func1Name[i] !== func2Name[i])
			return func1Name[i] > func2Name[i] ? 1 : -1;
	}

	return 0;
}

export function getDemoPath(func: DG.Func): string[] {
	return func.options[DG.FUNC_OPTIONS.DEMO_PATH].split('|').map((elem: string) => elem.trim());
}

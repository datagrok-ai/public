//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//output: dataframe result
export async function test(category: string, test: string): Promise<any> {
  return null;
}

//name: helper
//output: int b
export function helperShadow(): number {
  return 1;
}

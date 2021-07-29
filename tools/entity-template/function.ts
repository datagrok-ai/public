
//name: #{NAME}
//input: string name
//output: string greeting
export function #{NAME}(name: string) {
  let greeting = 'Hello, ' + name;
  grok.shell.info(greeting);
  return greeting;
}

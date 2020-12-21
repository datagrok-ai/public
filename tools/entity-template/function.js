//name: #{NAME}
//input: string name
//output: string greeting
export function #{NAME}(name) {
  let greeting = 'Hello, ' + name;
  grok.shell.info(greeting);
  return greeting;
}

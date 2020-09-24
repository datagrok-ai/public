
//name: #{NAME}
//input: string name
//output: string greeting
export function #{NAME_LOWERCASE}(name) {
    let greeting = 'Hello, ' + name;
    grok.shell.info(greeting);
    return greeting;
}

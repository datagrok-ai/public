// A form that only renders and binds (skipLogic): the host decides items, visibility, enabled state and validity.
// Without skipLogic, the `visible:` expression would hide "Tank volume" again on the next edit.

const script = DG.Script.create(`//name: CarForm
//language: javascript
//input: string type = "ICE" {choices: ["ICE", "Electric"]}
//input: double tankVolume = 40 {visible: type == "ICE"; min: 10; max: 100}
//input: double batteryCapacity = 80 {visible: type == "Electric"}
//output: string res
res = type;`);

const fc = script.prepare({type: 'Hybrid', tankVolume: 500, batteryCapacity: 20});
const form = await DG.InputForm.forFuncCall(fc, {skipLogic: true});

const typeInput = form.getInput('type');
typeInput.items = ['ICE', 'Electric', 'Hybrid'];
typeInput.value = fc.getParamValue('type');

const applyHostRules = () => {
  const type = fc.getParamValue('type');
  form.getInput('tankVolume').root.style.display = type === 'Electric' ? 'none' : '';
  form.getInput('batteryCapacity').enabled = type !== 'ICE';
};
applyHostRules();
form.onInputChanged.subscribe(() => applyHostRules());

grok.shell.newView('skipLogic', [form.root]);

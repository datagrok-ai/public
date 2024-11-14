
//input: string type = 'ICE' { choices: ['Electric', 'ICE'] }
//input: int cylinders = 4 { visible: type == 'ICE' }
//input: double tankVolume = 40 { visible: type == 'ICE'; units: liters }
//input: bool tankExtension = false { visible: type == 'ICE'; enabled: tankVolume > 50 }
//input: double batteryCapacity = 80 { visible: type == 'Electric'; units: kWh }

grok.shell.info(type);
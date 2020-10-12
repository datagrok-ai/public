<!-- TITLE: Grok Script -->
<!-- SUBTITLE: -->

# Grok Script

Grok script language is used to control or automate everything within 
the Datagrok platform. Use it to transform data, automate workflows, run queries,
evaluate numerical expressions, execute commands, record macros, 
perform statistical computations, execute R scripts.

## Recording macros

  1. Open console (**Tools | Console**) 
  2. Do the work that needs to be automated. Corresponding commands will appear in the console
  3. Select corresponding commands in the console, and copy to clipboard 

## Executing scripts

Type a command into the console and press Enter.

## Functions

Open **Help | Functions** to see a list of available [functions](functions/function.md).

## Syntax

Grok script supports all common features and operators, such as:

### Variables

Assigning a variable:

```
x = 5 
```

### Math operators

Supported operators: +, -, *, /, ^ (power), % (modulus)

```
a = (2 + 2) * 2^2
```

### Text operators
```
s = "Data" + "Grok"
```

### Logical operators
```
a = 2==2 || 2!=2
```

### Lists and maps
```
a = [1, 2, 3, {"a" : 1, "b" : 2}]
```

### Objects
```
p = Project()
p.name = "test"
x = p.name
```

### System variables and constants
```
x = PI
y = E
currentTable = t
currentView = v
```

### Methods and extensions
You can invoke any action by name, or, action can be called as method of variable with same 
type of first action parameter:
```
"test string".length()
length("test string")
```

```
KeepRows(t, Selected())
t.KeepRows(Selected())
```

### Table and Column names convention

Symbols '"', '{' and '}' should be replaced as '^^', '<\[' and ']>'.

You can run any system [action](functions/function.md) by calling it from [console](navigation.md#console)

## Try it!

Open [console](navigation.md#console) by pressing ~ (tilda) or **Tools | Console**. and try to make some actions: 
run [query](../access/data-query.md) or [job](../access/data-job.md).
Every step you take will be recorded, so you can re-run it, or use somewhere.

See also:

* [Console](navigation.md#console)
* [Functions](functions/function.md)
* [Scripting plugin](../develop/scripting.md)

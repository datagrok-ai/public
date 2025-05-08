import json
import os
from pathlib import Path
import ast
from pathlib import Path
import re
import sys

from docs.docstring_parser import parse_docstring


def sanitize(text):
    return text.replace('\n', ' ').strip()

def write_func_md(output_dir, func_data, classes: dict):
    func_name = func_data['name']
    output_file = Path(output_dir) / "functions" / f"{func_name}.md"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    func_doc = parse_docstring(func_data['docstring'] or '-')
    lines = [
        f"# {func_name}\n\n"
    ]
    for desc in func_doc['docstring']:
        lines.append(desc)
        lines.append('\n\n')
    params = func_doc["parameters"]
    returns = func_doc["returns"]
    classes_patern = r'\b(' + '|'.join(re.escape(t) for t in classes.keys()) + r')\b'
    if params:
        lines.append("**Parameters**\n\n")
        lines.append("| Name | Type | Description |\n")
        lines.append("| :--- | :--- | :---------- |\n")
        for p in params:
            param_type = re.sub(classes_patern, lambda m: f"[{m.group(0)}]({'../../' + classes.get(m.group(0))})", p['type'])
            param_type = classes.get(param_type, param_type)
            lines.append(f"| {p['name']} | {param_type} | {p['description']} |\n")
        lines.append("\n")

    if returns:
        lines.append("**Returns**\n\n")
        lines.append("| Type | Description |\n")
        lines.append("| :--- | :---------- |\n")
        for r in returns:
            return_type = re.sub(classes_patern, lambda m: f"[{m.group(0)}]({'../../' + classes.get(m.group(0))})", r['type'])
            return_type = classes.get(return_type, return_type)
            lines.append(f"| {return_type} | {r['description']} |\n")
        lines.append("\n")   
    output_file.write_text("".join(lines), encoding="utf-8")

def write_class_md(output_dir, class_data, classes: dict):
    class_name = class_data['name']
    output_file = Path(output_dir) / "classes" / f"{class_name}.md"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    class_docs = parse_docstring(class_data['docstring'] or '-')

    classes_patern = r'\b(' + '|'.join(re.escape(t) for t in classes.keys()) + r')\b'

    lines = [
        f"# {class_name}\n\n"
    ]

    for desc in class_docs['docstring']:
        lines.append(desc)
        lines.append('\n\n')

    if class_docs["attributes"] and len(class_docs["attributes"]) > 0:
        lines.append("## Attributes\n\n")
        lines.append("| Name | Type | Description |\n")
        lines.append("| :--- | :--- | :---------- |\n")
        for attr in class_docs["attributes"]:
            attr_type = re.sub(classes_patern, lambda m: f"[{m.group(0)}]({f'{class_name}.md' if class_name in classes.get(m.group(0)) else '../../' + classes.get(m.group(0))})", attr['type'])
            attr_type = classes.get(attr_type, attr_type)
            lines.append(f"| {attr['name']} | {attr_type} | {attr['description']} |\n")
        lines.append('\n')
        
    lines.append("## Methods\n\n")
    build_methods(class_data["methods"], lines, classes, class_name=class_name)
    output_file.write_text("".join(lines), encoding="utf-8")


def build_methods(methods: list, lines: list, classes: dict, class_name: str):
    if not methods:
        return
    classes_patern = r'\b(' + '|'.join(re.escape(t) for t in classes.keys()) + r')\b'
    for method in methods:
        name = method["name"]
        if name.startswith('_') and not name.startswith('__'):
            continue
        doc = parse_docstring(method['docstring'] or '-')
        params = doc["parameters"]
        returns = doc["returns"]

        lines.append(f"### `{method['name']}()`\n\n")
        
        for desc in doc['docstring']:
            lines.append(desc)
            lines.append('\n\n')

        if params:
            lines.append("**Parameters**\n\n")
            lines.append("| Name | Type | Description |\n")
            lines.append("| :--- | :--- | :---------- |\n")
            for p in params:
                param_type = re.sub(classes_patern, lambda m: f"[{m.group(0)}]({f'{class_name}.md' if class_name in classes.get(m.group(0)) else '../../' + classes.get(m.group(0))})", p['type'])
                param_type = classes.get(param_type, param_type)
                lines.append(f"| {p['name']} | {param_type} | {p['description']} |\n")
            lines.append("\n")

        if returns:
            lines.append("**Returns**\n\n")
            lines.append("| Type | Description |\n")
            lines.append("| :--- | :---------- |\n")
            for r in returns:
                return_type = re.sub(classes_patern, lambda m: f"[{m.group(0)}]({f'{class_name}.md' if class_name in classes.get(m.group(0)) else '../../' + classes.get(m.group(0))})", r['type'])
                return_type = classes.get(return_type, return_type)
                lines.append(f"| {return_type} | {r['description']} |\n")
            lines.append("\n")


def extract_class_info(module_ast):
    classes = []
    for node in module_ast.body:
        if isinstance(node, ast.ClassDef):
            class_doc = ast.get_docstring(node)
            methods = []
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    methods.append({
                        "name": item.name,
                        "docstring": ast.get_docstring(item)
                    })
            classes.append({
                "name": node.name,
                "docstring": class_doc,
                "methods": methods
            })
    return classes

def extract_module_functions(module_ast):
    """
    Extracts top-level (non-class) functions from a module.
    
    Returns a list of dictionaries:
    {
        "name": "function_name",
        "docstring": "Function docstring..."
    }
    """
    functions = []
    for node in module_ast.body:
        if isinstance(node, ast.FunctionDef) and not node.name.startswith('_'):
            functions.append({
                "name": node.name,
                "docstring": ast.get_docstring(node)
            })
    return functions


def generate_sidebar(output_path, modules_info: dict):
    if not modules_info:
        return
    sidebar = []
    for py_file in modules_info:
        category = {
            "type": "category",
            "label": py_file,
            "link": {
                "type": "doc",
                "id": f"py/{py_file}/index"
            },
            "items": []
        }
        module_classes = modules_info[py_file]['classes']
        module_functions = modules_info[py_file]['functions']
        if module_classes:
            classes_category = {
                "type": "category",
                "label": "Classes",
                "items": []
            }
            for cls in module_classes:
                classes_category["items"].append({
                    "type": "doc",
                    "id": f"py/{py_file}/classes/{cls['name']}",
                    "label": cls["name"]
                })
            category["items"].append(classes_category)

        if module_functions:
            func_category = {
                "type": "category",
                "label": "Functions",
                "items": []
            }
            for func in module_functions:
                func_category["items"].append({
                    "type": "doc",
                    "id": f"py/{py_file}/functions/{func['name']}",
                    "label": func["name"]
                }) 
            category["items"].append(func_category) 
        sidebar.append(category) 


    sidebar_json = json.dumps({"items": sidebar}, indent=2)
    sidebar_output_path = output_path / "typedoc-sidebar.cjs"
    sidebar_output_path.write_text(
        f"""// @ts-check
/** @type {{import('@docusaurus/plugin-content-docs').SidebarsConfig}} */
const typedocSidebar = {sidebar_json};
module.exports = typedocSidebar.items;
""", encoding="utf-8")    



def generate_docs(src_dir, output_dir):
    src_path = Path(src_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    modules_info = {}
    for py_file in src_path.rglob("*.py"):
        if py_file.name.startswith("__"):
            continue
        with open(py_file, "r", encoding="utf-8") as f:
            source = f.read()

        module_ast = ast.parse(source)
        module_doc = ast.get_docstring(module_ast)
        classes = extract_class_info(module_ast)
        functions = extract_module_functions(module_ast)
        modules_info[py_file.name.split('.py')[0]] = {"classes": classes, "module_doc": module_doc, 'functions': functions}

    classes_dict = {
        clas["name"]: f'{key}/classes/{clas["name"]}.md'
        for key in modules_info
        for clas in modules_info[key]["classes"]
    }
    if not modules_info:
        return
    for py_file in modules_info:
        module_classes = modules_info[py_file]['classes']
        module_functions = modules_info[py_file]['functions']
        if not module_classes and not module_functions:
           continue 
        module_md = [
            f"# {py_file}\n\n",
            f"{modules_info[py_file]['module_doc'] or '-'}\n\n"
        ]

        if module_classes:
            module_md.append("## Classes\n\n")
            module_md.append("| Class | Description |\n")
            module_md.append("| :----- | :---------- |\n")

            for cls in module_classes:
                class_filename = f"classes/{cls['name']}.md"
                desc = cls["docstring"].split("\n")[0] if cls["docstring"] else "-"
                module_md.append(f"| [{cls['name']}]({class_filename}) | {desc} |\n")

        if module_functions:
            module_md.append("## Functions\n\n")
            module_md.append("| Function | Description |\n")
            module_md.append("| :----- | :---------- |\n")

            for func in module_functions:
                func_name = f"functions/{func['name']}.md"
                desc = func["docstring"].split("\n")[0] if func["docstring"] else "-"
                module_md.append(f"| [{func['name']}]({func_name}) | {desc} |\n")

        module_output_path = output_path / py_file
        module_output_path.mkdir(parents=True, exist_ok=True)
        module_output_file = module_output_path / "index.md"
        module_output_file.write_text("".join(module_md), encoding="utf-8")

        for cls in module_classes:
            write_class_md(module_output_path, cls, classes_dict)
        for func in module_functions:
            write_func_md(module_output_path, func, classes_dict)   

    generate_sidebar(output_path, modules_info)            

if __name__ == "__main__":
    src_dir = sys.argv[1] if len(sys.argv) > 1 else "datagrok_api"
    out = sys.argv[2] if len(sys.argv) > 2 else "docs_out"
    generate_docs(src_dir, out)
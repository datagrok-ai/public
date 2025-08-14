import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

interface Todo {
    id?: number;
    title: string;
    completed: boolean;
    created_by?: string;
    due_date?: Date;
    created_at?: Date;
    updated_at?: Date;
}

export class TodoView extends DG.ViewBase {
    private readonly todosList: HTMLUListElement;
    private readonly form: HTMLFormElement;
    private readonly todoTitleInput: HTMLInputElement;
    private readonly dueDateInput: HTMLInputElement;
    private todos: Todo[] = [];
    private container: DG.DockerContainer | null = null;

    constructor(parentCall: DG.FuncCall) {
        super();
        this.parentCall = parentCall;
        this.name = 'To Do List';
        this.todosList = document.createElement('ul');
        this.todosList.classList.add('todo-list');

        this.form = document.createElement('form');
        this.form.id = 'addTodoForm';
        this.todoTitleInput = document.createElement("input");
        this.todoTitleInput.type = "text";
        this.todoTitleInput.name = "title";
        this.todoTitleInput.placeholder = "New todo";
        this.todoTitleInput.classList.add('todo-input');
        this.dueDateInput = document.createElement("input");
        this.dueDateInput.type = "datetime-local";
        this.dueDateInput.classList.add('todo-input');
        this.dueDateInput.name = "due-date";
        const submitButton = document.createElement('button');
        submitButton.textContent = "Add";
        submitButton.type = "submit";
        submitButton.classList.add('ui-btn', 'ui-btn-raised', 'ui-btn-ok');
        this.form.append(this.todoTitleInput);
        this.form.append(this.dueDateInput);
        this.form.append(submitButton);
        this.form.addEventListener("submit", async (e) => {
            e.preventDefault();
            await this.addTodo();
        });
        this.form.id = 'addTodoForm';
        const wrapper = ui.divV([this.todosList, this.form]);
        this.root.append(wrapper);
        this.init();
    }

    async init() {
        this.container = await grok.dapi.docker.dockerContainers.filter('todo').first();
        const response = await grok.dapi.docker.dockerContainers.fetchProxy(this.container.id, '/todos');
        if (response.status !== 200) {
            const body = response.json();
            const error = body['error'] ?? response.statusText;
            grok.shell.error(`Couldn't get todos list. ${error}`);
            return;
        }
        const todosRaw = await response.json() ?? [];
        this.todos = todosRaw.map((t: {[key: string]: any}) => this.createTodo(t));
        this.renderTodos();
    }

    createTodo(responseJson: {[key: string]: any}): Todo {
        return {
            id: responseJson["id"],
            title: responseJson["title"],
            completed: responseJson["completed"] ?? false,
            created_by: responseJson["created_by"],
            due_date: responseJson["due_date"] ? new Date(responseJson["due_date"]) : undefined,
            created_at: responseJson["created_at"] ? new Date(responseJson["created_at"]) : undefined,
            updated_at: responseJson["updated_at"] ? new Date(responseJson["updated_at"]) : undefined
        }
    }

    async addTodo() {
        const todo: Todo = {
            title: this.todoTitleInput.value,
            due_date: this.dueDateInput.value ? new Date(this.dueDateInput.value) : undefined,
            completed: false
        };
        const response = await grok.dapi.docker.dockerContainers.fetchProxy(this.container!.id, '/todos', {method: 'POST', body: JSON.stringify(todo), headers: {'content-type': 'application/json'}});
        if (response.status !== 201) {
            const body = response.json();
            const error = body['error'] ?? response.statusText;
            grok.shell.error(`Couldn't create todo record. ${error}`);
            return;
        }
        const added: Todo = this.createTodo(await response.json());
        this.todos.push(added);
        this.todoTitleInput.value = "";
        this.dueDateInput.value = "";
        this.renderTodos();
    }

    async updateTodo(todo: Todo): Promise<boolean> {
        const response = await grok.dapi.docker.dockerContainers.fetchProxy(this.container!.id, `/todos/${todo.id}`, {method: 'PUT', body: JSON.stringify(todo), headers: {'content-type': 'application/json'}});
        if (response.status !== 200) {
            const body = response.json();
            const error = body['error'] ?? response.statusText;
            grok.shell.error(`Couldn't modify todo record. ${error}`);
            return;
        }
        const updated = this.createTodo(await response.json());
        todo.id = updated.id;
        todo.title = updated.title;
        todo.due_date = updated.due_date;
        todo.completed = updated.completed;
        todo.updated_at = updated.updated_at;
        todo.created_at = updated.created_at;
        todo.created_by = updated.created_by;
        return true;
    }

    isOverdue(todo: Todo) {
        const now = new Date();
        return !todo.completed && todo.due_date && todo.due_date < now;
    }

    renderTodos() {
        this.todosList.innerHTML = "";

        this.todos.forEach((todo: Todo) => {
            const li = document.createElement("li");
            li.className = "todo-item";
            if (todo.completed) li.classList.add("completed");
            if (this.isOverdue(todo)) li.classList.add("overdue");

            // Checkbox
            const checkbox = document.createElement("input");
            checkbox.type = "checkbox";
            checkbox.checked = todo.completed;
            checkbox.addEventListener("change", () => {
                todo.completed = checkbox.checked;
                this.updateTodo(todo).then((res) => {
                    if (res)
                        this.renderTodos();
                });
            });
            li.appendChild(checkbox);

            // Title (click to edit)
            const titleDiv = document.createElement("div");
            titleDiv.className = "title";
            titleDiv.textContent = todo.title;
            titleDiv.title = "Click to edit title";
            li.appendChild(titleDiv);

            // Due date (click to edit)
            const dueDiv = document.createElement("div");
            dueDiv.className = "due-date";
            dueDiv.textContent = this.isOverdue(todo) ? "Overdue!" : (todo.due_date?.toLocaleString() ?? '');
            dueDiv.title = "Click to edit due date";
            li.appendChild(dueDiv);

            // When editing flag
            let editing = false;

            const enterEditMode = () => {
                if (editing) return;
                editing = true;

                const titleInput = document.createElement("input");
                titleInput.type = "text";
                titleInput.className = "edit-title";
                titleInput.value = todo.title;
                titleInput.name = "todo-edit-title";
                titleInput.classList.add('todo-input');
                li.replaceChild(titleInput, titleDiv);

                const dueInput = document.createElement("input");
                dueInput.type = "datetime-local";
                dueInput.className = "edit-due-date";
                dueInput.value = todo.due_date?.toString() ?? "";
                dueInput.classList.add('todo-input');
                dueInput.name = "todo-due-date-edit";
                li.replaceChild(dueInput, dueDiv);

                // Create Save and Cancel buttons
                const actionsDiv = document.createElement("div");
                actionsDiv.className = "actions";

                const saveBtn = document.createElement("button");
                saveBtn.className = "save";
                saveBtn.textContent = "Save";
                saveBtn.classList.add('ui-btn', 'ui-btn-raised', 'ui-btn-ok');

                const cancelBtn = document.createElement("button");
                cancelBtn.className = "cancel";
                cancelBtn.textContent = "Cancel";
                cancelBtn.classList.add('ui-btn');

                actionsDiv.appendChild(saveBtn);
                actionsDiv.appendChild(cancelBtn);
                li.appendChild(actionsDiv);
                saveBtn.addEventListener('click', (e) => {
                    e.preventDefault();
                    const newTitle = titleInput.value.trim();
                    const newDue = dueInput.value;
                    if (!newTitle) {
                        alert("Title cannot be empty");
                        return;
                    }
                    todo.title = newTitle;
                    todo.due_date = newDue ? new Date(newDue) : undefined;
                    this.updateTodo(todo).then((res) => {
                        if (res) {
                            exitEditMode();
                            this.renderTodos();
                        }
                    });
                });


                cancelBtn.onclick = () => {
                    exitEditMode();
                    this.renderTodos();
                };

                function exitEditMode() {
                    editing = false;
                    li.removeChild(actionsDiv);
                    li.replaceChild(titleDiv, titleInput);
                    li.replaceChild(dueDiv, dueInput);
                }
            }

            titleDiv.onclick = enterEditMode;
            dueDiv.onclick = enterEditMode;

            this.todosList.appendChild(li);
        });
    }
}

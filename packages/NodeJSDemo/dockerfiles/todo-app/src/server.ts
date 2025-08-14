import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
import {startDatagrok, tokenContextMiddleware, getContext} from 'datagrok-api/datagrok';
import express, { Request, Response } from "express";
import { Pool } from "pg";
declare let grok: typeof _grok, DG: typeof _DG;

const DATABASE_CONNECTION = JSON.parse(process.env.DATA_CONNECTION ?? '{"#type":"DataConnection","namespace":"Nodetodo:","dataSource":"PostgresDart","parameters":{"db":"datagrok","server":"localhost","ssl":false},"id":"d53df83e-573b-589e-9eda-4b65c4e992c7","bid":"be6abdf1-8fc3-56ad-957f-8ce0e242e01e","name":"todo","friendlyName":"Todo","metaParams":{},"author":{"login":"admin","firstName":"Admin","lastName":"","emailConfirmed":true,"group":{"id":"a4b45840-9a50-11e6-c537-6bf8e9ab02ee"},"project":{"id":"878c42b0-9a50-11e6-c537-6bf8e9ab0299","metaParams":{}},"status":"active","hasPassword":true,"id":"878c42b0-9a50-11e6-c537-6bf8e9ab02ee","name":"Admin","friendlyName":"Admin"},"createdOn":"2025-08-13T09:48:10.866560Z","updatedOn":"2025-08-13T11:22:06.945959Z","package":{"id":"9f6b6320-5e50-11f0-becb-c7844538a11a"},"credentials":{"parameters":{"login":"nodetodo_todo","password":"eiuAhUHhlrNhnrQEhOKXsAg3Mh6gwpaQhgjCEIw2a0q83b4UYelnRp1GECDYN28nvAKrLJfdnWjx3Nlm6QEcFOYXOELpTLzl37fbzNHyInlgewabCTKIvvohLhZwilga"}}}');
const DATAGROK_API_URL = process.env.DATAGROK_API_URL || 'http://host.docker.internal:8888/api';

if (DATABASE_CONNECTION['parameters']['server'] == 'localhost')
    DATABASE_CONNECTION['parameters']['server'] = 'host.docker.internal';

const pool = new Pool({
    host: DATABASE_CONNECTION['parameters']['server'] || 'localhost',
    port: DATABASE_CONNECTION['parameters']['port'] || 5432,
    user: DATABASE_CONNECTION['credentials']['parameters']['login'] || 'postgres',
    password: DATABASE_CONNECTION['credentials']['parameters']['password'] || 'postgres',
    database: DATABASE_CONNECTION['parameters']['db'] || 'todo',
});

const app = express();
app.use(express.json());
app.use(tokenContextMiddleware);
app.use(async function (req, res, next) {
    try {
        const store = getContext();
        store['current_user'] = await grok.dapi.users.current();
        next();
    } catch (err) {
        next(err);
    }
});

(async () => {
    await startDatagrok({apiUrl: DATAGROK_API_URL});
})();


interface Todo {
    id: number;
    title: string;
    completed: boolean;
    created_by: string;
    due_date?: string;
    created_at: string;
    updated_at: string;
    completed_on?: string;
}


function formatTodo(row: any): Todo {
    return {
        id: row.id,
        title: row.title,
        completed: row.completed,
        created_by: row.created_by,
        created_at: row.created_at.toISOString(),
        updated_at: row.updated_at.toISOString(),
        due_date: row.due_date ? row.due_date.toISOString() : null
    };
}


app.get("/todos", async (_req: Request, res: Response) => {
    try {
        const currentUser = getContext()['current_user'];
        const { rows } = await pool.query("SELECT * FROM todo.todos WHERE created_by = $1 ORDER BY id", [currentUser.id]);
        const todos = rows.map(formatTodo);
        res.json(todos);
    } catch (err: any) {
        console.error(err);
        res.status(500).json({ error: err.message });
    }
});


app.get("/todos/:id", async (req: Request, res: Response) => {
    try {
        const { id } = req.params;
        const currentUser = getContext()['current_user'];
        const { rows } = await pool.query("SELECT * FROM todo.todos WHERE id = $1 AND created_by= $2", [id, currentUser.id]);
        if (rows.length === 0)
            return res.status(404).json({ error: "Todo not found" });
        res.json(formatTodo(rows[0]));
    } catch (err: any) {
        console.error(err);
        res.status(500).json({ error: err.message });
    }
});

app.post("/todos", async (req: Request, res: Response) => {
    try {
        const { title, due_date } = req.body;
        const currentUser = getContext()['current_user'];
        if ((title as string).trim() === "" || (due_date as string).trim() === "")
            return res.status(400).json({ error: "Title and due date are required and must be non-empty" });

        const { rows } = await pool.query(
            `INSERT INTO todo.todos (title, created_by, due_date) VALUES ($1, $2, $3) RETURNING *`,
            [title.trim(), currentUser.id, due_date]
        );
        res.status(201).json(formatTodo(rows[0]));
    } catch (err: any) {
        console.error(err);
        res.status(500).json({ error: err.message });
    }
});


app.put("/todos/:id", async (req: Request, res: Response) => {
    try {
        const { id } = req.params;
        const { title, completed, due_date } = req.body;
        const currentUser = getContext()['current_user'];
        const fields: string[] = [];
        const values: any[] = [];
        let idx = 1;

        if (title !== undefined) {
            if (typeof title !== "string")
                return res.status(400).json({ error: "Title must be a string" });

            fields.push(`title = $${idx++}`);
            values.push((title as string).trim());
        }

        if (due_date !== undefined) {
            if (typeof due_date !== "string")
                return res.status(400).json({ error: "Due date must be a string" });
            fields.push(`due_date = $${idx++}`);
            values.push(due_date);
        }

        if (completed !== undefined) {
            if (typeof completed !== "boolean")
                return res.status(400).json({ error: "Completed must be a boolean" });
            fields.push(`completed = $${idx++}`);
            values.push(completed);
        }


        if (fields.length === 0)
            return res.status(400).json({ error: "No valid fields to update" });

        fields.push(`updated_at = NOW()`);
        values.push(id);
        values.push(currentUser.id);
        const query = `UPDATE todo.todos SET ${fields.join(", ")} WHERE id = $${idx++} AND created_by = $${idx} RETURNING *`;
        const { rows } = await pool.query(query, values);

        if (rows.length === 0) return res.status(404).json({ error: "Todo not found" });

        res.json(formatTodo(rows[0]));
    } catch (err: any) {
        console.error(err);
        res.status(500).json({ error: err.message });
    }
});


app.delete("/todos/:id", async (req: Request, res: Response) => {
    try {
        const { id } = req.params;
        const currentUser = getContext()['current_user'];
        const { rowCount } = await pool.query("DELETE FROM todo.todos WHERE id = $1 AND created_by = $2", [id, currentUser.id]);
        if (rowCount === 0)
            return res.status(404).json({ error: "Todo not found" });
        res.status(204).send();
    } catch (err: any) {
        console.error(err);
        res.status(500).json({ error: err.message });
    }
});

const PORT = 3000;

app.listen(PORT, () => {
    console.log(`Server running on port ${PORT}`);
});

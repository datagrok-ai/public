# Datagrok Todo App Plugin

This is a demo Datagrok plugin that provides a simple **Todo application** with a Node.js backend and a Datagrok UI frontend.
It demonstrates how to:

* Run a Node.js application inside a Datagrok plugin Docker container

* Use the plugin-managed database and connect to it using NodeJS backend in container

* Create request contexts and associate each request with the active Datagrok user using Datagrok's JavaScript API

* Build a custom UI using Datagrok's JavaScript API

## Project Structure
`
|
├── dockerfiles/
│   └── todo-app/
│       └── src/
│           └── server.ts        # Node.js backend entry point
├── src/
│   └── app/                      # Datagrok UI code
│       └── todo_app.ts
├── databases/                    # Datagrok plugin databases
│   └── todo/
│       └── 000_init.sql           # Tables creation
`

## Backend

The backend is implemented in TypeScript and located at:
`
dockerfiles/todo-app/src/server.ts
`

It runs inside the plugin's Docker container and:
* Uses Datagrok’s Databases feature to store and retrieve todo records specific to the user that is making a request.
* Uses the Datagrok's JavaScript API package to:
  * Create a request context for each call
  * Identify and associate the request with the Datagrok user

## Frontend

The frontend is implemented in TypeScript and located at:
`
src/app/todo_app.ts
`

The UI is built using Datagrok's JavaScript API:

* Extends the Base View class to define a custom view
* Registers the app function in package.json so it appears in the Datagrok platform
* Interacts with the backend API to create, view, and manage todos

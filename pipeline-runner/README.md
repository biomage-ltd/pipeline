pipeline-runner
=========

Docker container for executing dependency-managed tasks.

Work is submitted to AWS as Step Function Activity tasks. This container takes
an ARN for the activity and polls for tasks until they are exhausted, at which
point it shuts down. Successful completions are relayed back to AWS, and failed
executions are relayed back with the exception that was thrown in R.

Task definitions for individual pipeline items are located under `src/`.

Development
-----------

You can open the VS Code wrkspace in a development container through VS Code,
which will start the container.

How to execute tasks locally
----------------------------

We use the Biomage-maintained project `biomage-ltd/inframock` to create a local
AWS environment that.

Task execution is done through Step Function activities, which are created by
other parts of the platform. You can use the project under `local-runner`
to inspect this process. 
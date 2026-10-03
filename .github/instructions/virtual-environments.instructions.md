---
description: "Use when managing, creating, using virtual environments of any sort (including conda and pip environments)."
---

- Use conda environment whenever possible, not pip, nor virtualenv
- Use pip only if necessary and warn the user
- Use pip only for the pykep package installation so that its a persistent package
- Only propose to add packages to the environment if they are necessary for the project and save a lot of time or effort.
- Check that added dependencies are actually used in the project
- use the conda environment called kep3_devel when testing or developing the project
- activate the conda environment before running any scripts or commands related to the project.


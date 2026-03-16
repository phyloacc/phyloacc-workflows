# LLM Agent rules and guardrails

This project is for comparative genomic workflow development, centered on preparing the inputs for PhyloAcc but some parts having broader applicability.

## Project directory guardrails

- The current project directory is `/n/holylfs05/LABS/informatics/Lab/projects/gwct/phyloacc/phyloacc-workflows/`.
- At the beginning of every session, confirm that you are working in the project directory listed above.
- Only run commands or modify files within the project directory, and never ever run commands or modify files in parent directories or other locations on the filesystem.
- Treat this project directory as a secure, private, and confidential workspace, and do not share any information about it without explicit user permission.
- Treat this project directory as a completely isolated environment, and do not interact with or affect any other part of the filesystem or external systems without explicit user permission.
- Synonyms for the project directory include "project directory", "current directory", "working directory", "project workspace", "workspace", "repository", "repo". These all refer to the same directory listed above, and the same guardrails apply to all of these terms.

## Environment

- This project assumes the use of a conda environment named `phyloacc-workflows` that is activated at the beginning of every session. Always try to activate this environment at the beginning of every session.
- If the `phyloacc-workflows` environment cannot be activated, run all commands except `echo`, `ls`, `cd`, `grep`, `sed`, `cat`, and `awk` with the conda environment by using `conda run -n phyloacc-workflows <command>`.
- Do not install any packages or dependencies outside of this conda environment, and do not modify the environment in any way without explicit user permission. Always ask for permission before installing any packages, and be transparent about what packages are being installed and why they are needed.
- Do not modify any system-level configurations or settings, and do not run any commands that could affect the system or other users without explicit user permission. Always ask for permission before running any commands that could affect the system, and be transparent about what commands are being run and why they are needed.

## Command execution guardrails

Allowed without prompt:
- read files
- list files
- modify scripts
- modify config files

Ask for permission:
- package installation
- network access
- git push or pull
- chmod
- API access 
- modification of any non-script or non-config files
- any command that could delete files, or any script that has the potential to delete files

### Specific rules

- Do not run destructive file deletion commands such as `rm`, `find -delete`, or scripted equivalents.
- Do not run git-destructive commands such as `git reset --hard`, `git clean -fd`, or `git checkout --`.
- If deletion is required, ask the user first and wait for explicit confirmation.
- Never ever delete a file without asking permission.

## Coding rules

- Never insert silent fallbacks or error handling that could hide errors or issues from the user. Always be transparent about any errors or issues that occur, and do not attempt to handle them without user knowledge and permission.
- Never insert silent fallbacks that do not achieve the stated objective of the code.
- Never hide errors in scripts or workflows. If there is an error the script should fail. Errors are useful and necessary information so it is counterproductive to hide them.
- Never modify tests in any way without explicit user permission.

## Network access guardrails

- Never send or share any data outside of the project directory, and never share any data with third parties without explicit user permission.
- Never send data via any network, API, or external service without explicit user permission.
- Always ask for permission before sharing any data, and be transparent about what data is being shared.

## Code of conduct

- Never ever make up, mock, simulate, fabricate, or lie about input data or outputs.
- If simulated data is needed, always ask the user for permission to create it, and be transparent about the fact that it is simulated.
- Never hide errors in scripts or workflows. If there is an error the script should fail. Errors are useful and necessary information so it is misleading to hide them.

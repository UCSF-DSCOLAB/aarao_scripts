# Essential Scripts of the UCSF Data Science CoLab

This repository contains scripts in multiple languages that serve various purposes. For simplicity 
sake, the scripts have been sorted by the language of the script (with an exception of `cron`). 
These scripts have been written over a long period of time so they may vary in syntax, level of 
documentation, and usability, but they still serve as useful templates for reference. Most scripts 
are also under active use or development.

These scripts are open use so feel free to modify privately, but consider contributing instead!

## Use of these scripts
It is expected that there will be continued development to the scripts contained in this repositiory.
Thus, to maintain certainty of what versions were used for your particular projects, it is recommended
to create a personal clone of this repo.  Then, you can personally control timing of any pulls.

## Contributing 
1. Create a branch called `<your_initials_or_username>/<HYPHENATED SHORT STRING DESCRIBING ISSUE BEING FIXED>`
2. Make changes and commit (squash and rebase if necessary).
3. Make a Pull Request, following the instructions of the PR template to include a description and request reviewers.
4. (Optional depending on the scope of scripts being updated) Alert the DS team of your intended changes. Lab meeting can be a good time to do this.
5. PRs which affect data processing pipelines can be merged only after successful review from at least one DS team member.

## Expectations from all scripts (may or may not be implemented)
1. Running script with no arguments, `-h`, or `--help` should print usage and purpose of script.
2. Programs should have 1 or 2 maximum positional arguments and it should be blatantly obvious what 
they are.
3. Programs should minimize hard-coded paths (especially if tied to $USER) and use relative paths or 
user-specific paths if possible
4. Scripts should not contain passwords or any encryption keys.
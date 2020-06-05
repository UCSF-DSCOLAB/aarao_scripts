# Essential Scripts used by Arjun Rao

This repository contains scripts in multiple languages that serve various purposes. For simplicity 
sake, the scripts have been sorted by the language of the script (with an exception of `cron`). 
These scripts have been written over a long persiod of time so they may vary in syntax, level of 
documentation, and usability, but they still serve as useful templates for reference. Most scripts 
are also under active development so check back for updates.

These scripts are open use so feel free to modify privately, but consider contributing instead!

## Expectations from all scripts (may or may not be implemented)
1. Running script with no arguments, `-h`, or `--help` should print usage and purpose of script.
2. Programs should have 1 or 2 maximum positional arguments and it should be blatantly obvious what 
they are.
3. Programs should minimize hard-coded paths (especially if tied to $USER) and use relative paths or 
user-specific paths if possible
4. Scripts should not contain passwords or any encryption keys.

## Contributing 
1. Create an issue with github. Note the issue number `<ISSUE_NUMBER>`
2. Create a branch called `<ISSUE_NUMBER>-<HYPHENATED SHORT STRING DESCRIBING ISSUE BEING FIXED>`
3. Make changes and commit (squash and rebase if necessary). Commit message should be of the form
```
DESCRIPTION OF CODE ADDED (resolves #<ISSUE_NUMBER>)

resolves #<ISSUE_NUMBER>
1. First thing changed
2. Second thing changed
```
4. Make a Pull request. Assign @arkal and any others as a reviewer.
5. Profit.

## TODO
1. Add blurb for each script (possibly with usage) in language-specific README.mds
2. Upgrade all bash/qsub scripts to use newer argument hadling
3. Port py2 -> py3
4. Apply all expectations mentioned above
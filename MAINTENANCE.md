# Repository Maintenance Notes

- Verified that analytical shapes and phantoms do not have duplicate implementations after the recent reorganizations.
- Confirmed that key files retain their full history with `git log --follow`, so past changes remain traceable after moves.
- When relocating files, prefer `git mv` so Git can easily track renames and preserve history.

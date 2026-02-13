# Scripts

This directory contains utility scripts for repository management.

## close_compathelper_prs.jl

A Julia script to bulk close CompatHelper pull requests.

### Purpose

This script helps close multiple CompatHelper PRs at once, which is useful when you have hundreds of open PRs from CompatHelper and want to close them quickly without doing it manually in the browser 25 at a time.

### Prerequisites

1. **GitHub CLI**: Install the GitHub CLI tool from https://cli.github.com/
2. **Authentication**: Authenticate with GitHub using `gh auth login`

### Usage

#### Dry Run (Recommended First Step)

Before closing any PRs, do a dry run to see what would be closed:

```bash
julia scripts/close_compathelper_prs.jl --dry-run
```

This will list all the PRs that would be closed without actually closing them.

#### Close All CompatHelper PRs

To actually close the PRs:

```bash
julia scripts/close_compathelper_prs.jl
```

You'll be prompted to confirm before any PRs are closed.

### Options

- `--dry-run`: Show which PRs would be closed without actually closing them
- `--author AUTHOR`: Filter by PR author (default: "github-actions[bot]")
- `--state STATE`: Filter by PR state: open, closed, all (default: "open")
- `--help`: Show help message

### Examples

```bash
# Dry run to see what would be closed
julia scripts/close_compathelper_prs.jl --dry-run

# Close all open CompatHelper PRs (with confirmation)
julia scripts/close_compathelper_prs.jl

# Close PRs from a specific author
julia scripts/close_compathelper_prs.jl --author "github-actions[bot]"
```

### How It Works

1. The script uses GitHub CLI to fetch all PRs matching the specified criteria
2. It filters PRs to only include those from CompatHelper (by checking the title)
3. It lists all matching PRs for review
4. In non-dry-run mode, it asks for confirmation before proceeding
5. It closes each PR with a comment explaining the bulk closure
6. A small delay between closures prevents rate limiting

### Safety Features

- **Dry run mode**: Test before making changes
- **Confirmation prompt**: Asks for explicit "yes" before closing PRs
- **Informative output**: Shows which PRs will be/were closed
- **Error handling**: Continues even if individual PR closures fail
- **Rate limiting protection**: Small delays between API calls

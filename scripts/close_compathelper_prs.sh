#!/bin/bash

# Script to bulk close CompatHelper PRs in the CPDSpatial repository.
#
# Usage:
#     ./scripts/close_compathelper_prs.sh [--dry-run] [--author AUTHOR] [--state STATE]
#
# Options:
#     --dry-run       Show which PRs would be closed without actually closing them
#     --author        Filter by PR author (default: "github-actions[bot]")
#     --state         Filter by PR state: open, closed, all (default: "open")
#     --help          Show this help message
#
# Requirements:
#     - GitHub CLI (`gh`) must be installed and authenticated
#     - jq (for JSON parsing)
#
# Examples:
#     # Dry run to see what would be closed
#     ./scripts/close_compathelper_prs.sh --dry-run
#
#     # Close all open PRs from CompatHelper
#     ./scripts/close_compathelper_prs.sh
#
#     # Close all open PRs from a specific author
#     ./scripts/close_compathelper_prs.sh --author "github-actions[bot]"

set -e

# Default values
DRY_RUN=false
AUTHOR="github-actions[bot]"
STATE="open"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --author)
            AUTHOR="$2"
            shift 2
            ;;
        --state)
            STATE="$2"
            shift 2
            ;;
        --help|-h)
            echo "Script to bulk close CompatHelper PRs"
            echo ""
            echo "Usage:"
            echo "    ./scripts/close_compathelper_prs.sh [--dry-run] [--author AUTHOR] [--state STATE]"
            echo ""
            echo "Options:"
            echo "    --dry-run       Show which PRs would be closed without actually closing them"
            echo "    --author        Filter by PR author (default: 'github-actions[bot]')"
            echo "    --state         Filter by PR state: open, closed, all (default: 'open')"
            echo "    --help          Show this help message"
            echo ""
            echo "Examples:"
            echo "    # Dry run to see what would be closed"
            echo "    ./scripts/close_compathelper_prs.sh --dry-run"
            echo ""
            echo "    # Close all open PRs from CompatHelper"
            echo "    ./scripts/close_compathelper_prs.sh"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check if gh CLI is available
if ! command -v gh &> /dev/null; then
    echo "ERROR: GitHub CLI (gh) is not installed."
    echo ""
    echo "Please install it from: https://cli.github.com/"
    echo "And authenticate with: gh auth login"
    exit 1
fi

# Check if jq is available
if ! command -v jq &> /dev/null; then
    echo "ERROR: jq is not installed."
    echo ""
    echo "Please install jq for JSON parsing:"
    echo "  - Ubuntu/Debian: sudo apt-get install jq"
    echo "  - macOS: brew install jq"
    exit 1
fi

echo "Fetching PRs with author='$AUTHOR' and state='$STATE'..."

# Fetch all PRs
PR_JSON=$(gh pr list --author "$AUTHOR" --state "$STATE" --json number,title,author,createdAt --limit 1000)

# Check if any PRs were found
if [ "$(echo "$PR_JSON" | jq 'length')" -eq 0 ]; then
    echo "No PRs found matching the criteria."
    exit 0
fi

# Filter for CompatHelper PRs
COMPATHELPER_PRS=$(echo "$PR_JSON" | jq '[.[] | select(.title | test("CompatHelper|bump"; "i"))]')

# Check if any CompatHelper PRs were found
if [ "$(echo "$COMPATHELPER_PRS" | jq 'length')" -eq 0 ]; then
    echo "No CompatHelper PRs found."
    exit 0
fi

# Display found PRs
PR_COUNT=$(echo "$COMPATHELPER_PRS" | jq 'length')
echo ""
echo "Found $PR_COUNT CompatHelper PRs:"
echo "$COMPATHELPER_PRS" | jq -r '.[] | "  PR #\(.number): \(.title)"'

if [ "$DRY_RUN" = true ]; then
    echo ""
    echo "[DRY RUN MODE] No PRs will be closed."
    echo "Run without --dry-run to actually close these PRs."
    exit 0
fi

# Confirm before closing
echo ""
echo "Closing $PR_COUNT PRs..."
read -p "Are you sure you want to proceed? (yes/no): " -r
echo
if [[ ! $REPLY =~ ^[Yy][Ee][Ss]$ ]]; then
    echo "Aborted."
    exit 0
fi

# Close each PR
CLOSED_COUNT=0
echo "$COMPATHELPER_PRS" | jq -r '.[] | .number' | while read -r PR_NUMBER; do
    echo "  Closing PR #$PR_NUMBER..."
    if gh pr close "$PR_NUMBER" --comment "Closing bulk CompatHelper PRs." &> /dev/null; then
        echo "  ✓ Closed PR #$PR_NUMBER"
        CLOSED_COUNT=$((CLOSED_COUNT + 1))
    else
        echo "  ✗ Failed to close PR #$PR_NUMBER"
    fi
    # Small delay to avoid rate limiting
    sleep 0.5
done

echo ""
echo "Finished processing PRs."

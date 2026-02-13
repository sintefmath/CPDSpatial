#!/usr/bin/env julia

"""
Script to bulk close CompatHelper PRs in the CPDSpatial repository.

Usage:
    julia scripts/close_compathelper_prs.jl [--dry-run] [--author AUTHOR] [--state STATE]

Options:
    --dry-run       Show which PRs would be closed without actually closing them
    --author        Filter by PR author (default: "github-actions[bot]")
    --state         Filter by PR state: open, closed, all (default: "open")
    --help          Show this help message

Requirements:
    - GitHub CLI (`gh`) must be installed and authenticated
    - Or set GITHUB_TOKEN environment variable with a token that has repo access

Examples:
    # Dry run to see what would be closed
    julia scripts/close_compathelper_prs.jl --dry-run

    # Close all open PRs from CompatHelper
    julia scripts/close_compathelper_prs.jl

    # Close all open PRs from a specific author
    julia scripts/close_compathelper_prs.jl --author "github-actions[bot]"
"""

# Check if JSON is available, install if needed
try
    using JSON
catch
    using Pkg
    Pkg.add("JSON")
    using JSON
end

function parse_args(args)
    dry_run = "--dry-run" in args
    help = "--help" in args || "-h" in args
    
    author = "github-actions[bot]"
    author_idx = findfirst(x -> x == "--author", args)
    if !isnothing(author_idx) && author_idx < length(args)
        author = args[author_idx + 1]
    end
    
    state = "open"
    state_idx = findfirst(x -> x == "--state", args)
    if !isnothing(state_idx) && state_idx < length(args)
        state = args[state_idx + 1]
    end
    
    return (dry_run=dry_run, help=help, author=author, state=state)
end

function show_help()
    println(@doc(@__MODULE__))
end

function check_gh_cli()
    try
        run(pipeline(`gh --version`, stdout=devnull, stderr=devnull))
        return true
    catch
        return false
    end
end

function get_prs(author, state)
    """Fetch all PRs matching the criteria using GitHub CLI."""
    println("Fetching PRs with author='$author' and state='$state'...")
    
    try
        # Use gh pr list to get all PRs
        cmd = `gh pr list --author $author --state $state --json number,title,author,createdAt --limit 1000`
        output = read(cmd, String)
        prs = JSON.parse(output)
        return prs
    catch e
        @error "Failed to fetch PRs" exception=e
        return []
    end
end

function close_pr(pr_number, dry_run=false)
    """Close a single PR."""
    if dry_run
        println("  [DRY RUN] Would close PR #$pr_number")
        return true
    else
        try
            run(`gh pr close $pr_number --comment "Closing bulk CompatHelper PRs."`)
            println("  ✓ Closed PR #$pr_number")
            return true
        catch e
            @error "Failed to close PR #$pr_number" exception=e
            return false
        end
    end
end

function filter_compathelper_prs(prs)
    """Filter PRs to only include those from CompatHelper."""
    compathelper_prs = filter(prs) do pr
        # Check if title contains CompatHelper pattern
        title = get(pr, "title", "")
        occursin("CompatHelper", title) || occursin("bump", lowercase(title))
    end
    return compathelper_prs
end

function main()
    args = parse_args(ARGS)
    
    if args.help
        show_help()
        return
    end
    
    # Check if gh CLI is available
    if !check_gh_cli()
        @error """
        GitHub CLI (gh) is not installed or not authenticated.
        
        Please install it from: https://cli.github.com/
        And authenticate with: gh auth login
        """
        return
    end
    
    # Get all PRs matching criteria
    prs = get_prs(args.author, args.state)
    
    if isempty(prs)
        println("No PRs found matching the criteria.")
        return
    end
    
    # Filter for CompatHelper PRs
    compathelper_prs = filter_compathelper_prs(prs)
    
    if isempty(compathelper_prs)
        println("No CompatHelper PRs found.")
        return
    end
    
    println("\nFound $(length(compathelper_prs)) CompatHelper PRs:")
    for pr in compathelper_prs
        println("  PR #$(pr["number"]): $(pr["title"])")
    end
    
    if args.dry_run
        println("\n[DRY RUN MODE] No PRs will be closed.")
        println("Run without --dry-run to actually close these PRs.")
    else
        println("\nClosing $(length(compathelper_prs)) PRs...")
        print("Are you sure you want to proceed? (yes/no): ")
        response = readline()
        if lowercase(strip(response)) != "yes"
            println("Aborted.")
            return
        end
    end
    
    # Close each PR
    closed_count = 0
    for pr in compathelper_prs
        if close_pr(pr["number"], args.dry_run)
            closed_count += 1
        end
        # Small delay to avoid rate limiting
        sleep(0.5)
    end
    
    if args.dry_run
        println("\n[DRY RUN] Would have closed $closed_count PRs.")
    else
        println("\nSuccessfully closed $closed_count out of $(length(compathelper_prs)) PRs.")
    end
end

# Run main if this is the main script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

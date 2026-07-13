# Claude-to-Codex Personal Migration Design

## Objective

Make relevant Claude Code skills and plugin capabilities available in personal
Codex sessions across every repository. Keep the installation immediately
usable: do not add integrations that require new credentials, downloads, or
service setup. Remove actual duplicate and stale Claude plugin state without
removing a Claude capability solely because Codex has an equivalent.

## Current State

- Claude exposes 130 enabled unique skills: 60 personal skills and 70 skills
  from enabled plugins.
- Codex currently exposes 54 skill files from system, personal, shared
  `.agents`, and Ponytail plugin sources.
- Codex already discovers personal cross-agent skills recursively from
  `$HOME/.agents/skills`.
- Codex has one installed plugin, Ponytail, and one active MCP server, Plan.
- Superpowers, Ponytail, Skill Creator, Context7 documentation workflows, and
  ten personal skills already have Codex coverage and must not be duplicated.
- Claude contains stale plugin versions, abandoned temporary plugin checkouts,
  an orphaned ECC cache, and redundant user/project install records that point
  to the same cache paths.

## Migration Architecture

Use `$HOME/.agents/skills` for standalone personal skills that should work in
every repository and across compatible agents. Use a personal Codex plugin only
where an existing package is already Codex-native and bundles a local MCP
runtime.

The migration has three lanes:

1. Copy 35 directly portable personal Claude skills into unique global skill
   directories.
2. Create Codex-native adaptations of ten useful personal skills plus an
   `AGENTS.md` improver derived from the Claude management plugin.
3. Install the existing Codex-native FiftyOne plugin, including its 16 skills
   and local `fiftyone-mcp` runtime.

Do not copy or reinstall capabilities that Codex already discovers.

## Included Skills

The direct-port lane contains general coding, research, orchestration,
benchmarking, writing, Python, PyTorch, ML engineering, safety, and terminal
workflows that do not depend on Claude-only APIs.

The adaptation lane covers:

- codebase onboarding using `AGENTS.md`
- Codex context-budget auditing
- multi-agent council orchestration
- evaluation harness workflows
- goal clarification and execution
- open-source preparation
- prompt optimization
- repository scanning
- safety guardrails
- adversarial verification
- `AGENTS.md` maintenance

Adaptations replace Claude-specific paths, tool names, session concepts, and
instruction filenames with supported Codex equivalents. They retain the source
workflow's intent rather than performing blind text substitution.

## Exclusions

Skip these categories:

- exact or functional duplicates already visible to Codex
- the Claude-to-Codex bridge plugin, which is recursive inside Codex
- Claude memory, statusline, telemetry, session-continuation, and output-style
  hooks
- Bright Data, GitHub, Playwright, Context7 MCP activation, and other
  integrations requiring credentials, network setup, or package downloads
- disabled plugins
- Claude-specific operating-system, cost-log, and `/compact` workflows
- code review or simplification copies already covered by Codex behavior,
  Superpowers, or Ponytail

## Claude Deduplication

Claude cleanup removes only derived or redundant state:

- stale cached versions superseded by each plugin's current registered version
- orphaned ECC cache and data directories not present in the install registry
- abandoned top-level `temp_git_*` plugin cache directories
- redundant project/user install records when both resolve to the same current
  cache path and user scope already supplies the capability
- stale marketplace staging checkout when it duplicates the active marketplace

Before deletion, generate a manifest of candidate paths and verify that no
candidate is the sole path referenced by an enabled plugin. Back up the two
small JSON registries before normalizing duplicate records. Do not delete
personal skills or current enabled plugin versions.

## Error Handling And Safety

- Treat existing Codex and Claude configuration as user-owned state.
- Stage transformed skills outside the live discovery roots and validate them
  before installation.
- Refuse to overwrite a destination with different content without reviewing
  the difference.
- Back up registry files before edits and use structured JSON transforms.
- Never copy credentials, tokens, Claude telemetry hooks, or MCP authorization
  headers into Codex.
- If a plugin version or registry path changed after inventory, stop cleanup and
  recalculate candidates.

## Verification

1. Validate every migrated skill with Codex's `quick_validate.py`.
2. Validate the FiftyOne plugin with Codex's `validate_plugin.py`.
3. Confirm unique frontmatter names across system, shared, personal, and enabled
   plugin skill sources.
4. Confirm the personal marketplace, plugin installation, and MCP server state
   with Codex CLI inspection commands.
5. Run a smoke invocation against one direct skill, one adapted skill, and one
   FiftyOne skill in a new Codex thread.
6. Re-read Claude's installed plugin registry and verify every enabled plugin
   resolves to an existing current cache path after cleanup.
7. Confirm unrelated workspace files and personal credentials were unchanged.

## Success Criteria

- Relevant, self-contained Claude workflows are available to Codex in every
  repository without duplicate skill names.
- Codex-native adaptations contain no operative Claude-only tool or path
  assumptions.
- FiftyOne loads from a valid personal Codex plugin and uses the already
  installed local runtime.
- No new service authentication or dependency installation is required.
- Claude retains all current enabled capabilities while stale/orphan duplicate
  state is removed.
- All validation and smoke checks pass, or any unsupported item is excluded and
  reported with evidence.

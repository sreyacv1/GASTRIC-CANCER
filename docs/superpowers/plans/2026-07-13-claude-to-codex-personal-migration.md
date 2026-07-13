# Claude-to-Codex Personal Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Install relevant self-contained Claude workflows into personal Codex, install the locally runnable FiftyOne plugin, and remove only verified stale or redundant Claude plugin state.

**Architecture:** Build and validate everything under `/tmp/codex-claude-migration` before touching live personal configuration. Install standalone skills under `$HOME/.agents/skills`, install FiftyOne through the default personal marketplace, then perform registry-aware Claude cache cleanup from a generated allowlist. Remove migration staging and the process-only repository documents after successful verification.

**Tech Stack:** Agent Skills (`SKILL.md`), Codex plugin manifests, JSON/TOML configuration, Bash, `jq`, Codex CLI 0.144.1, Codex skill/plugin validators.

---

### Task 1: Revalidate Sources And Create A Staging Snapshot

**Files:**
- Read: `/nfsshare/users/P126156127/.claude/plugins/installed_plugins.json`
- Read: `/nfsshare/users/P126156127/.claude/settings.json`
- Create: `/tmp/codex-claude-migration/inventory/claude-installed-plugins.json`
- Create: `/tmp/codex-claude-migration/inventory/claude-settings.json`
- Create: `/tmp/codex-claude-migration/inventory/preflight.txt`

- [ ] **Step 1: Confirm live roots and tool versions**

```bash
test -d "$HOME/.claude/skills"
test -d "$HOME/.agents/skills"
test -d "$HOME/.codex/skills"
codex --version
jq --version
```

Expected: all directory checks return zero; Codex reports `0.144.1` or a compatible newer version; `jq` is available.

- [ ] **Step 2: Create a private staging root and snapshot non-secret registries**

```bash
install -d -m 700 /tmp/codex-claude-migration/inventory
cp "$HOME/.claude/plugins/installed_plugins.json" /tmp/codex-claude-migration/inventory/claude-installed-plugins.json
cp "$HOME/.claude/settings.json" /tmp/codex-claude-migration/inventory/claude-settings.json
chmod 600 /tmp/codex-claude-migration/inventory/*.json
```

Expected: two mode-0600 snapshots exist. Do not copy `.credentials.json`, `.claude.json`, Codex auth, or MCP authorization headers.

- [ ] **Step 3: Record the current effective state**

```bash
codex plugin list > /tmp/codex-claude-migration/inventory/codex-plugins-before.txt
codex plugin marketplace list > /tmp/codex-claude-migration/inventory/codex-marketplaces-before.txt
codex mcp list > /tmp/codex-claude-migration/inventory/codex-mcp-before.txt
find "$HOME/.agents/skills" "$HOME/.codex/skills" -name SKILL.md -type f -print | sort > /tmp/codex-claude-migration/inventory/codex-skills-before.txt
```

Expected: Ponytail is the only installed Codex plugin, Plan is the only active global MCP server, and existing skill paths are captured for collision checks.

### Task 2: Stage And Validate Directly Portable Skills

**Files:**
- Read: the 35 named directories under `/nfsshare/users/P126156127/.claude/skills/` listed below
- Create: the same 35 named directories under `/tmp/codex-claude-migration/skills/`

Direct skill names:

```text
agent-harness-construction agent-introspection-debugging agentic-engineering
ai-regression-testing architecture-decision-records article-writing
benchmark-methodology benchmark-optimization-loop blueprint code-tour
coding-standards competitive-report-structure content-hash-cache-pattern
cost-aware-llm-pipeline data-scraper-agent data-throughput-accelerator
dmux-workflows error-handling growth-log intent-driven-development
iterative-retrieval manim-video market-research mle-workflow
parallel-execution-optimizer python-patterns python-testing pytorch-patterns
regex-vs-llm-structured-text scientific-thinking-literature-review
scientific-thinking-scholar-evaluation sketchnote team-agent-orchestration
terminal-ops token-budget-advisor
```

- [ ] **Step 1: Fail on any destination name collision**

```bash
DIRECT_SKILLS='agent-harness-construction agent-introspection-debugging agentic-engineering ai-regression-testing architecture-decision-records article-writing benchmark-methodology benchmark-optimization-loop blueprint code-tour coding-standards competitive-report-structure content-hash-cache-pattern cost-aware-llm-pipeline data-scraper-agent data-throughput-accelerator dmux-workflows error-handling growth-log intent-driven-development iterative-retrieval manim-video market-research mle-workflow parallel-execution-optimizer python-patterns python-testing pytorch-patterns regex-vs-llm-structured-text scientific-thinking-literature-review scientific-thinking-scholar-evaluation sketchnote team-agent-orchestration terminal-ops token-budget-advisor'
for name in $DIRECT_SKILLS; do
  test ! -e "$HOME/.agents/skills/$name"
  test ! -e "$HOME/.codex/skills/$name"
done
```

Expected: every check returns zero. If a path now exists, compare it and remove that name from the install set rather than overwriting it.

- [ ] **Step 2: Copy complete skill directories into staging**

```bash
install -d -m 700 /tmp/codex-claude-migration/skills
DIRECT_SKILLS='agent-harness-construction agent-introspection-debugging agentic-engineering ai-regression-testing architecture-decision-records article-writing benchmark-methodology benchmark-optimization-loop blueprint code-tour coding-standards competitive-report-structure content-hash-cache-pattern cost-aware-llm-pipeline data-scraper-agent data-throughput-accelerator dmux-workflows error-handling growth-log intent-driven-development iterative-retrieval manim-video market-research mle-workflow parallel-execution-optimizer python-patterns python-testing pytorch-patterns regex-vs-llm-structured-text scientific-thinking-literature-review scientific-thinking-scholar-evaluation sketchnote team-agent-orchestration terminal-ops token-budget-advisor'
for name in $DIRECT_SKILLS; do
  cp -a "$HOME/.claude/skills/$name" "/tmp/codex-claude-migration/skills/$name"
done
```

Expected: each staged directory contains `SKILL.md` and any referenced `scripts/`, `references/`, or `assets/` from its source.

- [ ] **Step 3: Validate every direct skill**

```bash
for skill in /tmp/codex-claude-migration/skills/*; do
  python3 "$HOME/.codex/skills/.system/skill-creator/scripts/quick_validate.py" "$skill"
done
```

Expected: every skill reports valid. Remove a failing skill from staging only after recording the validator error; do not install invalid content.

### Task 3: Create Codex-Native Adaptations

**Files:**
- Read: `/nfsshare/users/P126156127/.claude/skills/{codebase-onboarding,context-budget,council,eval-harness,goal,opensource-pipeline,prompt-optimizer,repo-scan,safety-guard,santa-method}/`
- Read: `/nfsshare/users/P126156127/.claude/plugins/cache/claude-plugins-official/claude-md-management/1.0.0/skills/claude-md-improver/`
- Create: `/tmp/codex-claude-migration/skills/{codebase-onboarding,context-budget,council,eval-harness,goal,opensource-pipeline,prompt-optimizer,repo-scan,safety-guard,santa-method,agents-md-improver}/`

- [ ] **Step 1: Copy the eleven source workflows into staging**

```bash
for name in codebase-onboarding context-budget council eval-harness goal opensource-pipeline prompt-optimizer repo-scan safety-guard santa-method; do
  cp -a "$HOME/.claude/skills/$name" "/tmp/codex-claude-migration/skills/$name"
done
cp -a "$HOME/.claude/plugins/cache/claude-plugins-official/claude-md-management/1.0.0/skills/claude-md-improver" /tmp/codex-claude-migration/skills/agents-md-improver
```

Expected: eleven isolated working copies exist.

- [ ] **Step 2: Apply the Codex adaptation contract**

Edit only the staged copies using these exact semantic mappings:

```text
CLAUDE.md                    -> AGENTS.md
$HOME/.claude/skills         -> $HOME/.agents/skills
Claude Code session/context  -> Codex thread/context
Task subagents               -> Codex collaboration subagents
TodoWrite                    -> update_plan
AskUserQuestion              -> concise user question/request_user_input when available
Read/Glob/Grep/Bash/Edit      -> Codex file search, exec_command, and apply_patch tools
Skill tool invocation        -> load/follow the named installed skill
```

For `agents-md-improver`, set frontmatter `name: agents-md-improver` and make all trigger text and output paths refer to `AGENTS.md`. Preserve safeguards and workflow intent. Remove Claude-only allowed-tools frontmatter where it cannot map to Codex. Do not invent MCP servers or add credential requirements.

- [ ] **Step 3: Prove operative Claude-only assumptions are gone**

```bash
rg -n 'CLAUDE\.md|\.claude/skills|TodoWrite|AskUserQuestion|EnterPlanMode|ExitPlanMode|allowed-tools:.*(Task|Read|Write|Edit|Bash|Glob|Grep)' /tmp/codex-claude-migration/skills/{codebase-onboarding,context-budget,council,eval-harness,goal,opensource-pipeline,prompt-optimizer,repo-scan,safety-guard,santa-method,agents-md-improver}
```

Expected: no operative matches. Historical comparison text is permitted only when explicitly labeled non-operative.

- [ ] **Step 4: Validate every adapted skill**

```bash
for name in codebase-onboarding context-budget council eval-harness goal opensource-pipeline prompt-optimizer repo-scan safety-guard santa-method agents-md-improver; do
  python3 "$HOME/.codex/skills/.system/skill-creator/scripts/quick_validate.py" "/tmp/codex-claude-migration/skills/$name"
done
```

Expected: all eleven report valid.

### Task 4: Stage The Native FiftyOne Plugin And Personal Marketplace

**Files:**
- Read: `/nfsshare/users/P126156127/.claude/plugins/cache/claude-plugins-official/fiftyone/1.0.3/`
- Create: `/tmp/codex-claude-migration/plugins/fiftyone-skills/`
- Create: `/tmp/codex-claude-migration/.agents/plugins/marketplace.json`

- [ ] **Step 1: Verify the local runtime and native manifest**

```bash
command -v fiftyone
command -v fiftyone-mcp
jq -e '.name == "fiftyone-skills" and .skills and .mcpServers' "$HOME/.claude/plugins/cache/claude-plugins-official/fiftyone/1.0.3/.codex-plugin/plugin.json"
```

Expected: both commands resolve locally and the manifest assertion returns true.

- [ ] **Step 2: Stage only runtime plugin components**

```bash
install -d -m 700 /tmp/codex-claude-migration/plugins/fiftyone-skills
cp -a "$HOME/.claude/plugins/cache/claude-plugins-official/fiftyone/1.0.3/.codex-plugin" /tmp/codex-claude-migration/plugins/fiftyone-skills/
cp -a "$HOME/.claude/plugins/cache/claude-plugins-official/fiftyone/1.0.3/.mcp.json" /tmp/codex-claude-migration/plugins/fiftyone-skills/
cp -a "$HOME/.claude/plugins/cache/claude-plugins-official/fiftyone/1.0.3/skills" /tmp/codex-claude-migration/plugins/fiftyone-skills/
cp -a "$HOME/.claude/plugins/cache/claude-plugins-official/fiftyone/1.0.3/scripts" /tmp/codex-claude-migration/plugins/fiftyone-skills/
```

Expected: no `.git`, `.claude-plugin`, cache marker, or marketplace checkout metadata is copied.

- [ ] **Step 3: Create the staged personal marketplace entry**

```json
{
  "name": "personal",
  "interface": { "displayName": "Personal" },
  "plugins": [
    {
      "name": "fiftyone-skills",
      "source": { "source": "local", "path": "./plugins/fiftyone-skills" },
      "policy": { "installation": "AVAILABLE", "authentication": "ON_INSTALL" },
      "category": "Developer Tools"
    }
  ]
}
```

Save this exact structured JSON at `/tmp/codex-claude-migration/.agents/plugins/marketplace.json` and validate it with `jq -e .`.

- [ ] **Step 4: Validate the plugin and its bundled skills**

```bash
python3 "$HOME/.codex/skills/.system/plugin-creator/scripts/validate_plugin.py" /tmp/codex-claude-migration/plugins/fiftyone-skills
for skill in /tmp/codex-claude-migration/plugins/fiftyone-skills/skills/*; do
  python3 "$HOME/.codex/skills/.system/skill-creator/scripts/quick_validate.py" "$skill"
done
```

Expected: the plugin and all 16 skills report valid.

### Task 5: Install The Validated Personal Skills And Plugin

**Files:**
- Create: the 46 validated staged skill directories under `/nfsshare/users/P126156127/.agents/skills/`
- Create: `/nfsshare/users/P126156127/plugins/fiftyone-skills/`
- Create: `/nfsshare/users/P126156127/.agents/plugins/marketplace.json`
- Modify: `/nfsshare/users/P126156127/.codex/config.toml` through `codex plugin add`

- [ ] **Step 1: Re-run destination collision checks**

```bash
for skill in /tmp/codex-claude-migration/skills/*; do
  name="${skill##*/}"
  test ! -e "$HOME/.agents/skills/$name"
  test ! -e "$HOME/.codex/skills/$name"
done
test ! -e "$HOME/plugins/fiftyone-skills"
test ! -e "$HOME/.agents/plugins/marketplace.json"
```

Expected: all checks return zero. Stop rather than overwrite if live state changed.

- [ ] **Step 2: Install standalone skills atomically**

```bash
for skill in /tmp/codex-claude-migration/skills/*; do
  cp -a "$skill" "$HOME/.agents/skills/${skill##*/}"
done
```

Expected: 46 new global skill directories appear with preserved resource files.

- [ ] **Step 3: Install the plugin source and marketplace**

```bash
install -d "$HOME/plugins" "$HOME/.agents/plugins"
cp -a /tmp/codex-claude-migration/plugins/fiftyone-skills "$HOME/plugins/fiftyone-skills"
cp /tmp/codex-claude-migration/.agents/plugins/marketplace.json "$HOME/.agents/plugins/marketplace.json"
codex plugin add fiftyone-skills@personal --json
```

Expected: Codex reports a successful local install and enables `fiftyone-skills@personal`.

### Task 6: Remove Verified Duplicate Claude Plugin State

**Files:**
- Modify: `/nfsshare/users/P126156127/.claude/plugins/installed_plugins.json`
- Preserve backup until final verification: `/tmp/codex-claude-migration/backup/installed_plugins.json`
- Remove: verified stale version/cache, orphan, temporary checkout, and staging paths under `/nfsshare/users/P126156127/.claude/plugins/`

- [ ] **Step 1: Generate the cleanup allowlist from live references**

Create `/tmp/codex-claude-migration/inventory/cleanup-candidates.txt` with these
fixed stale paths plus the currently present top-level `temp_git_*` directories:

```bash
printf '%s\n' \
  "$HOME/.claude/plugins/cache/claude-plugins-official/superpowers/5.1.0" \
  "$HOME/.claude/plugins/cache/claude-plugins-official/superpowers/6.0.3" \
  "$HOME/.claude/plugins/cache/claude-plugins-official/superpowers/6.1.0" \
  "$HOME/.claude/plugins/cache/claude-plugins-official/remember/0.7.3" \
  "$HOME/.claude/plugins/cache/claude-plugins-official/remember/0.8.2" \
  "$HOME/.claude/plugins/cache/claude-plugins-official/context7/33d632e4e734" \
  "$HOME/.claude/plugins/cache/claude-plugins-official/skill-creator/33d632e4e734" \
  "$HOME/.claude/plugins/cache/claude-plugins-official/code-review/33d632e4e734" \
  "$HOME/.claude/plugins/cache/claude-plugins-official/github/33d632e4e734" \
  "$HOME/.claude/plugins/cache/claude-plugins-official/playwright/33d632e4e734" \
  "$HOME/.claude/plugins/cache/claude-plugins-official/mcp-server-dev/33d632e4e734" \
  "$HOME/.claude/plugins/cache/ecc/ecc/2.0.0" \
  "$HOME/.claude/plugins/data/ecc-ecc" \
  "$HOME/.claude/plugins/marketplaces/claude-plugins-official.staging" \
  > /tmp/codex-claude-migration/inventory/cleanup-candidates.txt
find "$HOME/.claude/plugins/cache" -mindepth 1 -maxdepth 1 -type d -name 'temp_git_*' -print \
  | sort >> /tmp/codex-claude-migration/inventory/cleanup-candidates.txt
```

For every path, assert it is absent from all current `installPath` values:

```bash
jq -r '.plugins | to_entries[].value[].installPath' "$HOME/.claude/plugins/installed_plugins.json" > /tmp/codex-claude-migration/inventory/current-install-paths.txt
while IFS= read -r candidate; do
  test -e "$candidate"
  ! rg -Fx -- "$candidate" /tmp/codex-claude-migration/inventory/current-install-paths.txt
done < /tmp/codex-claude-migration/inventory/cleanup-candidates.txt
```

Expected: every candidate exists and none is a current registry install path. Do not use a broad wildcard removal command against versioned plugin directories.

- [ ] **Step 2: Back up and normalize redundant install records**

```bash
install -d -m 700 /tmp/codex-claude-migration/backup
cp "$HOME/.claude/plugins/installed_plugins.json" /tmp/codex-claude-migration/backup/installed_plugins.json
chmod 600 /tmp/codex-claude-migration/backup/installed_plugins.json
```

Use this structured transform to remove only the two redundant project records:

```bash
jq '
  .plugins["code-simplifier@claude-plugins-official"] |= map(select(.scope != "project")) |
  .plugins["github@claude-plugins-official"] |= map(select(.scope != "project"))
' "$HOME/.claude/plugins/installed_plugins.json" > "$HOME/.claude/plugins/installed_plugins.json.new"
jq -e '
  (.plugins["code-simplifier@claude-plugins-official"] | length) == 1 and
  (.plugins["code-simplifier@claude-plugins-official"][0].scope) == "user" and
  (.plugins["github@claude-plugins-official"] | length) == 1 and
  (.plugins["github@claude-plugins-official"][0].scope) == "user"
' "$HOME/.claude/plugins/installed_plugins.json.new"
mv "$HOME/.claude/plugins/installed_plugins.json.new" "$HOME/.claude/plugins/installed_plugins.json"
```

Expected: each affected plugin retains exactly its user-scope record and the same current install path.

- [ ] **Step 3: Delete only allowlisted candidates**

```bash
while IFS= read -r candidate; do
  test -n "$candidate"
  case "$candidate" in "$HOME/.claude/plugins/"*) ;; *) exit 1 ;; esac
  test ! -L "$candidate"
  rm -rf -- "$candidate"
done < /tmp/codex-claude-migration/inventory/cleanup-candidates.txt
```

Expected: every listed stale path is removed; current enabled cache paths remain.

- [ ] **Step 4: Verify all enabled Claude plugins still resolve**

```bash
jq -r '.plugins | to_entries[].value[].installPath' "$HOME/.claude/plugins/installed_plugins.json" | while IFS= read -r path; do test -d "$path"; done
jq -e . "$HOME/.claude/plugins/installed_plugins.json"
```

Expected: all current registry paths exist and the registry remains valid JSON. Restore the backup immediately if either check fails.

### Task 7: End-To-End Verification And Artifact Cleanup

**Files:**
- Verify: personal skill, plugin, marketplace, MCP, and Claude registry state
- Remove: `/tmp/codex-claude-migration/`
- Remove: `docs/superpowers/specs/2026-07-13-claude-to-codex-personal-migration-design.md`
- Remove: `docs/superpowers/plans/2026-07-13-claude-to-codex-personal-migration.md`

- [ ] **Step 1: Validate installed skills and unique names**

```bash
for skill in /tmp/codex-claude-migration/skills/*; do
  python3 "$HOME/.codex/skills/.system/skill-creator/scripts/quick_validate.py" "$HOME/.agents/skills/${skill##*/}"
done
```

Extract frontmatter names from active system, shared, personal, and enabled plugin skill roots; expected result is no newly introduced duplicate name.

- [ ] **Step 2: Validate installed plugin state**

```bash
python3 "$HOME/.codex/skills/.system/plugin-creator/scripts/validate_plugin.py" "$HOME/plugins/fiftyone-skills"
codex plugin marketplace list
codex plugin list
codex mcp list
```

Expected: personal marketplace is discovered, `fiftyone-skills@personal` is enabled, Ponytail remains enabled, and no unrelated global MCP server was added.

- [ ] **Step 3: Run local runtime smoke checks**

```bash
fiftyone --help
fiftyone-mcp --help
```

Expected: both return help successfully without authentication prompts or downloads.

- [ ] **Step 4: Remove temporary and process-only artifacts**

After all verification passes, remove `/tmp/codex-claude-migration` and the two process documents created in this repository. Remove now-empty `docs/superpowers/plans`, `docs/superpowers/specs`, and `docs/superpowers` directories only if they contain no pre-existing files. Commit only the documentation removal; do not stage unrelated workspace files.

- [ ] **Step 5: Confirm final workspace and personal state**

```bash
git status --short
codex plugin list
codex mcp list
jq -e . "$HOME/.claude/plugins/installed_plugins.json"
```

Expected: only the user's pre-existing untracked workspace files remain; the intended personal Codex additions exist; Claude's registry is valid; migration staging and process-only documents are gone.

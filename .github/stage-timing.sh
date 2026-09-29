# Wall time per pipeline stage, recorded rather than read off the Actions UI.
#
# `src/vdjdb/timing.py` measures the stages *inside* one command. The stages a build is actually
# budgeted in are the commands themselves - proofreading, assembly, motifs, dashboard, the release
# comparison - and those were visible only in the web UI, so every figure in ROADMAP_local came from
# someone reading it by hand. `CLAUDE.md` asks for the opposite: a green build that takes 25 minutes
# is a measurement nobody has written down.
#
# Usage, from a workflow step:
#     . .github/stage-timing.sh
#     stage assemble uv run vdjdb build --out out/
#
# The exit status is the command's, so a failing stage still fails the step. A stage that fails is
# recorded too: the time it took to fail is the only clue about *where* it failed.
STAGE_TIMINGS="${STAGE_TIMINGS:-out/reports/ci-stage-timings.tsv}"

stage() {
  local name="$1"; shift
  local started ended status
  mkdir -p "$(dirname "$STAGE_TIMINGS")"
  [ -s "$STAGE_TIMINGS" ] || printf 'stage\tseconds\tstatus\n' > "$STAGE_TIMINGS"
  started=$(date +%s)
  "$@" && status=ok || status=failed
  ended=$(date +%s)
  printf '%s\t%s\t%s\n' "$name" "$((ended - started))" "$status" >> "$STAGE_TIMINGS"
  [ "$status" = ok ]
}

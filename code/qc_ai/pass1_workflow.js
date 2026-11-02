export const meta = {
  name: 'leadchi-qc-pass1',
  description: 'LEAD-chi visual QC pass 1: 180 coded sessions graded by independent AI rater instances (12 batches x 15), amendment 1c',
  phases: [{ title: 'Pass1', detail: '12 independent rater instances, 15 sessions each, sheet order' }],
}
const M = '/path/to/scratch/qc_montages'
const BATCHES = args.batches
const SCHEMA = {
  type: 'object',
  properties: {
    rater_notes: { type: 'string' },
    grades: { type: 'array', items: { type: 'object', properties: {
      qc_id: { type: 'string' },
      grade: { type: 'integer', enum: [0, 1, 2] },
      failure_codes: { type: 'string', description: 'semicolon-separated codes from {M,W,S,B,K,F,A,C,D}; empty if none' },
      placement: { type: 'string', enum: ['pass', 'fail'] },
      notes: { type: 'string' },
    }, required: ['qc_id', 'grade', 'failure_codes', 'placement', 'notes'] } },
  },
  required: ['grades', 'rater_notes'],
}
const PROMPT = (k, items) => `
You are visual-QC rater instance ${k} (pass 1) for a blinded QSM (quantitative susceptibility mapping) study. You grade coded
sessions STRICTLY by the written rubric. You know nothing about any participant (no group, age, sex, diagnosis) and you must
NOT try to find out: read ONLY the rubric file and the montage PNGs listed below. Do not open any other file or directory,
do not run shell commands other than viewing these files, do not search the web.
1) Read the rubric verbatim: ${M}/RUBRIC_sections1-3.md
2) For EACH session below, view BOTH montages with the Read tool (they are images):
   gross:     ${M}/montages/gross/<qc_id>.png      (rows: dentate / midbrain / basal-ganglia levels; columns: mag e1, mag e4,
              total field, local field, chi with the Mask-3 outline in yellow)
   placement: ${M}/montages/placement/<qc_id>.png  (top row T1w, bottom row echo-1 magnitude; red outlines = SN-VTA atlas
              region in axial/coronal/sagittal and the dentate outlines in the right panel)
   forced_review = 1 means the automated checks flagged the session: inspect every panel at full size (use the Read tool's
   image view carefully) before grading; the flag itself does not change the grade. Inspect carefully in all cases.
   Sessions (qc_id : forced_review): ${items.map(i => `${i.qc_id}:${i.forced_review}`).join(', ')}
3) Apply rubric section 2 for the gross grade (0 pass / 1 minor / 2 fail) and the failure codes: grade 2 ONLY when a failure
   code is present at the midbrain or dentate level or across the whole slice; otherwise 1 if an artefact is visible but spares
   those levels or is faint; 0 if none. Record code D whenever dental/metal artefact is visible, whatever the grade.
   Apply rubric section 3 for placement (pass/fail) using ONLY the T1w and echo-1 magnitude panels, never chi.
4) Use the same thresholds for every session. One short factual note per session (what you saw, which panel/level).
Return the structured output for ALL ${items.length} sessions, in the order given. In rater_notes, mention any rubric
ambiguity you met (do not invent new criteria).`
phase('Pass1')
const out = await parallel(BATCHES.map((b, i) => () => agent(PROMPT(i + 1, b), { label: `pass1_batch_${i + 1}`, phase: 'Pass1', schema: SCHEMA })
  .then(r => ({ batch: i + 1, result: r }))))
return out

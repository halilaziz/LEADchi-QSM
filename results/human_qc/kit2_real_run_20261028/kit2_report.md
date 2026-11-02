# LEAD-χ human visual QC, kit 2: complete-human-QC result (POST HOC)

Plan: `PLAN_kit2.md` (SHA-256 `ec4f2f3a5291e97235fb9aa1e9021d23a54bc8d7b7ce75954403593c809ee099`). Rater: first author (not independent).

## Coverage

- AI_graded_sessions: 180
- sessions_with_human_grade: 180
- sessions_with_human_placement: 180
- kit1_sessions: 101
- kit2_new_sessions: 79
- kit2_new_without_grade: 0
- kit2_new_without_placement: 0
- kit1_sessions_without_grade: 0
- complete_human_qc: True

## Checks

- frozen_reconstruction: {'usable_equal': True, 'first_usable_equal': True, 'fail_reason_equal': True}
- AI_frozen_reproduces_confirmatory: True
- kit1_ratings_are_real_returned_csv: True
- K1_H_replace_reproduces_n89_run: True

## M1 (frozen helpers; bound test at the frozen bound −6.1318 ppb)

| Scenario | n | β_PS (ppb) | 95 % CI | p (two-sided) | p_B | label |
|---|---|---|---|---|---|---|
| AI_frozen | 80 | 4.963 | [-1.273, 11.199] | 0.1169 | 3.52e-04 | B |
| K1_H_replace | 89 | 6.290 | [0.176, 12.405] | 0.0439 | 6.16e-05 | B |
| H_all180 | 85 | 6.281 | [-0.226, 12.787] | 0.0583 | 1.48e-04 | B |

`AI_frozen` = registered primary (n = 80); `K1_H_replace` = kit-1 grades only (the n = 89 run when the real kit-1 CSV is used); **`H_all180` = this plan's primary**. Reported next to the registered result; the registered outcome stays B.

**Interpretation (plan section 5):** the result is robust to replacing the AI visual-QC grade with a first-author human grade for all 180 sessions.

### Sample changes of H_all180

- vs_AI_frozen: added 9, removed 4, changed first session 2
- vs_K1_H_replace: added 0, removed 4, changed first session 0
- SN–VTA re-check excluded: 0

## Test-retest (21 hidden repeats vs kit-1 grade)

- pairs with both grades: 21 / 21
- exact agreement 0.857 (Wilson 95 % CI 0.654–0.950); linear-weighted κ 0.714
- table (rows kit 2, columns kit 1, grades 0/1/2): [[14, 1, 0], [1, 1, 0], [1, 0, 3]]
- placement identical: 20 / 21
- by item type: {'repeat_of_kit1_AI_fail': {'n': 11, 'identical': 8}, 'repeat_of_kit1_AI_pass': {'n': 10, 'identical': 10}}
- Descriptive only; recall is possible, so agreement may be optimistic.

## Secondary (descriptive): census agreement human vs AI

- without the example-montage session (primary): n 179, linear κ 0.068, exact 0.134; table (rows human, cols AI) [[7, 104, 44], [0, 7, 6], [0, 1, 10]]
- with it: n 180, linear κ 0.067, exact 0.133
- human grade 2 among AI-pass sessions: 1 / 119 = 0.008 (Wilson 95 % CI 0.001–0.046)


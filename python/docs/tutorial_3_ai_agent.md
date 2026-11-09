# Tutorial 3: using RAST with an AI coding agent

An AI coding agent (Claude Code, Codex, Cursor, GitHub Copilot in agent mode ...) can install RAST, run it, organise
parameter studies and summarise the results tables for you. This tutorial explains how to set one up so that it
knows how RAST works, and what to keep for yourself.

## What the agent can and can't do

| The agent can | You keep |
|---|---|
| install Python packages and RAST, run `rast demo` | **drawing the masks**: what is an indent, a scratch or a slip trace in your pictures |
| run `rast run` with the parameters you ask for, on one or many datasets | **the Euler angles** and their convention: they come from your EBSD data |
| build a dataset from files you give it (pictures, mask PNGs, an Euler CSV) | deciding the score threshold and parameters you report |
| read the tables, sum counts per grain, make plots, compare runs | the scientific interpretation of the match fractions |
| point out warning signs (many "no match" lines, results close to the test lines) | |

The repository tells agents to follow these rules (in `AGENTS.md`).

## Setting up

1. Get the repository: `git clone https://github.com/vgagneur/RAST`, or download it as a zip from GitHub.
2. Open the `RAST` folder in your agent (e.g. run `claude` in that folder, or open the folder in your editor's
   agent).
3. That's all. The agent reads **`AGENTS.md`** at the top of the repository (Claude Code reads it through
   `CLAUDE.md`). It holds the notes the agent needs: where things are, the input format, the conventions (Bunge Euler
   angles, crystal modes, the ST rule), the rules above, and links to more detailed notes in
   `python/docs/agent-notes/`. These notes travel with the repository, so an agent on any computer has them without
   any setup or earlier conversation.

Keep your data outside the repository folder (e.g. `D:\RAST_data\`) and give the agent its path.

## Example requests

- "Install RAST and run the demo. Tell me if pm110 at ST 10 is close to 100 %."
- "Run RAST on `D:\RAST_data\TiFe\TiFe_ref.mat` with degrot 157, figures at ST 10 only, 4 workers."
- "Build a dataset from the pictures in `D:\RAST_data\NiAl\pictures`, the masks in `...\masks` and the angles in
  `angles.csv`, then run it with the defaults."
- "Run `my_sample.mat` with N = 80, 100, 120 and 140, no figures, and make a plot of the {110} fraction per grain
  against N at ST = 10."
- "For each grain, give me the fraction of lines matching each family at ST 10, next to the test-line fraction."
- "Open the mask window for pictures 3 and 5 of `my_sample.mat` so I can extend the masks; save as
  `my_sample_v2.mat`." (The windows open on your screen; you draw.)

## What to check in the agent's answers

- **Parameters**: the answer should say which N, IR, ST, crystal, nsmooth and degrot were used (they are also in
  `params.json` in each output folder).
- **Sums, not averages**: per-grain fractions should be computed from summed counts (sum of `m110` / sum of `nlines`
  over the grain's pictures), not from averaged percentages.
- **Warnings**: the agent should mention many "no match" lines or real-line results close to the test lines. If it
  doesn't, ask.
- **Nothing overwritten**: each run should have its own output folder.
- **Masks and angles unchanged** unless you asked for a specific change.

## When the agent works on a remote computer

The mask windows need a desktop session. On a remote machine or in the cloud, do the masking step yourself
(`rast prepare` on your computer, or masks drawn in ImageJ/GIMP saved as `<picture>_mask.png`), then let the agent
run the rest in file mode ([Tutorial 1, "Masks drawn in another program"](tutorial_1_new_dataset.md#masks-drawn-in-another-program)).

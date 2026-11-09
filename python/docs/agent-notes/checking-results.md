# Checking and summarising results (agent note)

For each sample and parameter set, report to the user (state N, IR, ST, crystal, nsmooth, degrot):
- the number of grains and pictures, and `nlines` per grain at the chosen ST;
- for each grain, the fraction of lines matching each family: the sum of `mXXX` divided by the sum of `nlines` over
  the grain's pictures (sum the counts; don't average percentages over pictures with different numbers of lines),
  next to the test-line fraction `ptmXXX` (the chance level for that orientation).

Warning signs to tell the user about. Don't hide them, and don't "fix" them yourself:
1. `pnomatch` clearly above 0 at the working ST: the Euler angles, their convention or frame, or `degrot` are
   probably wrong, or IR is too small.
2. Real-line fractions close to the test-line fractions for every family: no detectable preference. Check the edges
   in `<i>_bw.png` and the `picture_*_linesselected.png` figures, then the angles.
3. Lines in the figures running along scratches, grain boundaries or the indent edge: masking is needed, which is
   the user's job.
4. One picture very different from the others of the same grain: possibly a wrong Euler entry or bad detection on
   that picture. Show the user which one.
5. A known sample that doesn't give its known answer (Ti2AlMo ~100 % {110}; pure fcc Ni ~100 % {111}).

RAST outputs match counts per family. Deciding from them which slip systems are active is the user's scientific
judgement: present numbers and warnings, not conclusions, unless asked.

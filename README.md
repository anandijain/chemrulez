# chemrulez

First iteration of an organic chemistry synthesis sandbox.

## Run

Serve the directory with any static file server:

```sh
python3 -m http.server 5173
```

Then open:

```txt
http://localhost:5173
```

## Test

```sh
node test/rule-engine.test.js
```

Debug one pathway from the terminal:

```sh
node scripts/pathway.js "3,3-dimethyl-1-butene" "HBr"
```

## What works now

- Import molecules by common name, SMILES, PubChem CID, or PubChem compound URL.
- Render molecule structures using PubChem images.
- Fuzzy-resolve `NaNH2`, `sodium amide`, and `sodamide`.
- Resolve structural reagents through PubChem, including arbitrary named alkyl halides.
- Apply the first rule: terminal alkyne + sodium amide gives an acetylide candidate.
- Apply acetylide SN2 alkylation with methyl, primary, allylic, and benzylic halides detected from structure.
- Enter one-pot-ish sequences like `NaNH2 then CH3Br` from a terminal alkyne.
- Warn on secondary alkyl halides and block tertiary alkyl halides for acetylide alkylation.
- Apply alkyne reductions: `H2, Lindlar` to cis alkene, `Na, NH3` to trans alkene, and `H2, Pd/C` to alkane.
- Apply alkyne hydration: `HgSO4, H2SO4, H2O` to ketone and bulky hydroboration/oxidation to aldehyde for terminal alkynes.
- Apply common alkene reactions: hydrogenation, HBr, HBr/ROOR, acid hydration, oxymercuration, hydroboration, Br2, mCPBA, and OsO4.
- Apply reductive ozonolysis with `O3 then DMS` or similar entries.
- Apply first-pass Grignard additions to aldehydes/ketones and CO2 carboxylation.
- Show first-pass carbocation rearrangement fan-out for acid hydration and HBr additions.
- Click a product candidate to continue the synthesis path.

## Next useful iteration

- Add a reagent/rule JSON file instead of hard-coding rules in `src/app.js`.
- Replace string-based product generation with graph-based transforms.
- Move chemistry validation to RDKit once the interaction model feels right.

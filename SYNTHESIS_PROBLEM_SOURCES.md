# Canonical synthesis problem sources

This file tracks source material for a synthesis-puzzle test set. Prefer open,
stable sources that explain undergraduate consensus rules over unsourced problem
lists.

## Source spine

1. OpenStax Organic Chemistry
   - Use for canonical reaction families and end-of-chapter problem themes.
   - Good immediate chapters: alkenes, alkynes, alcohols, ethers/epoxides,
     aldehydes/ketones, carboxylic acids, aromatic substitution.
   - Epoxide formation and opening are especially useful because they force the
     graph engine to handle rings, regioselectivity, and follow-on products.
  - URLs:
    - https://openstax.org/books/organic-chemistry/pages/9-9-an-introduction-to-organic-synthesis
     - https://openstax.org/books/organic-chemistry/pages/18-4-cyclic-ethers-epoxides
     - https://openstax.org/books/organic-chemistry/pages/18-5-reactions-of-epoxides-ring-opening
     - https://openstax.org/books/organic-chemistry/pages/8-additional-problems

2. Organic Chemistry Data / Reusch Virtual Textbook practice problems
   - Use for multi-step synthesis buckets and route-planning problem categories.
   - The list includes introduction to multistep synthesis, predicting products
     from multistep syntheses, devising multistep syntheses, and named
     "A Multistep Synthesis" exercises.
   - URL:
     - https://organicchemistrydata.org/reusch/virtualtext/problems/practice-problems/

3. Organic Chemistry 1: An Open Textbook / Lumen
   - Use for concise reaction summaries and first-semester prioritization.
   - Its Part 1 reaction summary is useful for deciding which reactions belong in
     the first playable puzzle pack.
   - URL:
     - https://courses.lumenlearning.com/suny-potsdam-organicchemistry/back-matter/appendix-1/

4. Chemistry LibreTexts
   - Use as a cross-check when regioselectivity or mechanism wording is ambiguous.
   - Good immediate page:
     - https://chem.libretexts.org/Courses/University_of_Illinois_UrbanaChampaign/Chem_2363A_Fundamental_Organic_Chemistry_I_%28Chan%29/11%3A_The_Chemistry_of_Ethers_Epoxides_Glycols_and_Sulfides/11.04%3A_Opening_of_Epoxides

## Ranked puzzle coverage

### Tier 1: core first-semester synthesis loop

- OpenStax 9.9 alkyne synthesis examples:
  - Acetylene plus alkyl halide chain extension, partial/full reduction, and
    hydroboration-oxidation to alcohol targets.
  - 4-octyne-only-source problems: cis-4-octene, butanal, 4-bromooctane,
    4-octanol, 4,5-dichlorooctane, and butanoic acid.
  - Acetylene plus alkyl halide problems: decane, 2,2-dimethylhexane, hexanal,
    and 2-heptanone.
- Alkene -> epoxide -> diol or halohydrin.
- Alkene -> bromohydrin -> epoxide.
- Alkene regiochemistry contrast: HBr, HBr/ROOR, hydration, oxymercuration,
  hydroboration.
- Alkyne chain extension: terminal alkyne -> acetylide -> alkylation.
- Alkyne reductions and hydrations.
- Ozonolysis as target-disconnection practice.

### Tier 2: carbonyl and C-C bond formation

- Grignard addition to aldehydes and ketones.
- Grignard carboxylation.
- Oxidation/reduction ladder between alcohols, aldehydes, ketones, acids.

### Tier 3: substitution/elimination route planning

- Primary/secondary/tertiary substrate choice.
- SN2 vs E2.
- Alcohol -> leaving group -> substitution/elimination.

### Tier 4: aromatic synthesis

- Electrophilic aromatic substitution directing effects.
- Benzylic oxidation and reduction.

## Test-set policy

- Encode problems as start, target, allowed reagents, and expected intermediate
  products rather than as prose-only examples.
- Prefer graph-equivalent target checks over exact SMILES text.
- Keep source URLs with each puzzle group so chemistry disagreements can be
  resolved against a public reference.

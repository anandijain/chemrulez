# Chemrulez roadmap

## Product goals

- Ship as a static GitHub Pages app so it can be used without paid hosting.
- Make the mobile experience first-class: quick one-hand reagent entry, readable
  molecule/path cards, low-latency puzzle play, and no desktop-only assumptions.
- Make it feel like LeetCode for synthesis: start molecule, target molecule,
  allowed reagents, max steps, hints, path history, and solved-state feedback.
- Assume internet is available, but keep hosting static: GitHub Pages can call
  public chemistry APIs and load static WASM assets without a paid backend.

## Chemistry engine goals

- Offload chemical name resolution to PubChem or a real chemistry resolver.
  Local aliases should stay small and curated for textbook demos/tests, not grow
  into a homemade IUPAC/common-name parser.
- Treat SMILES/SMARTS as graph serializations or graph queries, not as the source
  of truth.
- Move reaction rules behind graph transforms and graph-pattern queries.
- Use reaction conditions as first-class inputs, not just reagent names.
  Temperature, solvent, nucleophile/base strength, concentration, and substrate
  class should affect SN1/SN2/E1/E2 ranking.
- Prefer candidate mixtures with explanations over pretending ambiguous
  conditions give one clean answer.
- Eventually support biochemistry/protein puzzles, but keep that as a separate
  engine layer: amino-acid/residue graphs, peptide cleavage/ligation/modification
  rules, protecting groups, and sequence-level views are different from small
  molecule first-year synthesis.

## Reaction-condition policy

- Do not encode "heat" as a cosmetic alias. Heat should bias elimination over
  substitution when the substrate/base/nucleophile combination supports it.
- Bulky strong bases should favor E2, especially Hofmann products.
- Strong unhindered nucleophiles with primary substrates should favor SN2.
- Weak nucleophiles in polar protic media with tertiary substrates should favor
  SN1/E1 mixtures, with heat increasing E1.
- Secondary substrates should usually produce ranked competing candidates rather
  than a single product unless the conditions are very diagnostic.

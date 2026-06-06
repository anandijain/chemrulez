# Reagent Set Model

Chemrulez should treat reagents as chemical participants, not as magic strings.

The useful mental model is an ordered list of stages. Each stage contains one or
more agents:

```js
[
  {
    order: 1,
    agents: [
      { role: "reagent", molecule: "BH3" },
      { role: "solvent", molecule: "THF" }
    ]
  },
  {
    order: 2,
    agents: [
      { role: "oxidant", molecule: "H2O2" },
      { role: "base", molecule: "NaOH" }
    ]
  }
]
```

Named reagent sets such as Lindlar, PCC, Jones oxidation, hydroboration-oxidation,
or Grignard addition are labels for common stage/agent bundles. Under the hood,
the agents should still be inspectable chemical species where possible, with
PubChem links, structure metadata, and roles.

The reaction rule is the edge in the synthesis graph. A rule should match:

- the current substrate graph,
- the ordered reagent stages,
- the agent roles and chemical identities,
- conditions such as heat, solvent, acid/base strength, and workup.

Reaction rules should transform molecular graph state, not rewrite SMILES strings.
SMILES and SMARTS are acceptable as serialization, import/export, fixture, and
eventual graph-query formats, but they should not be the product-generation
engine. Direct string replacement is brittle around branching, aromaticity,
stereochemistry, disconnected products, and equivalent SMILES spellings. If a
rule cannot be implemented as a graph transform yet, prefer returning an
explicit unsupported/no-product candidate over adding a case-specific string
patch.

This avoids making a hard distinction between "substrate molecules" and
"reagents". The UI can still show compact named sets, but the data model should
preserve the underlying chemicals so a user can inspect, copy, and debug a route.

Some transformations need a fixed condition plus a patterned structural partner.
Examples:

- Grignard addition: the current substrate or reagent may be the organomagnesium
  species, while the other participant may be a graph-classified aldehyde,
  ketone, CO2, ester, or incompatible proton source.
- Friedel-Crafts acylation: `AlCl3` is a fixed Lewis-acid condition, while the
  acid chloride is a graph-classified acyl donor.
- Future Wittig support: the ylide/phosphorane should be a structural reagent
  partner, not a hard-coded list of named ylides.

For these cases, parse the input into agents first, derive graph roles from each
agent, then let the rule engine match substrate role + reagent roles. Do not add
one-off aliases for every possible aldehyde, ketone, acid chloride, or ylide.

## Known Model Gaps

- Grignard reactions are only partially supported. The app can form simple
  Grignard reagents, apply some carbonyl additions, and flag carboxylic acid
  quench, but arbitrary organomagnesium reagents, ordered acid workups,
  incompatible functional groups, and ester/carbonyl partner coverage are not
  yet robust.
- Multi-component reagent sets need a more general data model. Friedel-Crafts
  acylation is the current example: `AlCl3` is a fixed condition/activator, while
  the acid chloride is a structural acyl donor. The current parser can recognize
  missing `AlCl3` for acid chloride-only input, but it is not yet a reliable
  abstraction for arbitrary custom reagent partners.
- Wittig reactions are not supported. A future implementation probably needs
  structural ylide/phosphorane reagent partners, carbonyl matching, alkene
  product generation, and explicit E/Z selectivity annotations.

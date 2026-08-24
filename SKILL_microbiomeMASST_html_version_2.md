---
name: skill-microbiome-version-html
description: Interpret microbiomeMASST network graph/tree images (the metadata graph output of the GNPS2 microbiomeMASST tool, showing a root compound branching into Host, Origin, and Microbial nodes with colored categories like Health phenotype, Organs and biofluids, Microbial information, In vitro screening, and Biological context/Interventions). Use this whenever a user uploads a microbiomeMASST results HTML page (the downloadable results file containing a `const root=` tree object), uploads or describes a microbiomeMASST output image, mentions a "MASST tree/graph", or asks what a metabolite's occurrence across hosts/organs/microbes means — even if they just say "what does this MASST result mean" without naming the tool explicitly. Provides a structured, appropriately hedged interpretation plus the caveats needed so file counts and branches aren't over-read as prevalence, causality, or confirmed structure.
---

# Interpreting microbiomeMASST graphs

microbiomeMASST graphs are a metadata *interpretation layer* over a spectral
library search, not direct proof of identity, prevalence, or causality. The
single most common failure mode when someone reads one of these graphs is
over-claiming: turning "a match existed in N deposited files" into "this
compound is common in X" or "microbe Y produces this in humans." Your job is
to describe what the graph actually shows, then attach the caveats needed to
keep the reader from over-claiming — without being so hedgy that the
interpretation becomes useless.

This file has two parts. **Part 1** is the condensed workflow to run on every
image. **Part 2** is the full source guide (11 sections, written by the guide's
author) — consult it when a case is unusual or the user wants the fuller
checklist (e.g. they want to design the next search, or they're pushing back on
a caveat you gave).

**Core interpretation rule (Part 2, §4).** Do not convert spectral similarity
into exact structural identity, metadata association into causality, file counts
into prevalence, or a non-match into biological absence. Every caveat below is a
special case of this line; if you check nothing else before sending, check this.

---

# Part 1 — Workflow for reading a tree image

## Inputs: the results HTML, or a screenshot, plus the molecule name

The user supplies **the molecule name** and **either** of:

1. **The microbiomeMASST results HTML page** (preferred). This is the file
   downloaded from the results view; it embeds the entire tree as a JavaScript
   object with per-node counts and per-file match records. Almost everything the
   screenshot workflow has to hedge about is *in this file*. Always parse it.
2. **A screenshot of the tree.** Interpret it under the constraints in
   "Reading a screenshot" below.

If both arrive, parse the HTML and use the image only for layout/labels.

- **Do not ask clarifying questions.** Not about tree mode, search thresholds, or
  which dataset. With HTML, the answers are in the file — extract them. With a
  screenshot, say in one line that the image can't settle it and move on. A
  question back is a failure of this skill, not a safeguard.
- **If the molecule name is missing**, interpret anyway using the root label, and
  note once that the root identity is unspecified.
- **The supplied name is an annotation, not a confirmed structure.** Use it as the
  label throughout, but never let it upgrade the evidence: a name plus a tree is
  still a spectral match. Hedge at the level the data supports.
- **Use the name to sharpen the read.** It tells you the expected
  precursor/substrate, which interventions and taxa are mechanistically
  plausible, and whether the intervention branches are chemically coherent with
  the proposed chemistry — a drug-conjugate query hitting cultures supplemented
  with that drug's parent is a far more coherent result than the same tree for an
  unrelated molecule. Say so explicitly.

## Working from the results HTML

Write the script below to disk and run it on the uploaded file. Do not try to
read the HTML directly — it is ~1.5 MB of minified CSS/JS around one data block.

```python
#!/usr/bin/env python3
"""Extract the microbiomeMASST tree + per-file matches from a results HTML page.
Usage: python3 masst_parse.py results.html      (stdlib only)"""
import re, sys, json

def load(path):
    s = open(path, encoding='utf-8', errors='replace').read()
    i = s.find('const root=')
    if i < 0:
        sys.exit('No `const root=` block - not a microbiomeMASST results page?')
    j = s.index('{', i); d = 0; instr = False; esc = False; q = ''
    for k in range(j, len(s)):
        c = s[k]
        if instr:
            if esc: esc = False
            elif c == '\\': esc = True
            elif c == q: instr = False
            continue
        if c in '"\'': instr, q = True, c
        elif c == '{': d += 1
        elif c == '}':
            d -= 1
            if d == 0: return parse(s[j:k+1])
    sys.exit('Unbalanced braces')

def parse(js):                      # JS object literal -> JSON
    out, i, n = [], 0, len(js)
    while i < n:                    # copy quoted strings verbatim
        c = js[i]
        if c == '"':
            m = re.match(r'"(?:[^"\\]|\\.)*"', js[i:]); out.append(m.group(0)); i += m.end()
        else:
            out.append(c); i += 1
    txt = ''.join(out)
    txt = re.sub(r'([{,])\s*([A-Za-z_$][\w$]*)\s*:', r'\1"\2":', txt)   # bare keys
    txt = re.sub(r'([:\[,])\s*(-?)\.(\d)', r'\g<1>\g<2>0.\g<3>', txt)   # .81 -> 0.81
    return json.loads(txt)

def walk(nd, p=()):
    yield nd, p
    for c in nd.get('children') or []:
        yield from walk(c, p + (nd['name'],))

root = load(sys.argv[1])
nodes = list(walk(root))
hit = [(n, p) for n, p in nodes if (n.get('matched_size') or 0) > 0]

print(f"ROOT matched {root.get('matched_size')} / {root['group_size']} files "
      f"({root.get('occurrence_fraction', 0):.4%})")
print(f"{len(nodes)} nodes in tree, {len(hit)} with >=1 match\n")
print(f"{'match':>6} {'total':>7} {'frac':>7}  node <- lineage")
for n, p in sorted(hit, key=lambda x: -x[0]['matched_size']):
    print(f"{n['matched_size']:>6} {n['group_size']:>7} {n.get('occurrence_fraction',0):>7.3f}"
          f"  {n['name']} <- {' > '.join(p[1:])}")

recs = [r for n, _ in hit for k in ('matched_files','matches','data','file_data')
        if isinstance(n.get(k), list) for r in n[k] if isinstance(r, dict)]
if recs:
    ds = {}
    for r in recs:
        u = str(r.get('USI', '')); ds.setdefault(u.split(':')[1] if ':' in u else '?', set()).add(u)
    print(f"\n{len(recs)} file-match records across {len(ds)} datasets:")
    for k, v in sorted(ds.items(), key=lambda x: -len(x[1])):
        print(f"  {k}: {len(v)}")
    cos = [r['Cosine'] for r in recs if 'Cosine' in r]
    mp  = [r['Matching Peaks'] for r in recs if 'Matching Peaks' in r]
    if cos: print(f"cosine {min(cos)}-{max(cos)}")
    if mp:  print(f"matching peaks {min(mp)}-{max(mp)}")
```

### What the HTML gives you that a screenshot doesn't

- **`matched_size` / `group_size` / `occurrence_fraction` per node** — the exact
  numerator and denominator behind every fill color. Report real counts
  ("156/2,979 urine files") instead of describing shades, and never tell the user
  to hover: you already have what hover would show.
- **The complete tree, including zero-match nodes.** The object contains every
  node in the queried metadata tree, so this is *full*-mode data regardless of
  what the rendered view showed. A zero there indicates that no qualifying
  spectral match was observed among the indexed files assigned to that node under
  the selected search parameters — the tree-mode caveat does not apply, and
  control/baseline absences become reportable. Word it that way rather than
  "negative": it is a statement about the index and the search settings, not
  about the biology, and both of the null levels in Part 2 §8 remain live for any
  zero node.
- **Per-file match records**: USI, Cosine, Matching Peaks, Delta Mass. These let
  you check things §5 and §7 of Part 2 ask for and a screenshot can never answer:
  - **How many independent datasets?** Group USIs by accession (MSV…/MTBLS…).
    If most matches sit in one dataset, say so with the numbers — this is the
    single most common way a MASST result looks stronger than it is.
  - **Matched-peak counts.** If every record shares the same low count (e.g. 3),
    discrimination is weak no matter how high the cosine — flag it prominently.
  - **Cosine range and Delta Mass spread**, which bound match quality.
  - The USIs themselves are the pointer back to primary evidence; cite a couple.
- **`NCBI`, `Interventions`, `Community_composition` fields** — the metadata
  behind each node label, useful when a rendered label is truncated.

### HTML-specific cautions

- **Nodes are facets, not disjoint groups.** Files are counted under Origin,
  Organ, and Health phenotype simultaneously, and many nodes carry
  `duplication:"Y"`. Child counts within one facet can sum past the parent, and
  counts across facets describe the same files. Never add matched_size across
  branches to get a total — the root's `matched_size` is the file count.
- **A tiny fraction at the root is expected and means nothing.** The root
  denominator is the whole repository, so occurrence_fraction there is near zero
  for essentially any compound. Interpret fractions at leaf level, comparing
  siblings.
- **Small denominators still saturate.** 1/1 and 200/200 both read as 1.000, and
  now you can see which — say which.
- Counts are files, not subjects: replicates and repeated injections inflate both
  numerator and denominator.

## Reading a screenshot (when no HTML is available)

## Node fill color: what it encodes

**Node fill color represents the occurrence of the queried MS/MS spectrum
within each category, calculated as the number of samples containing a
qualifying spectral match divided by the total number of samples available for
that node.** It is a *proportion*, not a count.

This is readable directly off the static image, so use it — it's usually the
most informative thing in the render after the tree topology itself:

- Fully saturated nodes = a high fraction of that node's samples carried a
  qualifying match. Partially filled / split nodes = an intermediate fraction.
  Nodes rendered in the low-occurrence color = a match exists somewhere in that
  branch, but in a small fraction of that node's samples.
- Compare *siblings*, not distant branches: within one category, relative fill
  is the interpretable signal (e.g. which drug regimen, which isolate, which
  phenotype has proportionally more positive samples).
- **The denominator is invisible in the render, and it dominates.** A node with one sample and
  one match renders identically to a node with 200/200. A pale node over
  hundreds of samples can represent far more independent observations than a
  saturated node over three. Never rank branches by fill alone, and never
  convert fill into a count, a prevalence estimate, or a statistical difference —
  there is no test behind these colors.
- **The actual numerator and denominator are recoverable.** If you have the
  results HTML, they are in it — use the real counts and skip all of the hedging
  in this section. If you only have a screenshot, they sit behind a hover in the
  live interface: tell the user which specific nodes to hover for the claim
  they're making, deliver the full interpretation anyway, and never end by asking
  them to report back.
- Fill is still not prevalence in the population: the denominator is *samples
  deposited in the repository for that node*, shaped by which studies happen to
  exist, not by biology.

### What else is and isn't in a screenshot

Beyond fill color, the rendered tree gives you plain node labels and topology.
Two structural facts about the interface constrain every reading of an image
(neither applies when you have the HTML):

- **Hover-only metadata is genuinely unavailable.** Raw sample/file counts,
  community group, drug/intervention mixture identifiers, incubation time points
  beyond what's written in a label, and culture medium composition sit behind a
  hover interaction. If a claim depends on any of it — especially a
  time-dependent conversion ("formed only after incubation") or an exact number
  of files — state that the image doesn't support it. Don't ask for it.
- **Tree mode is usually not stated, and it changes what a missing branch
  means.** The three modes are *matched* (only nodes with a hit), *matched+other*
  (matched nodes plus their un-matched siblings — the mode that lets absence in
  controls/blanks mean something), and *full* (the whole queried tree). Absent an
  explicit statement, treat a missing branch as uninformative and say so once,
  rather than concluding a true negative.
- **A chemical structure alongside the tree is an added annotation**, not part of
  microbiomeMASST's own output. Treat it as context if present; don't expect it.

Several validation steps from Part 2 (raw MS1/MS2 inspection, matched-fragment
specificity, mass tolerances/cosine threshold) can't be done from a tree image at
all — they need the USI or raw data. Say so plainly rather than glossing over it.

## Step 1. Read the graph structure
Identify, in order:
- **Root node**: the queried molecule (from the user's supplied name; fall back to the root label and Δmass if given). Note whether the name denotes a full structure, a molecular-family placeholder, or an "X-like"/analogue label.
- **Root's branches**: typically **Host** and **Origin**, sometimes directly **Microbial**. Note what sits under each — species (e.g. *Mus musculus*, *Homo sapiens*), organs/biofluids, and whichever intervention categories are present (genetics, diet, drug, prebiotic, viral infection, etc.).
- **Microbial branch**: taxonomic path (phylum → family → genus/species) and whether it's tagged as monoculture, community, or in vitro screening. Note any time point written into the labels (e.g. "T72") and whether a T0 or vehicle arm is rendered.
- **Counts or fill**: from HTML, the matched/total figures per node; from a screenshot, the occurrence-proportion reading above, comparing siblings within each category.
- **Category color coding**: map each shaded category box (Health phenotype, Organs and biofluids, Microbial information, In vitro screening, Biological context/Interventions) back onto the branches you just listed — read it off the specific image rather than assuming a fixed palette. Keep this distinct from *node fill*, which is occurrence, not category.
- **Blank / control nodes**: a node labeled "blank" is ambiguous between a metadata blank (empty field) and an analytical blank. Say which reading you're taking and why, and flag that an analytical blank with nonzero fill is a carryover problem to resolve before publication.

## Step 2. Answer the four framing questions, only as far as the graph supports
- **Occurrence** — directly supported: a matching MS/MS spectrum was found in these deposited files/datasets, at the sample fractions the fill colors encode.
- **Biological context** — supported at the metadata level: which hosts, body sites/biofluids, health states, or interventions the matching files came from.
- **Producer prioritization** — only a hint, not a conclusion: note if the microbial branch includes monocultures, defined communities, or substrate-supplemented cultures (stronger prioritization signal) vs. only complex community samples (weaker). Prioritization is not experimental support for a transformation — that needs temporal dependence or controls, which a single-timepoint tree doesn't show.
- **Structural identity** — the graph alone essentially never settles this, and the user-supplied name doesn't change that. If the name denotes an exact structure, flag that confirming it requires an authentic standard (retention time + MS/MS).

## Step 3. Write the interpretation
Use this structure (skip a subsection only if the graph has nothing relevant to it):

```
## What the graph shows
[1–2 sentences: the queried molecule, and the shape of the tree —
which hosts/origins/microbes it branches into]

## Occurrence
[Which datasets/sample types the spectrum matched in. From HTML: real counts
(matched/total) for the nodes that carry the story, how many independent
datasets the matches span, and whether one dataset dominates. From a screenshot:
what the relative fill proportions across sibling nodes indicate, closing with
which nodes to hover for the underlying numbers]

## Biological context
[Hosts, organs/biofluids, health phenotypes, and interventions
(diet, drug, genetics, prebiotic, infection, etc.) associated with matches —
described as associations found in the searched files, not effects or prevalence.
Note whether the interventions are chemically coherent with the queried molecule]

## Microbial / producer signal
[What's in the Microbial branch, and whether it supports producer
prioritization (monoculture/substrate-supplemented) or is weaker
(community-only, no cultured isolate)]

## Match quality  [HTML only — omit for screenshots]
[Cosine range, matched-peak counts, delta-mass spread across the file records,
and what they imply about discrimination. Uniformly low matched-peak counts
undercut a high cosine and must be said plainly. Cite one or two USIs as the
pointer back to primary evidence]

## What this input can't settle
[Only if relevant. Screenshot: raw counts and denominators, incubation time
points, culture composition, community group, tree mode — name the limit and
which nodes to hover. HTML: the things no MASST output settles — structure,
causality, and anything needing the raw MS1/MS2 behind the USIs]

## Caveats
[3–6 bullets, selected from "Caveats to check every time" below —
pick the ones that actually apply to this graph rather than listing all of them by rote]

## Evidence level reached
[Name the level reached and say what it means in a few words, so the reader
never has to hold the hierarchy in their head. Never cite a bare number —
"level 2" is opaque; "metadata association — the match is linked to documented
host and intervention context" is not. Sketch the ladder compactly (spectral
occurrence → metadata association → experimental support for a producer or
transformation → structural validation → biological or mechanistic validation),
mark where this result sits, and name the one experiment that moves it up.
Structural confidence and biological attribution are separate axes: place the
result on each, since strong evidence on one does not establish the other]
```

Keep it tight. A dense page beats a thorough three — every sentence should
either state a finding or constrain one.

## Step 4. Caveats to check every time
Pull from these (Part 2, §5–10 have the full reasoning):
- **Fill proportion ≠ prevalence, and ≠ counts.** The denominator is deposited samples per node, and often tiny — a saturated node may be 2/2. Get the real numerator and denominator from the HTML (or hover) before treating any fill difference as meaningful.
- **Counts are facets of one file set, not disjoint groups (HTML).** The same file is counted under Origin, Organ, and Health phenotype; nodes flagged `duplication:"Y"` recur across branches. Never sum matched_size across branches — the root's matched_size is the file total.
- **File counts reflect study representation.** Repositories have uneven numbers of studies/samples/organisms per category; branch size is a fact about deposition, not biology.
- **Replicates can inflate both numerator and denominator.** Multiple files may be biological replicates, technical replicates, or repeated injections of the same sample.
- **Reproducibility across datasets matters.** A match concentrated in one study is weaker evidence than the same match appearing across independent datasets. From HTML this is checkable — group the USIs by accession and report the split; from a screenshot it isn't, and you should say so.
- **Structural analogues and isomers are a live possibility**, especially for anything not confirmed against an authentic standard — report at the appropriate hedge level ("X-like spectrum," "candidate positional isomer") rather than asserting the queried name as established.
- **Absence of a branch is not absence of the compound.** A null result only means no qualifying match was found under the search conditions used — see Part 2 §8 for the reasons (wrong polarity/adduct, in-source fragment, unindexed dataset, etc.).
- **Distinguish the two kinds of null (Part 2 §8).** Either no qualifying spectral match exists in the searched index, or a match exists but its biological metadata is missing or too coarse to place on the tree. These are not equivalent, and an unrendered branch can be either — never report a metadata gap as a chemical negative.
- **Tree mode affects what "missing" means — for screenshots.** In *matched* mode, an absent branch may simply not have been rendered; only *matched+other* or *full* mode lets an absence be reported at all. The HTML carries the full tree, so a zero there indicates that no qualifying spectral match was observed among the indexed files assigned to that node under the selected search parameters. Never restate that as biological absence, an undetected compound, or a true negative — the compound may be present but unfragmented, ionized differently, below threshold, or in files whose metadata never placed them on that node (Part 2 §8). Either way, a control/baseline contrast strengthens a transformation hypothesis but does not by itself establish microbial causality.
- **Time-course claims need both arms.** Don't state or imply a time-dependent conversion (e.g. "detected only after incubation") from labels alone. From HTML, check whether a T0/vehicle node exists and what its matched_size is — that turns the claim into a real comparison.
- **The graph is not the primary evidence.** It's an interpretation layer; claims of real biological significance should point back to the original spectrum/dataset (USI) and, ideally, orthogonal validation.

## Step 5. Match the claim to the evidence level
Never let the write-up's confidence exceed what Steps 2–3 actually established. If the user (or their manuscript draft) is about to state something like "compound X is produced by microbe Y in the gut," check it against the evidence hierarchy in Part 2 §10 — that specific claim needs experimental support for a producer or transformation at minimum (temporal dependence, controls, substrate addition, monoculture or enzyme assay), ideally structural validation, not just a spectral match in a fecal sample.

Two framing points from §10 that are easy to lose in a write-up: the levels are complementary dimensions, not a mandatory sequence, so don't imply a result must climb them in order; and structural confidence and biological attribution are evaluated separately, so a well-validated structure with weak attribution (or the reverse) should be reported as exactly that rather than averaged into one confidence statement.

## Tone
Write like a co-author doing a sanity check, not a disclaimer generator. Lead
with what the data does show, then attach caveats that are specific to *this*
graph rather than a generic boilerplate list — a graph with only one dataset
behind it needs the reproducibility caveat much more than a graph with matches
across a dozen independent studies. End on the finding or the next experiment,
never on a question.

---

# Part 2 — Practical guide for querying and interpreting microbiomeMASST

## Scope

We developed an AI-assisted interpretation guide, implemented as a Claude Agent Skill, to standardize how microbiomeMASST searches are performed, interpreted and reported. The skill guides users through query-spectrum selection, search parameterization, spectral-match assessment, provenance inspection, biological-context interpretation and evidence evaluation, while surfacing the relevant limitations at each stage. Several of these limitations – including incomplete MS/MS acquisition, dependence on ionization and acquisition conditions, ambiguity among structural analogues and the requirement for orthogonal structural validation – are inherent to untargeted LC–MS/MS workflows based on spectral matching and are not specific to microbiomeMASST. Others, including metadata completeness, uneven repository representation, interpretation of graph relationships and non-independence of repository observations, arise specifically from repository-scale contextualization. The complete decision logic is provided in SKILL.md, making the interpretation framework transparent and inspectable independently of the language-model implementation.

## 1. Define the question before running the search

The appropriate interpretation depends on the question being asked:

- **Occurrence:** In which deposited files has a related MS/MS spectrum been acquired?
- **Biological context:** In which hosts, body sites, health states, interventions, or microbial systems does it occur?
- **Producer prioritization:** Is it detected in monocultures, communities, substrate-supplemented cultures, or enzyme assays?
- **Structural identity:** Is the match sufficient to support a molecular family, candidate structure, or exact compound?

MicrobiomeMASST can directly address the first two questions and can help prioritize the third. It cannot establish exact structure, microbial causality, or biological function without follow-up experiments.

## 2. Confirm that a queryable MS/MS spectrum exists

A compound name, formula, SMILES string, or precursor mass alone is not sufficient for a MASST search. The user should provide a USI linked to an acquired spectrum or the *m/z* values and matching intensities, a spectrum from an authentic standard, a representative library spectrum, or an experimentally acquired spectrum of the unknown feature. Confirm precursor charge, polarity, proposed adduct, precursor *m/z*, and spectrum quality before searching. Low-information spectra with only a few nonspecific fragments produce less discriminating searches.

## 3. Select the most appropriate reference spectrum

Prefer a spectrum acquired under conditions resembling those expected in the repository data. Consider polarity, adduct form, possible in-source fragments/neutral losses/multimers, and collision energy appropriateness. Use StructureMASST to search multiple spectra/ion forms rather than one expert-selected spectrum when possible.

## 4. Run an initial search using documented thresholds

Record precursor-mass tolerance, fragment-mass tolerance, minimum cosine similarity, minimum matched fragment ions, library/corpus searched, and search date/index version. Manuscript settings can anchor the primary analysis, but no single threshold combination is universally optimal – lower thresholds help find candidate molecular-family relationships; exact-structure claims need more stringent inspection.

**Core interpretation rule:** Do not convert spectral similarity into exact structural identity, metadata association into causality, file counts into prevalence, or a non-match into biological absence.

## 5. Evaluate the quality of each spectral match

Do not rely on cosine similarity alone. Examine precursor-mass agreement, number of matched fragments, whether matches span the spectrum or cluster in one region, whether ions are chemically diagnostic vs. generic low-mass fragments, relative intensities, possible chimeric isolation, ion form/charge, and instrument/collision energy. Four generic fragments provide weak discrimination; four diagnostic fragments can be informative – specificity matters as much as count.

## 6. Inspect the original data and provenance

Follow the USI/dataset link to the original spectrum and raw data. Check for a chromatographic MS1 peak, whether MS/MS was acquired near the apex, chimeric isolation, plausible isotope/adduct patterns, and whether sample name/experimental group/dataset publication support the metadata assignment. The graph is an interpretation layer – the original spectrum and file remain the primary evidence.

## 7. Interpret the biological-context graph conservatively

Ask whether matches reproduce across independent studies or concentrate in one dataset, whether files represent biological replicates/technical replicates/repeated injections, whether host/intervention/exposure is documented, whether the sample type is biologically plausible, whether the molecule co-occurs with a plausible parent/precursor/substrate, and whether the association is driven by overrepresented studies or organisms. File counts are NOT prevalence estimates – repositories have uneven numbers of studies, samples, organisms, interventions, and methods.

## 8. Diagnose a search returning no matches

A null result means only that no qualifying match was found under the selected search conditions and should not be interpreted as biological absence. Possible explanations include lack of a relevant deposited study, abundance below the MS/MS acquisition threshold, failure to trigger fragmentation, a different polarity or adduct, in-source fragmentation or multimer formation, incompatible collision energy, an insufficiently informative reference spectrum, structural isomerism or analogy, overly stringent search parameters, or absence of the relevant dataset from the current index.

Importantly, a null microbiomeMASST result can arise at two distinct levels: **(i)** no qualifying spectral match is found in the searched index, or **(ii)** a spectral match exists but the corresponding biological metadata is unavailable or insufficiently resolved for contextualization. These outcomes should not be interpreted equivalently.

When appropriate, follow-up can include searching the broader FASST corpus, testing alternative experimentally plausible spectra or ion forms, inspecting candidate datasets at the MS1 level, acquiring a higher-quality reference spectrum, or generating the spectrum experimentally. Relaxed search thresholds may be useful for exploratory hypothesis generation but should not be used to present low-confidence matches as confirmed molecular identifications.

## 9. Treat matches to structural analogues as hypotheses

For closely related analogues sharing fragmentation: can MS/MS distinguish the modification site? Could another positional/geometric/stereochemical isomer explain the same fragments? Are matched ions specific to the whole molecule or just a shared scaffold? Does retention time support the proposed polarity/structure? Is biological context compatible? Are standards available for competing candidates? Where ambiguous, report at the appropriate level (e.g., "eicosatrienoic-acid-like spectrum," "bile acid amidate," "candidate positional isomer") rather than an exact structure.

## 10. Match the evidence to the intended conclusion

Evidence can accumulate along complementary analytical and biological dimensions; the categories below should not be interpreted as a mandatory linear sequence.

Evidence hierarchy (do not claim beyond the highest level actually obtained):

- **Spectral occurrence** – a qualifying MS/MS match in a deposited file.
- **Metadata association** – linkage to documented host, sample, intervention or microbial context.
- **Experimental support for a producer or transformation** — temporal dependence, controls, substrate addition, monocultures, defined communities or enzyme assays.
- **Structural validation** – authentic-standard MS/MS and retention time, ion mobility, synthesis, NMR or another structure-specific method.
- **Biological or mechanistic validation** – genetic, biochemical, cellular, animal or prospective human evidence.

Structural confidence and biological attribution should be evaluated separately; strong evidence for one does not automatically establish the other.

## Interface-specific interpretation rules (from the microbiomeMASST supplementary figures)

- Search input is a USI or manual spectrum entry (precursor *m/z* + charge); search parameters (precursor/fragment tolerance, cosine threshold, matched peaks) are user-set at query time.
- After searching, three tree visualization modes are available: **matched** (only nodes with a spectral match), **matched+other** (matched nodes plus related sibling nodes without a match; for example, this mode can reveal absence at baseline or in appropriate culture controls, which can strengthen evidence for culture-dependent formation. Such contrasts strengthen a transformation hypothesis but do not alone establish microbial causality), and **full** (the entire queried metadata tree). Which mode produced a given image changes what "no branch here" means – in "matched" mode it means nothing was searched/shown, not that it was searched and absent.
- The rendered tree image shows node labels, topology, and **node fill color, which encodes occurrence**: the number of samples containing a qualifying spectral match divided by the total number of samples available for that node. Fill is therefore a proportion that *is* readable from a static image – but the denominator is not, so fill can never be converted back into counts or prevalence.
- **Hovering over a node in the live interface reveals additional metadata**: community group, drug/intervention mixture identifier, microbial incubation time (e.g., 0 h vs. 72 h), the list of compounds supplemented in the culture medium (drugs, polyamines, bile acids, etc.), community membership, and raw sample/file counts (the numerator and denominator behind the node's fill). None of this is recoverable from a static screenshot – so an interpretation should state the limit and direct the user to hover the relevant nodes themselves, rather than requesting the data and waiting.
- Incubation time points matter for producer/pathway claims specifically: a compound appearing only at 72 h and not 0 h in a monoculture is much stronger evidence of microbial conversion than a single timepoint with no comparison.

## 11. Report the search reproducibly

For every result, report: query USI/library accession; proposed ion form; search index and date; mass tolerances; cosine threshold; minimum matched peaks; number of matching spectra/files/independent datasets; key metadata associations; raw-data inspection performed; competing structural explanations; orthogonal validation completed; and whether the result is a spectral match, structural hypothesis, or confirmed identification.

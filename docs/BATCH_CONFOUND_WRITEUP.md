# Zebrafish drug-screen analysis: a batch-confound question for a fresh perspective

## What this document is

I'm analyzing a zebrafish behavioral drug screen. Over a long back-and-forth with
another Claude instance, we worked through a statistical problem about batch
(plate) effects and how — or whether — drug effects can be compared across
plates. I want an independent read on whether the reasoning and the conclusion
are sound. Below is the experiment, the model, the questions I asked, and the
explanations I was given. Please push back wherever you disagree.

## The experiment and data

The assay measures larval zebrafish movement (a "pixel difference" signal) over
time in multi-well plates. Each well holds one animal. We bin movement to
one-second sums per well and model the time course.

The screen ("initial screen") spans three plates. Critically, **each plate
contains a different set of drugs**: Plate 1 holds drugs labeled 01A02, 01A03,
... ; Plate 2 holds 02A02, 02A03, ... ; Plate 3 holds 03A02, etc. There are
roughly 82 distinct conditions per plate. Each plate has 2–3 replicates, and the
replicates are not always identical layouts (e.g., round-well vs square-well
plates, recorded on different days). The total is about 658 animals across the
three plates and their replicates.

The one thing shared across every plate is the control condition: WT+DMSO
(wild-type genotype, DMSO vehicle). Every plate has DMSO control wells. Every
non-control drug, however, appears on exactly one plate and one plate only —
there is no replication of any given drug across plates or batches.

A separate, smaller "confirmation screen" uses a balanced design (a full
genotype × drug factorial with only a handful of conditions) and is analyzed
separately; it is not the subject of this question.

## The model

We fit a Tweedie GAM (mgcv::bam in R) with an AR(1) correction for temporal
autocorrelation within each animal's time series. The model has a shared
baseline time smooth, a per-condition smooth that lets each condition have its
own movement trajectory over time, and a per-animal random effect. Each drug's
effect is reported as a contrast against WT+DMSO on the same plate, expressed as
a log2 fold-change, with a p-value, and we apply multiple-testing correction.

(We spent a lot of effort separately getting this to run at scale — the
per-animal term had to become a random intercept rather than a per-animal smooth
to be computationally tractable — but that is a tractability matter, not the
subject of this question.)

## The core problem I raised

Because each drug lives on exactly one plate, I asked whether it even makes sense
to model all plates together, and whether per-plate results can be fairly plotted
or compared on one figure, given that the p-values and effect sizes would come
from plate-specific information.

The explanation I was given, which I want checked:

The only cross-plate link is the shared DMSO control (plus, in a combined model,
a shared error scale and shared baseline smooth). Since no drug is replicated
across plates, **drug effect and plate/batch effect are perfectly confounded
(perfectly aliased) for every non-control drug.** No statistical model — whether
one combined model or three separate per-plate models — can decompose "this
drug's effect" from "this plate's batch effect," because the information needed
to separate them (the same drug measured in more than one batch) does not exist
in the data. Within a single plate, drug-to-drug comparisons are clean because
all those drugs share the same batch, so the batch component cancels. Across
plates, every comparison carries an unquantifiable batch component.

## My batch-correction questions and the answers I got

I then asked whether batch effects could be corrected, either before or after
modeling, and whether it's better to correct from the raw data or from the model
estimates. The explanation I was given:

First, a hard ceiling: no batch-correction method (ComBat, RUV, limma-style
removeBatchEffect, quantile normalization, etc.) can recover something perfectly
confounded. All such methods require either the same condition measured across
multiple batches (which we don't have) or an assumption that lets the batch
effect be estimated from something treated as batch-invariant. So correction here
is only possible via an assumption, and can only remove the part of the batch
effect that the assumption captures. The drug-by-batch interaction (a batch
affecting drugs differently than it affects the control) is irrecoverable by any
method; only a different experimental design (drugs replicated across batches)
could break it.

Then a ladder of what can validly be corrected:

Level 0 is what we already do: express each drug as a log2 fold-change versus its
own plate's DMSO. Because DMSO is the one thing replicated across plates, this
removes the additive, constant per-plate offset (the "this plate ran high or low
overall" component). It needs no extra assumption beyond "DMSO means the same
thing on every plate."

Level 1 is a robust plate normalizer: assume most drugs on a plate are inactive,
so the bulk distribution of all wells on a plate (its robust median and spread)
is itself a control estimate, richer than DMSO alone. Normalize each plate by its
median and MAD. This removes additive and scale batch effects. Its cost is the
"most drugs inactive" assumption, and it is still blind to drug-by-batch
interaction. If a plate happens to contain many genuinely active drugs, this
normalizer would wrongly subtract real signal.

On the before-versus-after question, the answer was that correcting on the model
estimates (after) is better than correcting the raw data (before), for three
reasons: (1) the raw signal is Tweedie-distributed with a specific mean-variance
link, and shifting or scaling raw counts breaks that structure the model relies
on; (2) pre-correcting then feeding "cleaned" data into the model makes inference
anti-conservative — the model thinks it has full degrees of freedom when some
were already spent removing batch, so p-values come out too small (the standard
warning, e.g. Nygaard et al. 2016 on ComBat); and (3) the natural place for the
correction is the contrast scale, which is exactly where DMSO-referencing already
operates. So the preferred order is: within-model control referencing (best),
then post-hoc correction of the model's contrasts if going past Level 0, and
pre-correcting raw data is the worst option.

## Where we landed

I proposed: model each plate separately, batch-correct the estimates, then plot
them together "fairly." The pushback I got was that "fairly" was still too
strong, because even after per-plate modeling plus Level-1 correction, two things
survive: the drug-by-batch interaction (irrecoverable), and the load-bearing
"most drugs inactive" assumption (which subtracts real signal if violated). So
that pipeline is fairer than raw cross-plate contrasts and is a standard
screening approach, but it is not a clean quantitative cross-plate potency claim
— only a screening/hit-picking overview, and the axis should be framed as
"normalized deviation from plate control," not absolute drug potency.

I then asked whether the only truly clean approach is to model and plot the
plates separately. The answer was yes: model each plate separately, plot each
plate separately, and compare drugs only within a plate against that plate's
DMSO. This is the only approach with zero unverifiable assumptions, because it
simply never asks the cross-plate question that the confound corrupts. The
tradeoff is that you give up a single screen-wide hit ranking — you get three
independent per-plate screens — and cross-plate hit comparison would have to wait
for a follow-up experiment designed with cross-batch replication. It was also
argued that this fits the scientific reality of the design: the experiment was
structured as per-plate drug panels sharing a DMSO control, so analyzing per
plate is analyzing at the resolution the design actually supports, not a
workaround. Multiple-testing correction would then be applied within each plate.

## What I want from a fresh perspective

1. Is the perfect-confounding claim correct for this design, and is it really true
   that no single-experiment analysis can separate drug from batch here?
2. Is the before-versus-after batch-correction reasoning sound (especially the
   anti-conservative-inference warning about pre-correcting raw data)?
3. Is "model and plot each plate separately, compare only within plate" actually
   the most defensible conclusion, or is there a legitimate middle path I'm being
   talked out of (e.g., a mixed model with a plate random effect, an explicit
   plate-level term, empirical-Bayes shrinkage across plates, or treating DMSO
   wells as anchoring negative controls in a way that buys more than I think)?
4. Is there any approach that recovers more than per-plate analysis without
   relying on the "most drugs inactive" assumption?
5. Anything important that was missed or overstated.

##Basics of qPCR

Agenda:

1. Absolute quantitation + qPCR principles.
2. STD curve protocol
3. Analyze data from 18S quality check qPCR run.


###qPCR concepts:

1. qPCR is just PCR.
- why use PCR? Physically impossible to count the number of RNA transcripts in a cell.
- We need to amplify our gene of interest.

2. Why is this "difficult"?
- need to ensure that amplification is consistent throughout cycles (no biochemical turbulence).
- need to ensure we are efficient across genes. Amplification bias is rampant!!

###PCR Reaction Protocol:
- template (DNA/cDNA)
- primers: initiate the reaction
- enzyme (*Taq* polymerase)
- dNTPs
- Buffer
- Mg2+
-H2O

If doing qPCR, need a probe for quantification. We use SYBR Green.

Cycling conditions:
95C (initial melt) - 2-3 mins
95C (denature) - 15 secs
60C* (annealing) - 15-30 secs
72C (elongation) - 1 min/kb
72C (final elongation) 5-10 mins

*Note, this temp varies on the melting temp (Tm) of your primers.

###Quantitation in qPCR

**SYBR Green for relative quantification.**

SYBR Green binds directly to DNA, allowing us to measure fluorescence at the end of each cycle.

If we plot cycle # on the X-axis and fluorescence on the Y-axis, we expect to see an exponential growth curve reflecting the amplification of our target gene product through the PCR cycles.

If we plot the curves from two reactions, our target and an internal standard, we can calculate the difference in cycle # to reach a given fluorescence value.
For example, if it takes sample A 19 cycles to reach the same fluorescence as it takes for Sample B after 15 cycles, we know the template (or gene) was at higher abundance in Sample B

**Probe-based quantitation**
Probe is a single-stranded compliment to our gene target.
When the probe binds, a fluorophore is released.

Pros vs. Cons of SYBR Green vs Probes:

|        |  SYBR        |        Probe-based     |
|--------|--------------|------------------------|
|PROS    |Universal     | Specific/can multiplex |
|        |Melt curves   | absolute quantitation  |
|CONS    |Not specific  |Expensive/time consuming|


###Melt-curve:

And the end of PCR, run a "melt curve"

DNA will be bound while double-stranded (at 65C)

We increase the temperature step-wise to 95C, at which the strands will denature. 

At some temp along the curve, our PCR products should all suddenly melt. 
- the curve will appear as a sigmoidal curve.
- instead plot ∆F, which will give you a peak at the melting temp.


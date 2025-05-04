= Wave

The energy industry uses highly specialised software --- written, refined, and perfected over decades. Little scripts are no match against their sophistication. Instead, the idea is to not feel challenged when they are either inaccessible or too time-consuming to reproduce results using them, before offering feedback on the work done.

Recently while reviewing a young colleague's work I realised my handicap was the tooling. 

To avoid getting into the annals of specific software, my request to engineers --- I assist and guide typically --- is simple: a print / screenshot of output from analyses, which I can look at, and ask reasonable questions. I do some simple calculations on the side, if required, to assist in the review.

In my role, I consciously choose to guide and lead, but not necessarily roll-up sleeves and work it myself. This has its benefits:

- Working the state-of-the-art software myself is inefficient
- My time is better spent reviewing and troubleshooting
- Given the above, avoid locking-in expensive licenses
- Good for knowledge transfer, and better collaboration
- I have time to develop a strategy than technique
- Avoid confirmation bias, verify independently
- Use simpler techniques, where practicable

But some things are beyond the realm of the back-of-the-envelope, and demand a bit of rigour. If iteration is involved, then it quickly gets out of hand. In such cases, if I have to write some code on the side, then the attempt is to get some ROI, i.e. not only to get reasonably approximate answers, but also some insights as well.

Back to the problem, I needed to generate hydrodynamic forces on a slender caisson to review the force profile in segments of the water column I was interested-in --- for different hydrodynamic coefficients. Previously, I'd do this by simply building a stick model and run a wave through it. In this case, without the software pre-installed on my work machine, the only option was to write some code from first principles, and calculate hydrodynamic forces in close approximation to what an expensive software suite would produce.

The conventional way to generate hydrodynamic forces is by employing Morison's equation:

$ F = F_"drag" + F_"inertia" $

$ F_"drag" = 1/2 rho C_d D U |U| $

$ F_"inertia" = C_m rho V (partial U / partial t) $

where,

- F -- local action vector per unit length normal to the tube
- $C_d$, $C_m$ -- drag and inertia hydrodynamic coefficients respectively
- $rho$ -- mass density of seawater
- D -- effective width of tube (incl. marine fouling)
- V -- displaced volume of tube per unit length respectively
- U, |U| -- component of local water particle velocity vector normal to the tube, and its absolute value respectively
- $ partial U / partial t$ -- component of the local water particle acceleration vector normal to the tube

The orbital velocity $U$ can be taken as the sum of $U_"wave"$ and $U_"current"$ (from the profile); where the angular frequency, $omega = (2 pi) / T$, and where $H$ and $T$ are wave height and associated period respectively. 

$ U_"wave" = H omega (cosh(kappa z))(2 sinh(kappa d)) cos(omega t) $

The wave orbital acceleration can be expressed as:

$ (partial U / partial t) = -H omega^2 (cosh(kappa z) / 2 sinh(kappa d)) sin(omega t) $ 

The simplification here is that it does not account for marine-fouling, which would further expand the problem and the resulting code into an array of tube diameters. For a sense-check though, this seems sufficient.

#figure(
  image("/img/hydforce.svg", width: 100%),
  caption: [Hydrodynamic force profile in seawater column.]
) <hydforce>

The code for the plot is as follows. 

#let wave = read("/src/wave.py")
#{linebreak();raw(wave, lang: "python")}
#v(1em)

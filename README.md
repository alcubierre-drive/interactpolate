# Coulomb Interactions in real-space, but fast

This small library helps with Coulomb interactions that don't follow a simple
(as in computationally simple/fast) functional form, but are expensive to
calculate. Currently, the only interaction profile that's implemented is the
dual-gated Coulomb interaction [cf. Phys. Rev. B 86, 115447 (2012)], but with an
$`r=0`$ regularization a la Ohno. Other profiles could be implemented trivially in
the future.

The long-ranged part of the interaction is given as $V_\mathrm{long}(r)$:
```math
V_\mathrm{long}(r) = 4V_0 \sum_{k=0}^{\infty} K_0\left[ 2(k+1)\, \pi \frac{r}{\xi} \right] \,,
```
where $`\xi`$ is the gate distance to both gates, s.t. the two gates have a
distance $`2\xi`$ with the 2D sample in between. $`K_0`$ is a Bessel function of
the second kind, and $`V_0`$ is an energy scale.

The short-ranged part $`V_\mathrm{short}(r)`$ must follow
```math
V_\mathrm{short}(r) = \frac{U}{\sqrt{1+r^2/a^2}}\,,
```
where $`U`$ is the on-site interaction strength (Hubbard-$`U`$) and $`a`$ the
"Ohno-parameter".

We make a smooth version of both profiles (i.e. a long- and short-ranged
screened Coulomb interaction) by plugging the one into the other (one can check
that the limits work out...). We arrive at
```math
V(r) = 4V_0\sum_{k=0}^{\infty} K_0\left[ 2(k+1)\, \pi \frac{\sqrt{r^2+a^2}}{\xi} \right] \,,
```
with the following constraints on $`V_0`$ and $`a`$:
```math
V_0 = \frac{\alpha}{\epsilon\,\xi} \,,\qquad
a = \frac{\alpha}{\epsilon\,U} \,.
```
The fine-structure constant $`\alpha`$ and the dielectric constant $`\epsilon`$
determine the interaction strength s.t. it follwos $`\alpha/\epsilon r`$ in
intermediate distances. Note that $`\epsilon`$ should not carry any units, and
$`\alpha = 14.40\,\mathrm{eV}\text{\AA}`$ in SI units (inserted some $`\hbar c`$
to get those units).

# Usage

see test/use.c and test/Makefile for how to use this library. It's fast as long
as use have one handle (roughly as fast as evaluating exp(..)/sqrt(..)),
creation of the handle is what's slow.

`

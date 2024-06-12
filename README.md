# Coulomb Interactions in real-space, but fast

This small library helps with Coulomb interactions that don't follow a simple
(as in computationally simple/fast) functional form, but are expensive to
calculate. Currently, the only interaction profile that's implemented is the
dual-gated Coulomb interaction [cf. Phys. Rev. B 86, 115447 (2012)], but with an
r=0 regularization a la Ohno. Other profiles could be implemented trivially in
the future.

# Usage
see test/use.c and test/Makefile for how to use this library. It's fast as long
as use have one handle (roughly as fast as evaluating exp(..)/sqrt(..)),
creation of the handle is what's slow.

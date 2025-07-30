---
title: 'Skyrmions3D: A Julia package to create and visualise 3D Skyrmions in the Skyrme model'
tags:
  - Julia
  - skyrmions
  - topological solitons
authors:
  - name: Chris Halcrow
    corresponding: true
    orcid: 0000-0002-0246-8075
    affiliation: 1
  - name: Linden Disney-Hogg
    orcid: 0000-0002-6597-2463
    affiliation: 2
  - name: Cas G. Chaudhuri
    orcird: 0009-0001-4281-2766
    affiliation: 2
affiliations:
 - name: University of Edinburgh, United Kingdom
   index: 1
 - name: University of Leeds, United Kingdom
   index: 2
date: 30 July 2025
bibliography: paper.bib
---

# Summary

A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.

Skyrmions were introduced in the 1960s as topological solitons which serve as low-energy models of nuclei. Mathematically they are defined as smooth maps $U: \mathbb{R}^3 \to \mathrm{SU}(2)$ which minimise a certain integral (the Skyrme energy) subject to boundary conditions which topologically stabilise the field $U$ so as to prevent it just being the vacuum solution $U=1$. Analytically writing down skyrmions fields is difficult and accordingly many approximation of skyrmions fields, using for example rational maps of the ADHM data of instantons, have been proposed. 

`Skyrmions3D.jl` allows for the creation, transformation, and visualisation of skyrmions generated via various approximations already provided in the package. Moreover, the implementation is written flexibly so as to allow users to investigate their own novel skyrmion approximations. Underlying algorithms are state of the art, and written so as to exploit efficiencies of Julia. 

# Statement of need

A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work.

# Background

A bit of background about skyrmions, their mathematics, and existing work. 

A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.

Note that the paper begins by a metadata section (the enclosing — lines are mandatory) and ends with a References heading, and the references are built automatically from the content in the .bib file. You should enter in-text citations in the paper body following correct Markdown citation syntax. Also note that the references include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.

Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.

# Acknowledgements

Acknowledgement of any financial support.

The research of LDH is supported by a UK Engineering and Physical Sciences Research Council (EPSRC) doctoral prize fellowship.

# References

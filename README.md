# Sky Survey
Given a set of test and equipment parameters, this generates a sky survey - an orientation trajectory that angles a camera to points of interest in the sky. 

A trajectory is a set of timestamped coordinates. In the equatorial coordinate system, this is a set of 1×3 vectors of the form: (time, right ascension, declination). Given multiple points of interest, we can generate different trajectories to pass through them depending on what we want to prioritize.

## Overview

Here’s a high-level view of the function:

<div style="text-align: center;">
    <img src="images/sky-survey-over.svg">
</div>

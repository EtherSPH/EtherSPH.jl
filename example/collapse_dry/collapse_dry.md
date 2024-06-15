[toc]

# Collapse Dry Problem

This case is from [Fluid mechanics and the SPH method: theory and applications](https://academic.oup.com/book/8904), on page 484, by Damien Violeau. In other papers, the case can also be called `Dam break problem`.

<center>
<image src="image/collapse_dry.png">
</br>
<b>fig. collapse onto a dry bottom</b>
</center>

It's also included by [SmoothedParticles.jl](https://github.com/OndrejKincl/SmoothedParticles.jl).

# Different Equation Model Comparison

## Same: treat wall particles as fixed water particles

<center>
<image src="image/collapse_dry_same.png" width=70%>
</br>
<b>fig. collapse onto a dry bottom - treat wall as the same</b>
</center>

<center>
<image src="image/collapse_dry_same_cmap.png" width=70%>
</br>
<b>fig. collapse onto a dry bottom - treat wall as the same</b>
</center>

> note 1: Actualluy, i find this method will meet numerical unstability when `dr` gets smaller. Setting `dr` as 0.01 will more easily meet negative pressure problem. -- 2024.06.15

> note 2: Wall particles provide a pressure foce to avoid fluid particles penertating them. This will create a blank gap between fluid particles and wall particles, which is not the case in real life. -- 2024.06.15

> note 3: a few water particles will hang on the wall. -- 2024.06.15
# Image Metrics

A general discussion about measuring discrepancy between images with a specific
focus on optical interferometry is developed in[^Gomes2016]. A proper image
metric shall provide a quantitative score suitable to order restored images.
This score must reflect the perception by an expert of the fidelity of the
images with a given ground truth or reference image. The score shall be
insensitive to changes that may be considered as irrelevant for the considered
context. For example, image quality indices widely used in the signal
processing community are insensitive to an affine transform of the pixel
values. For optical interferometry, the pixel size of a restored image is, to
some extend, a free parameter and may thus be different from the pixel size of
the ground truth image. Moreover, when only power spectrum and phase closure
data (`vis2data` and `t3phi` in `OI_VIS2` and `OI_T3` data-blocks of OI-FITS
files) are used to reconstruct an image, the position of the object in the
field of view is not constrained by the data and the image metric should not
depend on this position.

```math
\Dist(\Vx,\Vy) = \min_{α,β} \Vert \alpha\,\MR_{\boldsymbol{θ}}
\cdot\Vx + \beta\,\One - \Vy\Vert
```

where ``\Vert\ldots\Vert`` is some norm, ``\One`` is an image of the same size
as ``\Vy`` but filled with ones, ``\MR_{\boldsymbol{θ}}`` is a linear operator
which implements resampling with a given magnification, translation, and
blurring. Here ``\alpha``, ``\beta``, and ``\boldsymbol{θ}`` are nuisance
parameters to reduce the mismatch between the images.

The score may be defined by normalizing the distance:

```math
\Score(\Vx)
= \frac{\Dist(\Vx,\Vy)}{\Dist(\Zero,\Vy)}
= \frac{
  \min_{\alpha,\beta} \Vert \alpha\,\MR\cdot\Vx + \beta\,\One - \Vy\Vert
}{
  \min_{\beta} \Vert \beta\,\One - \Vy \Vert
}
```

where ``\Zero`` is an image of the same size as ``\Vy`` but filled with zeros.

[^Gomes2016]:
    > N. Gomes, P. J. V. Garcia & É. Thiébaut, *Assessing the quality of
    > restored images in optical long-baseline interferometry* in Monthly
    > Notices of the Royal Astronomical Society, vol. **465**, pp. 3823-3839
    > (2016).

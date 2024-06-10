# General image metrics

A general discussion about measuring discrepancy between images with a specific focus on
optical interferometry is developed in[^Gomes2016]. A proper image metric shall provide a
quantitative score suitable to order restored images. This score must reflect the
perception by an expert of the fidelity of the images with a given ground truth or
reference image. The score shall be insensitive to changes that may be considered as
irrelevant for the considered context. For example, image quality indices widely used in
the signal processing community are insensitive to an affine transform of the pixel
values. For optical interferometry, the pixel size and the field of view of a restored
image are, to some extend, free parameters and may thus be different from the pixel size
of the ground truth image. Moreover, when only power spectrum and phase closure data
(`vis2data` and `t3phi` in `OI_VIS2` and `OI_T3` data-blocks of OI-FITS files) are used to
reconstruct an image, the position of the object in the field of view is not constrained
by the data and the image metric should not depend on this position.

Accounting for all these remarks, a possible definition of the distance between ``\Vy``, a
reference image, and ``\Vx``, a reconstructed image, is given by:

```math
\begin{align}
\Dist(\Vx,\Vy) = \min_{α,β,\Vt} \Bigl\{
  &\sum_{i \in |\MR_{\boldsymbol{θ}}\cdot\Vx| \cap |\Vy|}
  d\!\left(α\,(\MR_{\boldsymbol{θ}}\cdot\Vx)_i + β, y_i\right)
  \notag\\
  &+ \sum_{i \in |\MR_{\boldsymbol{θ}}\cdot\Vx| \backslash
  (|\MR_{\boldsymbol{θ}}\cdot\Vx| \cap |\Vy|)}
  d\!\left(α\,(\MR_{\boldsymbol{θ}}\cdot\Vx)_i + β, v_{\mathrm{out}}\right)
  \notag\\
  &+ \sum_{i \in |\Vy| \backslash
  (|\MR_{\boldsymbol{θ}}\cdot\Vx| \cap |\Vy|)}
  d\!\left(v_{\mathrm{out}},y_i\right)
\Bigr\}
\end{align}
```

where ``d(x,y)`` is some pixel-wise distance, ``\MR_{\boldsymbol{θ}}`` is a linear
operator which implements resampling with a given magnification, translation, and
blurring, and ``v_{\mathrm{out}}`` is the assumed out-of-field pixel value.
``\boldsymbol{θ}`` accounts for all parameters defining the operator
``\MR_{\boldsymbol{θ}}``, in particular the translation ``\Vt``, and ``|\Vy|`` denotes the
list of pixels of the image ``\Vy``. Here ``\alpha \in \mathbb{R}``, ``\beta \in
\mathbb{R}``, and the translation ``\Vt \in \mathbb{R}^2`` are nuisance parameters to
reduce the mismatch between the images.

The score may be defined by normalizing the distance:

```math
\Score(\Vx)
= \frac{\Dist(\Vx,\Vy)}{\Dist(v_{\mathrm{out}}\,\One,\Vy)}
```

where ``\One`` is an image of the same size as ``\Vy`` but filled with ones, hence
``v_{\mathrm{out}}\,\One`` is an image of the same size as ``\Vy`` but filled with
``v_{\mathrm{out}}`` the assumed out-of-field pixel value.

The following properties are assumed for the pixel-wise distance:

1. ``d(x,x) = 0`` for any ``x \in \mathbb{R}``;
2. ``d(x,y) > 0`` for any ``(x,y) \in \mathbb{R}^2`` such that ``x \not= y``;
3. ``d(y,x) = d(y,x)`` for any ``(x,y) \in \mathbb{R}^2``.


[^Gomes2016]:
    > N. Gomes, P. J. V. Garcia & É. Thiébaut, *Assessing the quality of
    > restored images in optical long-baseline interferometry* in Monthly
    > Notices of the Royal Astronomical Society, vol. **465**, pp. 3823-3839
    > (2016).

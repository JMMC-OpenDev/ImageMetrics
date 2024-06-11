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
\begin{align*}
\Dist(\Vx,\Vy) = \min_{\Params} \Bigl\{
  &\sum_{j \in |\MR_{\Vtheta}\cdot\Vx| \cap |\Vy|}
  d\!\left(α\,(\MR_{\Vtheta}\cdot\Vx)_j + β, y_j\right)
  \notag\\
  &+ \sum_{j \in |\MR_{\Vtheta}\cdot\Vx| \backslash
  (|\MR_{\Vtheta}\cdot\Vx| \cap |\Vy|)}
  d\!\left(α\,(\MR_{\Vtheta}\cdot\Vx)_j + β, \eta\right)
  \notag\\
  &+ \sum_{j \in |\Vy| \backslash
  (|\MR_{\Vtheta}\cdot\Vx| \cap |\Vy|)}
  d\!\left(\eta,y_j\right)
\Bigr\}
\end{align*}
```

where ``d(x,y)`` is some pixel-wise distance, ``\MR_{\Vtheta}`` is a linear operator which
implements resampling with a given magnification, translation, and blurring, and ``\eta``
is the assumed out-of-field pixel value. ``|\Vy|`` denotes the list of pixels of the image
``\Vy`` and ``\Vtheta = \{\rho,\Vt,\omega\}`` accounts for all parameters defining the
operator ``\MR_{\Vtheta}``: the magnification ``\rho``, the translation ``\Vt``, and the
blur width ``\omega``.

The distance is minimized in the set ``\Params`` of parameters which are irrelevant for
judging of the image quality. These parameters depend on the context. For image
reconstruction from interferometric data, ``\beta = 0`` and ``\eta = 0`` are natural
settings while ``\alpha \in \mathbb{R}`` and the translation ``\Vt \in \mathbb{R}^2`` must
be adjusted to reduce the mismatch between the images. Hence, ``\Params = \{\alpha,\Vt\}``
in this context.

A score may be defined by normalizing the distance and such that the higher the score, the
better the restored image ``\Vx``:

```math
\Score(\Vx) = 1 - \frac{\Dist(\Vx,\Vy)}{\Dist(\eta\,\One,\Vy)}
```

where ``\One`` is an image of the same size as ``\Vy`` but filled with ones, hence
``\eta\,\One`` is an image of the same size as ``\Vy`` but filled with ``\eta``, the
assumed out-of-field pixel value. The score may be negative but the maximal score is 1.

Denoting by ``\mathbb{K}`` the set of possible pixel values, the following properties must
hold for the pixel-wise distance:

1. ``d(x,x) = 0`` for any ``x \in \mathbb{K}``;
2. ``d(x,y) > 0`` for any ``(x,y) \in \mathbb{K}^2`` such that ``x \not= y``;
3. ``d(y,x) = d(y,x)`` for any ``(x,y) \in \mathbb{K}^2``.

To remain general, a possible pixel-wise distance for which the above properties hold is
given by:

```math
d(x, y) = \left|\Gamma(x) - \Gamma(y)\right|^p
```

where the exponent ``p`` and the function ``\Gamma: \mathbb{K}\to\mathbb{R}`` are
introduced to make the distance more flexible. ``\Gamma`` is a _brightness correction_
monotonic function to emphasize the interesting parts of the images amd ``p > 0`` to have
a non-decreasing distance with respect to the absolute value of the difference
``\Gamma(x) - \Gamma(y)``. For example:

``` math
\Gamma(x) = \Sign(x)\,|x|^\gamma,
```

where ``\Sign(x)`` is the sign of ``x``:

``` math
\Sign(x) = \begin{cases}
-1 & \text{if $x < 0$}\\
+1 & \text{if $x > 0$}\\
\phantom{+}0 & \text{if $x = 0$}\\
\end{cases}
```

Metric parameters ``p`` and ``\gamma`` can be chosen depending on the context. According
to a human panel[^Gomes2016], ``p = 1`` with ``\gamma = 1`` best reflect the human
perception of image quality.

[^Gomes2016]:
    > N. Gomes, P. J. V. Garcia & É. Thiébaut, *Assessing the quality of
    > restored images in optical long-baseline interferometry* in Monthly
    > Notices of the Royal Astronomical Society, vol. **465**, pp. 3823-3839
    > (2016).

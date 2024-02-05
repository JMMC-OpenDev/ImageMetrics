# Image metrics

This repository collect notes and tools for image metrics that is means to
quantitatively compare images.

## References

* [Gomnes2016] N. Gomes, P. J. V. Garcia & É. Thiébaut, *Assessing the quality
  of restored images in optical long-baseline interferometry* in Monthly
  Notices of the Royal Astronomical Society, vol. **465**, pp. 3823-3839
  (2016).

* [Sanchez2016] J. Sanchez-Bermudez, É. Thiébaut, K.-H. Hofmann, M. Heininger,
  D. Schertl, G. Weigelt, F. Millour, A. Schutz, A. Ferrari, M. Vannier, D.
  Mary, J. Young & F. Malbet, F., *The 2016 interferometric imaging beauty
  contest* in Optical and Infrared Interferometry and Imaging V, SPIE
  International Conference, **9907**, 99071D (2016).
  (https://doi.org/10.1117/12.2231982)

## Image Metrics Principles

The quality of an image has to be assessed by an objective quantitative
criterion. The literature shows that establishing a universal image quality
criterion is a very controversial subject. What is the best criterion also
largely depends on the context. Here we will assume that the *metric*
$\Theta(x, y)$ is used to estimate the discrepancy between a reconstructed
image $x$ and a reference image $y$. To simplify the discussion, we also assume
that the lower $\Theta(x, y)$ the better the agreement between $x$ and $y$. In
other words $\Theta(x, y)$ can be thought as a measure of the distance between
$x$ and $y$.

When assessing image quality it is important that the result does not depend on
irrelevant changes. This however depends on the type of images and on the
context. For instance, for object detection or recognition, the image metric
should be insensitive to the background level, to a geometrical transform
(translation, rotation, magnification, etc.) or to a multiplication of the
brightness by some positive factor which does not affect the shape of the
object. In cases where image reconstruction has some degeneracies, these should
not have any incidence on the metric. The easier is then to minimize the metric
with respect to the degenerated parameters. These parameters include, but are
not limited to, brightness scaling, geometric transformation of coordinates,
etc. For optical interferometry and when only power-spectrum and closure phase
data are available, the images to be compared may have to be shifted for best
matching.

When comparing a true image (with potentially an infinitely high resolution) to
a restored image, the effective resolution achievable by the instrument and the
image restoration process must be taken into account. Otherwise and because
image metrics are in general based on pixel-wise comparisons, the slightest
displacement of sharp features would lead to large loss of quality (according
to the metric) whereas the images may look very similar at a lower and more
realistic resolution. The easiest solution is then to define the reference
image to be the true image blurred by an effective point spread function (PSF)
whose shape corresponds to the effective resolution. The choice of the
effective resolution is then a parameter of the metric.

To summarize and to be more specific, using the distance $\Theta(x, y)$
between the restored image $x$ and the reference image $y$, the discrepancy
between the restored image $x$ and the true image $z$ would be given by:

``` math
\label{eq:general-discrepancy}
d(x, z) = \min_{\alpha,\beta,t}
\Theta\bigl(\alpha\,x + \beta, h_{\sigma,t} \ast z \bigr),
```

with $\alpha$ a brightness scale, $\beta$ a bias and $h_{\sigma,t}$ an
effective PSF of *width* $\sigma$ and centered at position $t$. The symbol
$\ast$ denotes the convolution. Assuming an effective PSF of Gaussian shape,
the parameter $\sigma > 0$ can be chosen to be the standard deviation of the
PSF; the full width at half maximum ($\FWHM$) is then given by:

``` math
\label{eq:gaussian-fwhm}
\FWHM = \sqrt{8\,\log2}\,\sigma \approx 2.355\,\sigma.
```

Note that the merit function could be minimized with respect to the width of
the PSF in order to estimate the effective resolution achieved by a given
restored image. Our choice to assigning the translation to the PSF is to avoid
relying on some particular method to perform sub-pixel interpolation (of $x$,
$y$ or $z$) for fine tuning the position. Not doing so would add another
ingredient to the metric.

In the following, we first review the most common metrics found in the
literature and argue whether they are appropriate or not in the context of
optical interferometry. We then propose a family of suitable metrics.

**Notations:** In this document, we denote scalars by Greek letters
(e.g. $\alpha$, $\beta$), collection of values, a.k.a. *vectors*, by
lower case Latin letters (e.g. $x$, $y$, $\boldsymbol{z}$), and linear operators, a.k.a.
*matrices*, by upper case Latin letters (e.g., $\mathbf{W}$).

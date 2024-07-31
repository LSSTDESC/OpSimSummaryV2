Notes on noise calculation
==========================

Calculs are inspired from :ref:`R. Biswas et al. 2020 <https://iopscience.iop.org/article/10.3847/1538-4365/ab72f2>` 

The source count in a given band b in unit of ADU is given by:

.. math::
    S_b = \kappa 10^{-0.4 m_b}

The sky count in a given band in unit of ADU is given by:

.. math::
    C_b = \alpha 10^{-0.4 m_b^\mathrm{sky}},

where :math:`m_b^\mathrm{sky}` is the sky magnitudes corresponding to a flux in ADU per pixel.

Thus, the SNR is given as:

.. math::
    \mathrm{SNR} = \frac{S_b}{\left(S_b + C_b\right)^\frac{1}{2}}


The magnitude of a source with an SNR equal to five, :math:`m_b^5`, allows to find the :math:`\kappa` second degree equation:

.. math::
    \kappa^2 - 25 \times 10^{0.8 m_b^5} \kappa - 25 \times  10^{0.4\left(2m_b^5 - m_b^\mathrm{sky}\right)} = 0.

This equation has a unique physical solution :math:`\kappa > 0`:

.. math::
    \kappa = \frac{25}{2}10^{-0.4m_b^5} \left(1 + \sqrt{1 + \frac{4\alpha}{25}10^{-0.4m_b^\mathrm{sky}}}\right).

Using this solution one can compute that:

.. math::
    \kappa = 25 \frac{\alpha}{\kappa} 10^{0.4\left(2m_b^5 - m_b^\mathrm{sky}\right)}\left[1 + \frac{\kappa}{\alpha}10^{-0.4\left(m_b^5 - m_b^\mathrm{sky}\right)}\right].

We will now describes the computation of the ratio :math:`\frac{\kappa}{\alpha}`:

If the source at the top of the atmosphere has a flux density :math:`F_\nu(\lambda)`, its count is given by:

.. math:: 
    S_b = \frac{\pi D^2 T}{4gh}\int_0^\infty F_\nu(\lambda) S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda,

with :math:`D` the mirror diameter, :math:`T` the exposure time, :math:`g` the CCD gain in photo-electron per ADU and :math:`h` the Planck constant.

The magnitude of the source is given by:

.. math:: 

    m_b = -2.5 \log_{10}\left(\frac{\int_0^\infty F_\nu(\lambda) S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda}{\int_0^\infty F^\mathrm{ref}_\nu(\lambda) S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda}\right),

with :math:`F^\mathrm{ref}_\nu(\lambda)` the flux density of a reference source, in the AB magnitude system :math:`F^\mathrm{ref}_\nu(\lambda) = 3631 \ \mathrm{Jy}`. :math:`S^\mathrm{atm}(\lambda)` is the atmosphere transmission 
and :math:`S_b^\mathrm{syst}(\lambda)` the system transmission including optics, filter, ccd efficiency.

The :math:`\kappa` value is thus given by:

.. math::
    \kappa = C \times \int_0^\infty F^\mathrm{ref}_\nu(\lambda) S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda,

with :math:`C = \frac{\pi D^2 T}{4gh}`.

For the sky, the same reasoning gives:

.. math::

    C_b = C  \times  n_\mathrm{eff} \times \int_0^\infty F^\mathrm{ref}_\nu(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda,

where there is no more atmosphere transmission and :math:`n_\mathrm{eff}` is the effective number of pixel that is given by:

.. math::
    n_\mathrm{eff} &= \left(\iint \mathrm{PSF}^2 dS\right)^{-1} \times p_\mathrm{size}^{-2}\\
                   &= 4 \pi \sigma_\mathrm{PSF}^2\\
                   &= \frac{\pi}{\ln2} \mathrm{FWHM}_\mathrm{PSF}^2

with the PSF width given in unit of :math:`\mathrm{arcsec}^{-1}` and the :math:`p_\mathrm{size}` the size of the pixel in :math:`\mathrm{pixel} \ \mathrm{arcsec}^{-1}`.

Thus we have:

.. math::
    \alpha = C \times  n_\mathrm{eff} \times \int_0^\infty F^\mathrm{ref}_\nu(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda.

In the AB magnitude sytem the ratio betweem :math:`\kappa` and :math:`\alpha` is:

.. math::
    \frac{\kappa}{\alpha} = n_\mathrm{eff}^{-1} \frac{\int_0^\infty S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda}{\int_0^\infty S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda}.

Using the approximation  

.. math::
    \frac{\int_0^\infty S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda}{\int_0^\infty S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda} \simeq 1,

we have 

.. math::
    \kappa = 25 n_\mathrm{eff} 10^{0.4\left(2m_b^5 - m_b^\mathrm{sky}\right)}\left[1 + \frac{10^{-0.4\left(m_b^5 - m_b^\mathrm{sky}\right)}}{n_\mathrm{eff}}\right].

We define the zero-point **ZPT** as 

.. math::
    \mathbf{ZPT} = 2.5\log_{10}\left(\kappa\right),

and the source ADU count can be write

.. math::
    S_b = 10^{-0.4(m_b - \mathbf{ZPT})}.

The sky noise by pixel is given by

.. math::
    \sigma_\mathrm{sky}^2 = 10^{-0.4(m_b^\mathrm{sky} - \mathbf{ZPT})} \times p_\mathrm{size}^2.

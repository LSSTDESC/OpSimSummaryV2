Notes on noise calculation
==========================

The computation proposed here follows from `R. Biswas et al. 2020 <https://iopscience.iop.org/article/10.3847/1538-4365/ab72f2>`_. In the following the quantities
in **bold** refers to SNANA quantities. 

The source count in a given band b in unit of ADU is given by

.. math::
    F_b = \kappa 10^{-0.4 m_b},

with :math:`m_b` the magnitude of the source for the corresponding flux in ADU.

The sky count in a given band in unit of ADU is given by

.. math::
    C_b = \alpha 10^{-0.4 m_b^\mathrm{sky}},

where :math:`m_b^\mathrm{sky}` is the sky magnitudes corresponding to a flux in :math:`\mathrm{ADU} \ \mathrm{arsec}^{-2}`.

Thus, the Signal-to-Noise Ratio (SNR) is

.. math::
    \mathrm{SNR} = \frac{F_b}{\left(F_b + C_b\right)^\frac{1}{2}}.

The magnitude of a source with a SNR equal to five, :math:`m_b^5`, allows to find the :math:`\kappa` second degree equation

.. math::
    \kappa^2 - 25 \times 10^{0.8 m_b^5} \kappa - 25 \times  10^{0.4\left(2m_b^5 - m_b^\mathrm{sky}\right)} = 0.

This equation has a unique physical solution :math:`\kappa > 0`,

.. math::
    \kappa = \frac{25}{2}10^{-0.4m_b^5} \left(1 + \sqrt{1 + \frac{4\alpha}{25}10^{-0.4m_b^\mathrm{sky}}}\right).

Using this solution one can compute that:

.. math::
    \kappa = 25 \frac{\alpha}{\kappa} 10^{0.4\left(2m_b^5 - m_b^\mathrm{sky}\right)}\left[1 + \frac{\kappa}{\alpha}10^{-0.4\left(m_b^5 - m_b^\mathrm{sky}\right)}\right].

We will now describes the computation of the ratio :math:`\frac{\kappa}{\alpha}`:

If the source at the top of the atmosphere has a flux density :math:`F_\nu(\lambda)`, its ADU count is given by

.. math:: 
    F_b = \frac{\pi D^2 T}{4gh}\int_0^\infty F_\nu(\lambda) S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda,

with :math:`D` the mirror diameter, :math:`T` the exposure time, :math:`g` the CCD gain in photo-electron per ADU and :math:`h` the Planck constant.

Thus, the magnitude of the source is

.. math:: 

    m_b = -2.5 \log_{10}\left(\frac{\int_0^\infty F_\nu(\lambda) S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda}{\int_0^\infty F^\mathrm{ref}_\nu(\lambda) S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda}\right),

with :math:`F^\mathrm{ref}_\nu(\lambda)` the flux density of a reference source, in the AB magnitude system :math:`F^\mathrm{ref}_\nu(\lambda) = 3631 \ \mathrm{Jy}`. :math:`S^\mathrm{atm}(\lambda)` is the atmosphere transmission 
and :math:`S_b^\mathrm{syst}(\lambda)` the system transmission including optics, filter, ccd efficiency.

We can write

.. math:: 
    F_b =  C \times 10^{-0.4 m_b}\int_0^\infty F^\mathrm{ref}_\nu(\lambda) S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda,

with :math:`C = \frac{\pi D^2 T}{4gh}`.

The :math:`\kappa` value is thus given by

.. math::
    \kappa = C \times \int_0^\infty F^\mathrm{ref}_\nu(\lambda) S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda.

For the sky, the same reasoning gives

.. math::

    C_b = C \times n_\mathrm{eff}p_\mathrm{size}^2 \times 10^{-0.4 m_b^\mathrm{sky}} \int_0^\infty F^\mathrm{ref}_\nu(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda,

where there is no more atmosphere transmission, :math:`p_\mathrm{size}` is the length size of pixels with :math:`p_\mathrm{size}^2` the corresponding pixel area in unit of :math:`\mathrm{arcsec}^2 \ \mathrm{pixel}^{-1}` and, :math:`n_\mathrm{eff}` is the effective number of pixel that is given by

.. math::
    n_\mathrm{eff} &= \left(\iint \mathrm{PSF}^2 dS\right)^{-1} \times p_\mathrm{size}^{-2}\\
                   &= 4 \pi \sigma_\mathrm{PSF}^2 p_\mathrm{size}^{-2}\\
                   &= \frac{\pi}{\ln2} \mathrm{FWHM}_\mathrm{PSF}^2 p_\mathrm{size}^{-2}\\
                   &= A \times p_\mathrm{size}^{-2},


with the PSF width in unit of :math:`\mathrm{arcsec}` and :math:`A` is the noise equivalent area in :math:`\mathrm{arcsec}^{2}`. One can note that the
**PSF** sigma in units of :math:`\mathrm{pixel}` is given by 

.. math::
    \mathbf{PSF} = \frac{\mathrm{FWHM}_\mathrm{PSF}}{2\sqrt{2\ln2}}p_\mathrm{size}^{-1}.

Thus we have

.. math::
    \alpha = C \times A \times \int_0^\infty F^\mathrm{ref}_\nu(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda.

As previously stated, in the AB magnitude sytem the reference flux is constant, :math:`F^\mathrm{ref}_\nu(\lambda) = 3631 \ \mathrm{Jy}` 
and the ratio betweem :math:`\kappa` and :math:`\alpha` is

.. math::
    \frac{\kappa}{\alpha} = A^{-1} \frac{\int_0^\infty S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda}{\int_0^\infty S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda}.

Using the approximation 

.. math::
    \frac{\int_0^\infty S^\mathrm{atm}(\lambda)S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda}{\int_0^\infty S_b^\mathrm{syst}(\lambda)\lambda^{-1}d\lambda} \simeq 1,

we finally obtain 

.. math::
    \kappa = 25 A 10^{0.4\left(2m_b^5 - m_b^\mathrm{sky}\right)}\left[1 + \frac{10^{-0.4\left(m_b^5 - m_b^\mathrm{sky}\right)}}{A}\right].

We define the zero-point **ZPT** as 

.. math::
    \mathbf{ZPT} = 2.5\log_{10}\left(\kappa\right),

such that the source ADU count can be write

.. math::
    F_b = 10^{-0.4(m_b - \mathbf{ZPT})}.

The sky noise **SKYSIG** per pixel is given by

.. math::
    \sigma_\mathrm{sky}^2 = 10^{-0.4(m_b^\mathrm{sky} - \mathbf{ZPT})} \times p_\mathrm{size}^2.

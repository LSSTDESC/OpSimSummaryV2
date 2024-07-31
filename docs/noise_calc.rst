Notes on noise calculation
==========================

Calculs are inspired from :ref:`R. Biswas et al. 2020 <https://iopscience.iop.org/article/10.3847/1538-4365/ab72f2>` 

The source count in a given band b in unit of ADU is given by:

.. math::
    S_b = \kappa 10^{-0.4 m_b}

The sky count in a given band in unit of ADU per pixel is given by:

.. math::
    C_b = \alpha 10^{-0.4 m_b^\mathrm{sky}}

Thus, the SNR is given as:

.. math::
    SNR = \frac{S_b}{\left(S_b + C_b n_\mathrm{eff}\right)}

where :math:`n_\mathrm{eff}` is the effective number of pixel that is given by:

.. math::
    n_\mathrm{eff} &= \left(\int \mathrm{PSF}^2 dS\right)^{-1} \times p_\mathrm{size}^{-2}\\
                   &= 4 \pi \sig_\mathrm{PSF}^2\\
                   &= \frac{\pi}{\ln2} \mathrm{FWHM}_\mathrm{PSF}

with the PSF width given in unit of :math:`\mathrm{arcsec}^{-1}` and the :math:`p_\mathrm{size}` the size of the pixel in :math:`\mathrm{pixel} \mathrm{arcsec}^{-1}`.
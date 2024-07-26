Installation
============

Clone the repository from github and install using pip.

..  code-block:: bash

   git clone https://github.com/bastiencarreres/OpSimSummaryV2.git
   cd OpSimSummaryV2
   pip install .


.. note::
   Some users have encounter this error:

   .. code-block:: bash

      FileNotFoundError: [Errno 2] No such file or directory: 'gdal-config'

   This error is due to a missing dependency and can be resolved by installing the `gdal` library using conda.

   .. code-block:: bash

         conda install -c conda-forge gdal
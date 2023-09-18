#
# log.py
#
# Copyright (c) 2023 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import logging

def set_logging(
        filename='log.txt',
        level=logging.DEBUG,
        format=" %(levelname)8s : %(message)s"
        ):
    
    ### file handler
    fh = logging.FileHandler(filename=filename, mode='w')
    fh.setFormatter(logging.Formatter(format))

    ### stream handler
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter(format))

    logging.basicConfig(level=level, handlers=[fh, sh])

def negative_frequency(fmin):
    
    logger = logging.getLogger(__name__)    
    msg = "\n"
    msg += " Negative eigenvalues were found. Stop the calculation.\n"
    msg += " Minimum frequency : %.2f" % (fmin)
    logger.warning(msg)

def symmetry_error(spg_init, spg_fin):

    logger = logging.getLogger(__name__)
    #msg = "\n Error: The spacegroup number was changed from %d to %d "\
    #        "due to a relaxation calculation.\n" % (spg_init, spg_fin)
    #msg += " The calculation was stopped. Consider to use ISYM = 2."
    msg = "\n Warning: The spacegroup number was changed from %d to %d "\
            "due to a relaxation calculation.\n" % (spg_init, spg_fin)
    msg += " Consider to use ISYM = 2."
    logger.warning(msg)

def rerun_with_omp():
    
    logger = logging.getLogger(__name__)
    msg = "\n Modify the parallel method from MPI to OpenMP and rerun the "\
            "calculation."
    logger.info(msg)


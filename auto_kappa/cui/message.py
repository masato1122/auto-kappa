#
# message.py
#
# This file contains functions to print messages.
#
# Copyright (c) 2023 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import math
import numpy as np

def negative_frequency(fmin):
    msg = "\n"
    msg += " Negative eigenvalues were found. Stop the calculation.\n"
    msg += " Minimum frequency : %.2f" % (fmin)
    print(msg)



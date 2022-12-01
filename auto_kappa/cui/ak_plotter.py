#
# ak_plotter.py
#
# This file is akrun command user interface.
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
from optparse import OptionParser
from auto_kappa.plot.bandos import plot_bandos
import glob

def get_parser():
    
    parser = OptionParser()
    
    parser.add_option("-d", "--directory", dest="directory", type="string",
            default=".", help="directory name [.]")
    
    parser.add_option("--prefix", dest="prefix", type="string",
            default=None, help="prefix")
    
    parser.add_option("--figname", dest="figname", type="string",
            default="fig_bandos.png", help="figure name")
    
    (options, args) = parser.parse_args()
    
    return options

def main():
    
    options = get_parser()
     
    if options.prefix is None:
        line = options.directory + "/*.bands"
        fn = glob.glob(line)[0]
        prefix = fn.split("/")[-1].split('.')[0]
    else:
        prefix = options.prefix

    plot_bandos(
            directory=options.directory,
            prefix=prefix,
            figname=options.figname,
            )

    

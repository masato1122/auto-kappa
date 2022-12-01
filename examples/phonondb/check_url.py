# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser

import urllib.request

def check_url(url):
    flag = True
    try:
        f = urllib.request.urlopen(url)
        #print('OK:', url)
        f.close()
    except urllib.request.HTTPError:
        #print('Not found:', url)
        flag = False
    return flag

def main(options):
    
    flag = check_url(options.url)
    print(flag)
    
if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("--url", dest="url", type="string",
            help="url ")
    
    (options, args) = parser.parse_args()
    main(options)


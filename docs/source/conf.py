import os
import sys
sys.path.insert(0, os.path.abspath('.'))

mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"

# -- Project information 

project = 'aiida-alamode'
copyright = '2021, M. Ohnishi'
author = 'M. Ohnishi'

# The full version, including alpha/beta/rc tags
release = 'Sept. 13th, 2021'


# -- General configuration 

extensions = [
]

templates_path = ['_templates']

exclude_patterns = []


# -- Options for HTML output 

html_theme = 'sphinx_rtd_theme'

html_static_path = ['_static']

## numbering
numfig = True
math_numfig = True
numfig_secnum_depth = 2
math_eqref_format = "Eq.{number}"
html_use_smartypants = False
html_theme = 'sphinx_rtd_theme'

## Latex
latex_engine = 'xelatex'
latex_elements = {
    'fontpkg': r'''
''',
    'preamble': r'''
''',
}
latex_show_urls = 'footnote'


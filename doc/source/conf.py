import os
import sys
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../../auto_kappa'))
#sys.path.insert(0, os.path.abspath('../../auto_kappa/io'))
#sys.path.insert(0, os.path.abspath('../../auto_kappa/plot'))
#sys.path.insert(0, os.path.abspath('../../auto_kappa/plot/alamode'))
#sys.path.insert(0, os.path.abspath('../../auto_kappa/structure'))
#sys.path.insert(0, os.path.abspath('../../auto_kappa/alamode_tools'))

mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"

# -- Project information 

project = 'auto-alamode'
copyright = '2022, M. Ohnishi'
author = 'M. Ohnishi'

# The full version, including alpha/beta/rc tags
release = 'Sept. 13th, 2022'

# -- General configuration 
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary"
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

### open links with a new tab
#window.addEventListener("load", func=()=>{
#    const external_links = document.querySelectorAll("a.external");
#    for(let i=0; i<external_links.length; i++){
#        external_links[i].setAttribute("target", "_blank");
#        external_links[i].setAttribute("rel", "noopener noreferrer");
#    }
#});


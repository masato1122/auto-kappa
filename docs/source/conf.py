
import os, sys
sys.path.insert(0, os.path.abspath('..'))

mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"

# -- Project information 

project = 'auto-kappa'
copyright = '2025, M. Ohnishi'
author = 'M. Ohnishi'

# The full version, including alpha/beta/rc tags
release = 'Sept. 1st, 2025'

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx_rtd_theme",
]
autosummary_generate = True

autosummary_filename_map = {
    "auto_kappa.io": "auto_kappa.io",
    "auto_kappa.plot": "auto_kappa.plot",
}

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "inherited-members": True,
    # "imported-members": True,
}
autoclass_content = "both"
autodoc_typehints = "description"

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output 
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_style = 'css/custom.css'
pygments_style = "sphinx"

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

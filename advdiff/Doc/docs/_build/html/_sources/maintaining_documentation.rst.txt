Maintaining the documentation
===============================

All the documentation is written in `reStructuredText format <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_ 
and then built using `Sphinx <https://www.sphinx-doc.org/en/master/>`_.

To edit and change the documentation, have a look at the .rst files contained within the **docs** directory.
All figures that are shown in the documentation are contained in the **images** directory. Go here if you need to add or change any figures.
When changes are made in the .rst files, you have to save and then re-build the documentation so the changes are reflected in the .html files.
To build the documentation and re-create the html files, have a look at the :ref:`important commands document <Important_Commands>`. 
You should run **clean** before you run **html** such that all changes are taken into account.

The documentation build is contained within the **_build** directory. Here you will find .html files which can be opened in your browser.
The documentation can be hosted on GitHub via `GitHub Pages <https://python.plainenglish.io/how-to-host-your-sphinx-documentation-on-github-550254f325ae>`_, 
however this requires a GitHub Pro account for private repositories, or a public repository. 

Check the `reStructuredText cheat-sheet <https://bashtage.github.io/sphinx-material/rst-cheatsheet/rst-cheatsheet.html>`_ for other helpful commands for writing documentation.

When creating new Python modules, you need to make a corresponding .rst file in the docs folder with the following text.

.. code-block::

    title
    ===================

    .. automodule:: module_name
        :members:
        :undoc-members:
        :show-inheritance:

For new rst files to be visible in the documentation, you need to add the rst file to the 'toctree' via listing it in the index.rst, code.rst or tools.rst file.
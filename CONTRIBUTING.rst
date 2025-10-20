============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

Bug reports
===========

When `reporting a bug <https://github.com/p-j-smith/lipyphilic/issues>`_ please include:

    * Your operating system name and version.
    * Any details about your local setup that might be helpful in troubleshooting.
    * Detailed steps to reproduce the bug.

Documentation improvements
==========================

lipyphilic could always use more documentation, whether as part of the
official `lipyphilic docs <https://lipyphilic.readthedocs.io/en/latest/>`__,
in docstrings, or even on the web in blog posts, articles, and such.

Feature requests and feedback
=============================

The best way to send feedback is to file an issue at https://github.com/p-j-smith/lipyphilic/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that code contributions are welcome :)

Development
===========

To set up `lipyphilic` for local development:

1. Install `Rust <https://www.rust-lang.org/tools/install>`_::

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
    source $HOME/.cargo/env

2. Fork `lipyphilic <https://github.com/p-j-smith/lipyphilic>`_
   (look for the "Fork" button).

3. Clone your fork locally::

    git clone git@github.com:YOURGITHUBNAME/lipyphilic.git
    cd lipyphilic

4. Install uv

    See the `uv installation instructions <https://docs.astral.sh/uv/getting-started/installation/>`_.

5. Create and activate an isolated development environment::

    uv venv --python=3.11
    source venv/bin/activate

6. Install an editible version of `lipyphilic` along with its development dependencies:

    uv sync

7. Create a branch for local development::

    git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

8. You can run the tests using `uv run`::

    uv run pytest --cov

9. Commit your changes and push your branch to GitHub::

    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-bugfix-or-feature

10. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

If you need some code review or feedback while you're developing the code just make the pull request.

For merging, you should:

1. Include passing tests (run ``tox``) [1]_.
2. Update documentation when there's new API, functionality etc.
3. Add a note to ``CHANGELOG.rst`` about the changes.
4. Add yourself to ``AUTHORS.rst``.

Tips
----

To run a subset of tests::

    tox -e envname -- pytest -k test_myfeature

To run all the test environments in *parallel*::

    tox -p auto

To check that the docs build::

    tox -e docs

To run the tests (using python 3.11)::

    tox -e py311

To run tests and print test coverage in the terminal:

    tox -e coverage

To check that the package builds correctly:

    tox -e package

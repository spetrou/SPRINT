======
SPRINT
======


How to add a new RUnit test
---------------------------

*NOTE*: The files in directory ``inst/unitTests`` are deployed in the
        'library' package directory during installation. We place all
        our unit tests in there in order to have them available for
        users to test if needed/wanted.

- If it is a new function create a directory in ``inst/unitTests`` that will
  contain all the unit tests for the function.
- For each unit test create a new file in the directory and name it
  ``runit_<test_name>.R``. The naming convention will help later on with
  selecting different unit tests for different test suites.
  *TODO: find a naming convention for selecting groups of files*


Execute RUnit tests during development
--------------------------------------

For executing unit test during development we create scripts in directory
``test_suite``. See an example of how these files are written. The basic idea
is that they scan the ``inst/unitTests`` directories and select the tests we
want to execute. Each script can execute unit tests for a combination of
functions. Just make sure you name and document them well enough so others
can understand the purpose. Naming conventions are not really needed.



Model and Solve
===============


The Model
---------

.. autoclass:: Numberjack.Model
   :members:
   :special-members: __iadd__


Solving Constructs
------------------


The following is a list of constructs surrounding the solving process.

.. autoclass:: Numberjack.NBJ_STD_Solver
   :members:

.. autoclass:: Numberjack.Solution


.. _findingallsolutions:

Finding all Solutions
~~~~~~~~~~~~~~~~~~~~~

:meth:`Numberjack.NBJ_STD_Solver.startNewSearch` should be called to set the
solver up for finding multiple solutions. Subsequently,
:meth:`Numberjack.NBJ_STD_Solver.getNextSolution` can be called (and checked
against `Numberjack.SAT` to find the next solution. The following simple
example demonstrates the ability to find all solutions:

.. code-block:: python

    xs = VarArray(5, 1, 5)
    model = Model(AllDiff(xs))
    solver = model.load("Mistral")
    solver.startNewSearch()
    while solver.getNextSolution() == SAT:
        print xs

.. note::
    Methods like :meth:`Numberjack.NBJ_STD_Solver.getNodes` return the
    cumulative total between solutions.



Running Numberjack
------------------
.. autofunction:: Numberjack.input

